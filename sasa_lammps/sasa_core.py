"""
SASA computation module using efficient C extension.

This module provides solvent accessible surface area calculations
using Monte Carlo sampling on atomic surfaces.
"""

import os
import numpy as np
import sasa_ext

def parse_xyz_file(xyz_file_path):
    """
    Parse XYZ file to extract coordinates and atom types.

    Parameters
    ----------
    xyz_file_path : str
        Path to the XYZ file

    Returns
    -------
    coords : numpy.ndarray
        (N, 3) array of atomic coordinates in Angstroms
    radii : numpy.ndarray
        (N,) array of VdW radii for each atom
    """
    with open(xyz_file_path, 'r') as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("Empty XYZ file")

    try:
        n_atoms = int(lines[0].strip())
    except (ValueError, IndexError) as e:
        raise ValueError(f"Could not read atom count from first line: {e}")

    if len(lines) < 2:
        raise ValueError("XYZ file must have at least 2 lines (atom count + comment)")

    coords = []
    elements = []

    for i in range(2, 2 + n_atoms):  # Skip header lines
        if i >= len(lines):
            raise ValueError(f"XYZ file claims {n_atoms} atoms but only has {i-2} atom lines")

        parts = lines[i].strip().split()
        if len(parts) < 4:
            raise ValueError(f"Line {i+1}: Expected at least 4 fields (element x y z), got {len(parts)}")

        element = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError as e:
            raise ValueError(f"Line {i+1}: Could not parse coordinates: {e}")

        elements.append(element)
        coords.append([x, y, z])

    coords = np.array(coords, dtype=np.float32)
    radii = np.array([get_vdw_radius(elem) for elem in elements], dtype=np.float32)

    return coords, radii

def get_vdw_radius(element):
    """
    Get VdW radius for element (in Angstroms).

    Uses standard VdW radii from authoritative literature sources.
    Values are hierarchically selected from:
    1. Bondi, A. (1964). "van der Waals Volumes and Radii".
       J. Phys. Chem. 68(3), 441-451. DOI: 10.1021/j100785a001
    2. Batsanov, S.S. (2001). "Van der Waals radii of elements".
       Inorg. Mater. 37(9), 871-885.
    3. Alvarez, S. (2013). "A cartography of the van der Waals territories".
       Dalton Trans. 42(24), 8617-8636. DOI: 10.1039/C3DT50599E

    Priority: Bondi > Batsanov > Alvarez for maximum compatibility.

    Parameters
    ----------
    element : str
        Element symbol (case-insensitive, digits removed)

    Returns
    -------
    float
        Van der Waals radius in Angstroms
    """
    # VdW radii in Angstroms - hierarchical selection: Bondi > Batsanov > Alvarez
    radii_map = {
        # Period 1
        'H': 1.20,   # Hydrogen (Bondi)
        'He': 1.40,  # Helium (Bondi)

        # Period 2
        'Li': 1.81,  # Lithium (Bondi)
        'Be': 1.9,   # Beryllium (Batsanov)
        'B': 1.8,    # Boron (Batsanov)
        'C': 1.70,   # Carbon (Bondi)
        'N': 1.55,   # Nitrogen (Bondi)
        'O': 1.52,   # Oxygen (Bondi)
        'F': 1.47,   # Fluorine (Bondi)
        'Ne': 1.54,  # Neon (Bondi)

        # Period 3
        'Na': 2.27,  # Sodium (Bondi)
        'Mg': 1.73,  # Magnesium (Bondi)
        'Al': 2.1,   # Aluminum (Batsanov)
        'Si': 2.10,  # Silicon (Batsanov)
        'P': 1.80,   # Phosphorus (Bondi)
        'S': 1.80,   # Sulfur (Bondi)
        'Cl': 1.75,  # Chlorine (Bondi)
        'Ar': 1.88,  # Argon (Bondi)

        # Period 4
        'K': 2.75,   # Potassium (Bondi)
        'Ca': 2.4,   # Calcium (Batsanov)
        'Sc': 2.3,   # Scandium (Batsanov)
        'Ti': 2.15,  # Titanium (Batsanov)
        'V': 2.05,   # Vanadium (Batsanov)
        'Cr': 2.05,  # Chromium (Batsanov)
        'Mn': 2.05,  # Manganese (Batsanov)
        'Fe': 2.05,  # Iron (Batsanov)
        'Co': 2.0,   # Cobalt (Batsanov)
        'Ni': 1.63,  # Nickel (Bondi)
        'Cu': 1.40,  # Copper (Bondi)
        'Zn': 1.39,  # Zinc (Bondi)
        'Ga': 1.87,  # Gallium (Bondi)
        'Ge': 2.1,   # Germanium (Batsanov)
        'As': 1.85,  # Arsenic (Bondi)
        'Se': 1.90,  # Selenium (Bondi)
        'Br': 1.85,  # Bromine (Bondi)
        'Kr': 2.02,  # Krypton (Bondi)

        # Period 5
        'Rb': 2.9,   # Rubidium (Batsanov)
        'Sr': 2.55,  # Strontium (Batsanov)
        'Y': 2.4,    # Yttrium (Batsanov)
        'Zr': 2.3,   # Zirconium (Batsanov)
        'Nb': 2.15,  # Niobium (Batsanov)
        'Mo': 2.1,   # Molybdenum (Batsanov)
        'Tc': 2.05,  # Technetium (Batsanov)
        'Ru': 2.05,  # Ruthenium (Batsanov)
        'Rh': 2.0,   # Rhodium (Batsanov)
        'Pd': 1.63,  # Palladium (Bondi)
        'Ag': 1.72,  # Silver (Bondi)
        'Cd': 1.58,  # Cadmium (Bondi)
        'In': 1.93,  # Indium (Bondi)
        'Sn': 2.17,  # Tin (Bondi)
        'Sb': 2.2,   # Antimony (Batsanov)
        'Te': 2.06,  # Tellurium (Bondi)
        'I': 1.98,   # Iodine (Bondi)
        'Xe': 2.16,  # Xenon (Bondi)

        # Period 6
        'Cs': 3.0,   # Cesium (Batsanov)
        'Ba': 2.7,   # Barium (Batsanov)
        'La': 2.5,   # Lanthanum (Batsanov)
        'Ce': 2.88,  # Cerium (Alvarez)
        'Pr': 2.92,  # Praseodymium (Alvarez)
        'Nd': 2.95,  # Neodymium (Alvarez)
        'Pm': 2.90,  # Promethium (Alvarez)
        'Sm': 2.90,  # Samarium (Alvarez)
        'Eu': 2.87,  # Europium (Alvarez)
        'Gd': 2.83,  # Gadolinium (Alvarez)
        'Tb': 2.79,  # Terbium (Alvarez)
        'Dy': 2.87,  # Dysprosium (Alvarez)
        'Ho': 2.81,  # Holmium (Alvarez)
        'Er': 2.83,  # Erbium (Alvarez)
        'Tm': 2.79,  # Thulium (Alvarez)
        'Yb': 2.80,  # Ytterbium (Alvarez)
        'Lu': 2.74,  # Lutetium (Alvarez)
        'Hf': 2.25,  # Hafnium (Batsanov)
        'Ta': 2.2,   # Tantalum (Batsanov)
        'W': 2.1,    # Tungsten (Batsanov)
        'Re': 2.05,  # Rhenium (Batsanov)
        'Os': 2.0,   # Osmium (Batsanov)
        'Ir': 2.0,   # Iridium (Batsanov)
        'Pt': 1.72,  # Platinum (Bondi)
        'Au': 1.66,  # Gold (Bondi)
        'Hg': 1.55,  # Mercury (Bondi)
        'Tl': 1.96,  # Thallium (Bondi)
        'Pb': 2.02,  # Lead (Bondi)
        'Bi': 2.3,   # Bismuth (Batsanov)
        'Po': 1.97,  # Polonium (Alvarez)
        'At': 2.02,  # Astatine (Alvarez)
        'Rn': 2.20,  # Radon (Alvarez)

        # Period 7 - Actinides (mostly Alvarez, Th from Batsanov)
        'Fr': 3.48,  # Francium (Alvarez)
        'Ra': 2.83,  # Radium (Alvarez)
        'Ac': 2.8,   # Actinium (Alvarez)
        'Th': 2.4,   # Thorium (Batsanov)
        'Pa': 2.88,  # Protactinium (Alvarez)
        'U': 1.86,   # Uranium (Bondi)
        'Np': 2.82,  # Neptunium (Alvarez)
        'Pu': 2.81,  # Plutonium (Alvarez)
        'Am': 2.83,  # Americium (Alvarez)
        'Cm': 3.05,  # Curium (Alvarez)
        'Bk': 3.4,   # Berkelium (Alvarez)
        'Cf': 3.05,  # Californium (Alvarez)
        'Es': 2.7,   # Einsteinium (Alvarez)

        # Superheavy elements (estimated)
        'Fm': 2.45,  # Fermium (estimated)
        'Md': 2.46,  # Mendelevium (estimated)
        'No': 2.46,  # Nobelium (estimated)
        'Lr': 2.46,  # Lawrencium (estimated)
        'Rf': 2.30,  # Rutherfordium (estimated)
        'Db': 2.30,  # Dubnium (estimated)
        'Sg': 2.30,  # Seaborgium (estimated)
        'Bh': 2.30,  # Bohrium (estimated)
        'Hs': 2.30,  # Hassium (estimated)
        'Mt': 2.30,  # Meitnerium (estimated)
        'Ds': 2.30,  # Darmstadtium (estimated)
        'Rg': 2.30,  # Roentgenium (estimated)
        'Cn': 2.30,  # Copernicium (estimated)
        'Nh': 2.30,  # Nihonium (estimated)
        'Fl': 2.30,  # Flerovium (estimated)
        'Mc': 2.30,  # Moscovium (estimated)
        'Lv': 2.30,  # Livermorium (estimated)
        'Ts': 2.30,  # Tennessine (estimated)
        'Og': 2.30,  # Oganesson (estimated)
    }

    # Clean element symbol (remove digits, make uppercase)
    clean_element = ''.join(c for c in element if c.isalpha()).capitalize()

    return radii_map.get(clean_element, 1.70)  # Default to carbon radius

def compute_sasa_from_xyz(xyz_file_path, srad=1.4, samples=500, points=True):
    """
    Compute SASA from XYZ file using Monte Carlo sampling.

    Parameters
    ----------
    xyz_file_path : str
        Path to XYZ coordinate file
    srad : float
        Probe radius (solvent radius)
    samples : int
        Number of Monte Carlo sample points per atom
    points : bool
        Whether to return surface point coordinates

    Returns
    -------
    total_sasa : float
        Total solvent accessible surface area
    surface_points : numpy.ndarray or None
        (N, 3) array of surface point coordinates if points=True, else None
    """

    # Parse coordinates and radii from XYZ file
    coords, radii = parse_xyz_file(xyz_file_path)

    # Compute SASA using C extension
    # Use fixed seed for reproducible results
    total_sasa, surface_points = sasa_ext.compute_sasa(
        coords, radii,
        probe_radius=srad,
        n_samples=samples,
        seed=38572111  # Fixed seed for reproducibility
    )

    if points:
        return total_sasa, surface_points.tolist()  # Convert to list for compatibility
    else:
        return total_sasa, None

def _create_sasa_xyz(path, xyz_file, srad, samples):
    """
    Create van der Waals surface points for molecular analysis.

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    xyz_file : str
        Name of the xyz-file to use
    srad : float
        Probe radius: Effectively a scaling factor for the vdW radii
    samples : int
        Maximum points on the atomic vdW sphere to generate per atom

    Returns
    -------
    sasa_points : numpy.ndarray
        (N, 3) Array of coordinates on the SAS.
        N is loosely determined by 'samples' argument.
    """
    from sasa_lammps.constants import SASAXYZ

    xyz_file_path = os.path.join(path, xyz_file)

    # Compute SASA surface points
    _, sasa_points = compute_sasa_from_xyz(xyz_file_path, srad=srad, samples=samples, points=True)

    # Convert to numpy array
    sasa_points = np.array(sasa_points)

    # Export points in expected format
    export_points = np.array(sasa_points, dtype=str)
    export_points = np.insert(export_points, 0, "He", axis=1)

    header = f"{len(export_points)}\n "
    np.savetxt(
        os.path.join(path, SASAXYZ),
        export_points,
        header=header,
        comments="",
        fmt="%s",
    )

    return sasa_points