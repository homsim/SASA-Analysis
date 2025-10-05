"""
SASA computation module using efficient C extension.

This module provides solvent accessible surface area calculations
using Monte Carlo sampling on atomic surfaces.
"""

import os
import numpy as np
try:
    import sasa_ext
    SASA_EXT_AVAILABLE = True
except ImportError:
    SASA_EXT_AVAILABLE = False

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

    Uses standard VdW radii from literature.
    """
    # VdW radii in Angstroms - standard values
    radii_map = {
        'H': 1.20,   # Hydrogen
        'C': 1.70,   # Carbon
        'N': 1.55,   # Nitrogen
        'O': 1.52,   # Oxygen
        'S': 1.80,   # Sulfur
        'P': 1.80,   # Phosphorus
        'F': 1.47,   # Fluorine
        'Cl': 1.75,  # Chlorine
        'Br': 1.85,  # Bromine
        'I': 1.98,   # Iodine
        'Na': 2.27,  # Sodium
        'K': 2.75,   # Potassium
        'Ca': 2.31,  # Calcium
        'Mg': 1.73,  # Magnesium
        'Fe': 2.04,  # Iron
        'Zn': 1.39,  # Zinc
        'He': 1.40,  # Helium (used as dummy in some test files)
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
    if not SASA_EXT_AVAILABLE:
        raise ImportError("SASA extension not available. Install package with C extension support.")

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