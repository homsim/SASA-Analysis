"""
SASA computation module using efficient C extension.

This module provides solvent accessible surface area calculations
using Monte Carlo sampling on atomic surfaces.
"""

import numpy as np
from pathlib import Path
import sasa_ext

from sasa_lammps.constants import RADII_MAP, SAS_SEED

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
    # Clean element symbol (remove digits, make uppercase)
    clean_element = ''.join(c for c in element if c.isalpha()).capitalize()

    return RADII_MAP.get(clean_element, 1.70)  # Default to carbon radius

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
        seed=SAS_SEED
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

    xyz_file_path = Path(path) / xyz_file

    # Compute SASA surface points
    _, sasa_points = compute_sasa_from_xyz(xyz_file_path, srad=srad, samples=samples, points=True)

    # Convert to numpy array
    sasa_points = np.array(sasa_points)

    # Export points in expected format
    export_points = np.array(sasa_points, dtype=str)
    export_points = np.insert(export_points, 0, "He", axis=1)

    header = f"{len(export_points)}\n "
    np.savetxt(
        Path(path) / SASAXYZ,
        export_points,
        header=header,
        comments="",
        fmt="%s",
    )

    return sasa_points