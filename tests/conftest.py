"""
Pytest configuration and shared fixtures for SASA tests.
"""

import pytest
import numpy as np
import json
import tempfile
import os
from pathlib import Path

import sasa_ext


@pytest.fixture
def reference_data():
    """Load VMD reference data if available."""
    ref_file = Path(__file__).parent / "resources" / "vmd_reference_data.json"
    if ref_file.exists():
        with open(ref_file, 'r') as f:
            return json.load(f)
    return {}

@pytest.fixture
def simple_test_cases():
    """Create simple test cases with known geometry."""
    return {
        'single_atom': {
            'coords': np.array([[0.0, 0.0, 0.0]], dtype=np.float32),
            'radii': np.array([1.5], dtype=np.float32),
            'expected_area_factor': 4 * np.pi  # For probe_radius=0
        },
        'two_separated': {
            'coords': np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float32),
            'radii': np.array([1.5, 1.5], dtype=np.float32),
            'probe_radius': 1.4,
            'expected_behavior': 'independent'  # Should act like two isolated atoms
        },
        'two_touching': {
            'coords': np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32),
            'radii': np.array([1.5, 1.5], dtype=np.float32),
            'probe_radius': 1.4,
            'expected_behavior': 'reduced'  # SASA should be less than 2x single atom
        },
        'buried_atom': {
            'coords': np.array([
                [0.0, 0.0, 0.0],   # Central atom
                [2.0, 0.0, 0.0],   # Surrounding atoms
                [-2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, -2.0, 0.0],
                [0.0, 0.0, 2.0],
                [0.0, 0.0, -2.0]
            ], dtype=np.float32),
            'radii': np.array([1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], dtype=np.float32),
            'probe_radius': 1.4,
            'expected_behavior': 'central_buried'  # Central atom should be mostly buried
        }
    }

@pytest.fixture
def test_xyz_files(tmp_path):
    """Create test XYZ files for parsing tests."""
    files = {}

    # Single atom
    single_atom = tmp_path / "single_atom.xyz"
    single_atom.write_text("1\nSingle carbon atom\nC 0.0 0.0 0.0\n")
    files['single_atom'] = single_atom

    # Two atoms
    two_atoms = tmp_path / "two_atoms.xyz"
    two_atoms.write_text("2\nTwo carbon atoms\nC 0.0 0.0 0.0\nC 5.0 0.0 0.0\n")
    files['two_atoms'] = two_atoms

    # Mixed elements
    mixed = tmp_path / "mixed_elements.xyz"
    mixed.write_text("""4
Mixed element molecule
C 0.0 0.0 0.0
H 1.0 0.0 0.0
N 0.0 1.0 0.0
O 0.0 0.0 1.0
""")
    files['mixed_elements'] = mixed

    # Large molecule (for performance tests)
    large = tmp_path / "large_molecule.xyz"
    n_atoms = 100
    content = f"{n_atoms}\nLarge test molecule\n"
    for i in range(n_atoms):
        # Random coordinates in a 20x20x20 box
        x, y, z = np.random.rand(3) * 20 - 10
        element = np.random.choice(['C', 'N', 'O', 'H'])
        content += f"{element} {x:.6f} {y:.6f} {z:.6f}\n"
    large.write_text(content)
    files['large_molecule'] = large

    return files

@pytest.fixture
def element_radii():
    """Expected VdW radii for testing radius assignment."""
    return {
        'H': 1.20,
        'C': 1.70,
        'N': 1.55,
        'O': 1.52,
        'S': 1.80,
        'P': 1.80,
        'F': 1.47,
        'Cl': 1.75,
        'Unknown': 1.70  # Default
    }

@pytest.fixture
def sasa_parameters():
    """Standard SASA parameter sets for testing."""
    return {
        'standard': {'probe_radius': 1.4, 'n_samples': 500, 'seed': 38572111},
        'small_probe': {'probe_radius': 1.0, 'n_samples': 500, 'seed': 38572111},
        'large_probe': {'probe_radius': 2.0, 'n_samples': 500, 'seed': 38572111},
        'few_samples': {'probe_radius': 1.4, 'n_samples': 100, 'seed': 38572111},
        'many_samples': {'probe_radius': 1.4, 'n_samples': 1000, 'seed': 38572111},
        'zero_probe': {'probe_radius': 0.0, 'n_samples': 500, 'seed': 38572111},
    }

@pytest.fixture
def convergence_tolerances():
    """Tolerance values for different types of comparisons."""
    return {
        'area_relative': 0.005,      # 0.5% for total SASA area
        'area_absolute': 1.0,        # 1.0 absolute tolerance
        'point_count_relative': 0.05, # 5% for point counts
        'reproducibility': 1e-10,    # Exact for same seed
        'convergence': 0.02,         # 2% for sample convergence
        'geometric': 0.01            # 1% for geometric relationships
    }

def create_test_molecule_xyz(coords, elements, filename):
    """Helper function to create XYZ files from coordinates and elements."""
    content = f"{len(coords)}\nTest molecule\n"
    for elem, coord in zip(elements, coords):
        content += f"{elem} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n"

    with open(filename, 'w') as f:
        f.write(content)

    return filename

def load_xyz_file(filepath):
    """Helper function to load coordinates and elements from XYZ file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    coords = []
    elements = []

    for i in range(2, 2 + n_atoms):
        parts = lines[i].strip().split()
        elements.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(coords, dtype=np.float32), elements

@pytest.fixture
def xyz_loader():
    """Provide the XYZ loading function as a fixture."""
    return load_xyz_file
