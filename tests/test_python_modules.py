"""
Tests for Python module components.

Tests the Python replacement modules, XYZ parsing, and integration
with the main sasa_lammps package.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path


class TestSASACoreModule:
    """Test the SASA core Python module."""

    def test_xyz_parsing(self, test_xyz_files, element_radii):
        """Test XYZ file parsing functionality."""
        from sasa_lammps.sasa_core import parse_xyz_file

        # Test single atom
        coords, radii = parse_xyz_file(str(test_xyz_files['single_atom']))
        assert coords.shape == (1, 3), "Should parse single atom coordinates"
        assert radii.shape == (1,), "Should assign single radius"
        assert radii[0] == element_radii['C'], "Should assign carbon radius"
        np.testing.assert_array_almost_equal(coords[0], [0.0, 0.0, 0.0])

        # Test two atoms
        coords, radii = parse_xyz_file(str(test_xyz_files['two_atoms']))
        assert coords.shape == (2, 3), "Should parse two atom coordinates"
        assert radii.shape == (2,), "Should assign two radii"
        expected_coords = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=np.float32)
        np.testing.assert_array_almost_equal(coords, expected_coords)

        # Test mixed elements
        coords, radii = parse_xyz_file(str(test_xyz_files['mixed_elements']))
        assert coords.shape == (4, 3), "Should parse four atoms"
        expected_radii = np.array([
            element_radii['C'], element_radii['H'],
            element_radii['N'], element_radii['O']
        ], dtype=np.float32)
        np.testing.assert_array_almost_equal(radii, expected_radii)

    def test_vdw_radius_assignment(self, element_radii):
        """Test VdW radius assignment for different elements."""
        from sasa_lammps.sasa_core import get_vdw_radius

        # Test known elements
        for element, expected_radius in element_radii.items():
            if element != 'Unknown':  # Skip the default case
                assert get_vdw_radius(element) == expected_radius, \
                    f"Wrong radius for {element}"

        # Test case insensitivity
        assert get_vdw_radius('c') == element_radii['C']
        assert get_vdw_radius('carbon') == element_radii['C']  # Should fall back to 'C'

        # Test unknown elements
        assert get_vdw_radius('Unobtainium') == element_radii['Unknown']
        assert get_vdw_radius('X') == element_radii['Unknown']

        # Test element symbols with numbers (common in PDB files)
        assert get_vdw_radius('C1') == element_radii['C']
        assert get_vdw_radius('N2') == element_radii['N']

    def test_compute_sasa_from_xyz(self, test_xyz_files):
        """Test the complete XYZ to SASA computation pipeline."""
        from sasa_lammps.sasa_core import compute_sasa_from_xyz

        # Test with single atom
        total_sasa, surface_points = compute_sasa_from_xyz(
            str(test_xyz_files['single_atom']),
            srad=1.4, samples=500, points=True
        )

        assert total_sasa > 0, "Should compute positive SASA"
        assert surface_points is not None, "Should return surface points"
        assert len(surface_points) > 0, "Should generate surface points"

        # Test without points
        total_sasa2, surface_points2 = compute_sasa_from_xyz(
            str(test_xyz_files['single_atom']),
            srad=1.4, samples=500, points=False
        )

        assert abs(total_sasa - total_sasa2) < 1e-6, "SASA should be same regardless of points flag"
        assert surface_points2 is None, "Should not return points when points=False"

    def test_create_sasa_xyz_function(self, test_xyz_files, tmp_path):
        """Test the SASA XYZ creation function."""
        from sasa_lammps.sasa_core import _create_sasa_xyz
        from sasa_lammps.constants import SASAXYZ

        # Copy test file to temporary directory
        test_file = test_xyz_files['two_atoms']
        temp_xyz = tmp_path / "test.xyz"
        temp_xyz.write_text(test_file.read_text())

        # Call SASA function
        sasa_points = _create_sasa_xyz(
            str(tmp_path), "test.xyz", srad=1.4, samples=500
        )

        # Check return value
        assert isinstance(sasa_points, np.ndarray), "Should return numpy array"
        assert sasa_points.shape[1] == 3, "Should have 3D coordinates"
        assert len(sasa_points) > 0, "Should generate surface points"

        # Check that SASA file was created
        sasa_file = tmp_path / SASAXYZ
        assert sasa_file.exists(), "Should create SASA XYZ file"

        # Check file format
        content = sasa_file.read_text().strip().split('\n')
        n_points = int(content[0])
        assert n_points == len(sasa_points), "File should contain all surface points"

        # Check that points in file match returned points
        for i, line in enumerate(content[2:]):  # Skip header lines
            parts = line.split()
            assert parts[0] == "He", "Should use He as dummy element"
            file_coords = [float(parts[j]) for j in range(1, 4)]
            np.testing.assert_array_almost_equal(file_coords, sasa_points[i], decimal=5)

    def test_malformed_xyz_handling(self, tmp_path):
        """Test handling of malformed XYZ files."""
        from sasa_lammps.sasa_core import parse_xyz_file

        # Empty file
        empty_file = tmp_path / "empty.xyz"
        empty_file.write_text("")

        with pytest.raises((ValueError, IndexError)):
            parse_xyz_file(str(empty_file))

        # Wrong atom count
        wrong_count_file = tmp_path / "wrong_count.xyz"
        wrong_count_file.write_text("5\nShould have 5 atoms but only has 2\nC 0 0 0\nH 1 0 0\n")

        with pytest.raises((ValueError, IndexError)):
            parse_xyz_file(str(wrong_count_file))

        # Missing coordinates
        missing_coords_file = tmp_path / "missing_coords.xyz"
        missing_coords_file.write_text("1\nMissing coordinates\nC\n")

        with pytest.raises((ValueError, IndexError)):
            parse_xyz_file(str(missing_coords_file))

        # Non-numeric coordinates
        bad_coords_file = tmp_path / "bad_coords.xyz"
        bad_coords_file.write_text("1\nBad coordinates\nC abc def ghi\n")

        with pytest.raises(ValueError):
            parse_xyz_file(str(bad_coords_file))


class TestConversionModuleIntegration:
    """Test integration with the modified conversion module."""

    def test_conversion_module_imports(self):
        """Test that conversion module imports correctly."""
        try:
            from sasa_lammps.conversion import _create_sasa_xyz
            # Should import without error
            assert callable(_create_sasa_xyz), "Should import SASA function"
        except ImportError as e:
            # May fail if dependencies like ovito or sasa_ext are missing - that's expected
            assert any(dep in str(e).lower() for dep in ["ovito", "sasa_ext", "sasa_core"]), \
                f"Unexpected import error: {e}"

    @pytest.mark.skipif(True, reason="Requires full package installation with dependencies")
    def test_create_sasa_xyz_integration(self, test_xyz_files, tmp_path):
        """Test the integrated _create_sasa_xyz function."""
        from sasa_lammps.conversion import _create_sasa_xyz

        # Copy test file to temporary directory
        test_file = test_xyz_files['single_atom']
        temp_xyz = tmp_path / "test.xyz"
        temp_xyz.write_text(test_file.read_text())

        # Call integrated function
        sasa_points = _create_sasa_xyz(str(tmp_path), "test.xyz", srad=1.4, samples=500)

        # Should return valid results
        assert isinstance(sasa_points, np.ndarray)
        assert sasa_points.shape[1] == 3
        assert len(sasa_points) > 0


class TestPackageStructure:
    """Test overall package structure and imports."""

    def test_sasa_ext_import_structure(self):
        """Test SASA extension import structure."""
        try:
            import sasa_ext

            # Should have main function
            assert hasattr(sasa_ext, 'compute_sasa'), "Should have compute_sasa function"
            assert callable(sasa_ext.compute_sasa), "compute_sasa should be callable"

            # Test function signature (C extensions may not have introspectable signatures)
            import inspect
            try:
                sig = inspect.signature(sasa_ext.compute_sasa)
                param_names = list(sig.parameters.keys())

                expected_params = ['coords', 'radii']  # At minimum
                for param in expected_params:
                    assert param in param_names, f"Missing parameter: {param}"
            except ValueError:
                # C extensions often don't have introspectable signatures
                # Just verify we can call it with expected arguments
                import numpy as np
                test_coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
                test_radii = np.array([1.5], dtype=np.float32)

                # This should work without error
                result = sasa_ext.compute_sasa(test_coords, test_radii)
                assert len(result) == 2, "Should return (sasa, points) tuple"
                assert isinstance(result[0], float), "First result should be float"
                assert hasattr(result[1], '__len__'), "Second result should be array-like"

        except ImportError:
            pytest.skip("SASA extension not available")

    def test_sasa_core_import_structure(self):
        """Test SASA core module structure."""
        from sasa_lammps import sasa_core

        # Should have key functions
        expected_functions = [
            'parse_xyz_file',
            'get_vdw_radius',
            'compute_sasa_from_xyz',
            '_create_sasa_xyz'
        ]

        for func_name in expected_functions:
            assert hasattr(sasa_core, func_name), f"Missing function: {func_name}"
            assert callable(getattr(sasa_core, func_name)), f"{func_name} not callable"

    def test_constants_availability(self):
        """Test that required constants are available."""
        from sasa_lammps.constants import SASAXYZ

        assert isinstance(SASAXYZ, str), "SASAXYZ should be a string"
        assert SASAXYZ.endswith('.xyz'), "SASAXYZ should be an XYZ filename"

    def test_package_version_consistency(self):
        """Test that package version is consistently defined."""
        # Check setup.py for version
        setup_py_path = Path(__file__).parent.parent / "setup.py"
        if setup_py_path.exists():
            setup_content = setup_py_path.read_text()
            assert 'version=' in setup_content, "setup.py should have version"

        # Check pyproject.toml for version
        pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
        if pyproject_path.exists():
            pyproject_content = pyproject_path.read_text()
            assert 'version =' in pyproject_content, "pyproject.toml should have version"


class TestErrorPropagation:
    """Test that errors are properly propagated through the Python layers."""

    def test_c_extension_error_propagation(self):
        """Test that C extension errors are properly caught and re-raised."""
        import sasa_ext

        # Test invalid input that should cause C extension to fail
        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(
                np.array([]), np.array([1.0])  # Mismatched shapes
            )

    def test_xyz_parsing_error_propagation(self, tmp_path):
        """Test that XYZ parsing errors are properly propagated."""
        from sasa_lammps.sasa_core import parse_xyz_file

        # Non-existent file
        with pytest.raises((FileNotFoundError, IOError)):
            parse_xyz_file("nonexistent_file.xyz")

        # Malformed file
        bad_file = tmp_path / "bad.xyz"
        bad_file.write_text("not a valid xyz file")

        with pytest.raises((ValueError, IndexError)):
            parse_xyz_file(str(bad_file))

    def test_sasa_computation_error_propagation(self, tmp_path):
        """Test error propagation in SASA computation pipeline."""
        from sasa_lammps.sasa_core import compute_sasa_from_xyz

        # Non-existent file
        with pytest.raises((FileNotFoundError, IOError)):
            compute_sasa_from_xyz("nonexistent.xyz")

        # Invalid parameters should be caught
        valid_file = tmp_path / "valid.xyz"
        valid_file.write_text("1\nTest\nC 0 0 0\n")

        with pytest.raises((ValueError, RuntimeError)):
            compute_sasa_from_xyz(str(valid_file), srad=1.4, samples=-1)  # Negative samples