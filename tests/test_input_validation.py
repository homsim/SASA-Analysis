"""
Input validation and error handling tests.

Tests edge cases, error conditions, and input validation for the SASA extension.
"""

import pytest
import numpy as np


class TestInputValidation:
    """Test input validation and error handling."""

    def test_array_shape_validation(self, sasa_ext_available):
        """Test validation of input array shapes."""
        import sasa_ext

        valid_coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        valid_radii = np.array([1.5], dtype=np.float32)

        # Test mismatched array sizes
        wrong_radii = np.array([1.5, 2.0], dtype=np.float32)
        with pytest.raises(ValueError, match="same number of atoms"):
            sasa_ext.compute_sasa(valid_coords, wrong_radii)

        # Test wrong coordinate shape (2D instead of 3D)
        wrong_coords_2d = np.array([[0.0, 0.0]], dtype=np.float32)
        with pytest.raises(ValueError, match="shape \\(n_atoms, 3\\)"):
            sasa_ext.compute_sasa(wrong_coords_2d, valid_radii)

        # Test wrong coordinate shape (1D)
        wrong_coords_1d = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        with pytest.raises(ValueError, match="(shape \\(n_atoms, 3\\)|too small depth|shape|2D)"):
            sasa_ext.compute_sasa(wrong_coords_1d, valid_radii)

        # Test wrong radii shape (2D instead of 1D)
        wrong_radii_2d = np.array([[1.5]], dtype=np.float32)
        with pytest.raises(ValueError, match="(1D array|too deep|shape)"):
            sasa_ext.compute_sasa(valid_coords, wrong_radii_2d)

    def test_array_type_conversion(self, sasa_ext_available):
        """Test that arrays are properly converted to correct types."""
        import sasa_ext

        # Test with different dtypes - should be automatically converted
        coords_int = np.array([[0, 0, 0]], dtype=np.int32)
        radii_float64 = np.array([1.5], dtype=np.float64)

        # Should work without error
        sasa, points = sasa_ext.compute_sasa(coords_int, radii_float64)
        assert sasa > 0
        assert len(points) > 0

        # Test with Python lists - should be converted
        coords_list = [[0.0, 0.0, 0.0]]
        radii_list = [1.5]

        sasa, points = sasa_ext.compute_sasa(coords_list, radii_list)
        assert sasa > 0
        assert len(points) > 0

    def test_empty_arrays(self, sasa_ext_available):
        """Test behavior with empty arrays."""
        import sasa_ext

        empty_coords = np.array([], dtype=np.float32).reshape(0, 3)
        empty_radii = np.array([], dtype=np.float32)

        # Should handle empty arrays gracefully
        sasa, points = sasa_ext.compute_sasa(empty_coords, empty_radii)
        assert sasa == 0.0
        assert len(points) == 0

    def test_parameter_validation(self, sasa_ext_available):
        """Test validation of algorithm parameters."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        # Test negative probe radius - should work (physical meaning: atom shrinkage)
        sasa, points = sasa_ext.compute_sasa(coords, radii, probe_radius=-0.5)
        assert sasa >= 0  # SASA should still be non-negative

        # Test zero samples - should fail or handle gracefully
        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(coords, radii, n_samples=0)

        # Test negative samples - should fail
        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(coords, radii, n_samples=-100)

    def test_extreme_coordinates(self, sasa_ext_available):
        """Test behavior with extreme coordinate values."""
        import sasa_ext

        # Very large coordinates
        large_coords = np.array([[1e6, 1e6, 1e6]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        sasa, points = sasa_ext.compute_sasa(large_coords, radii)
        assert sasa > 0
        assert len(points) > 0

        # Very small coordinates
        small_coords = np.array([[1e-6, 1e-6, 1e-6]], dtype=np.float32)
        sasa, points = sasa_ext.compute_sasa(small_coords, radii)
        assert sasa > 0
        assert len(points) > 0

    def test_extreme_radii(self, sasa_ext_available):
        """Test behavior with extreme radii values."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)

        # Very small radius
        tiny_radii = np.array([1e-6], dtype=np.float32)
        sasa, points = sasa_ext.compute_sasa(coords, tiny_radii, n_samples=100)
        assert sasa >= 0
        # Very small atoms might have no exposed points depending on probe radius

        # Very large radius
        huge_radii = np.array([1e3], dtype=np.float32)
        sasa, points = sasa_ext.compute_sasa(coords, huge_radii, n_samples=100)
        assert sasa > 0
        assert len(points) > 0

        # Zero radius
        zero_radii = np.array([0.0], dtype=np.float32)
        sasa, points = sasa_ext.compute_sasa(coords, zero_radii)
        assert sasa >= 0  # Should handle gracefully

    def test_nan_and_inf_handling(self, sasa_ext_available):
        """Test handling of NaN and infinite values."""
        import sasa_ext

        # NaN coordinates
        nan_coords = np.array([[np.nan, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(nan_coords, radii)

        # Infinite coordinates
        inf_coords = np.array([[np.inf, 0.0, 0.0]], dtype=np.float32)
        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(inf_coords, radii)

        # NaN radii
        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        nan_radii = np.array([np.nan], dtype=np.float32)
        with pytest.raises((ValueError, RuntimeError)):
            sasa_ext.compute_sasa(coords, nan_radii)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_point_sampling(self, sasa_ext_available):
        """Test with minimal sampling (n_samples=1)."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        # Should still run, though inaccurate
        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=1)
        assert sasa >= 0
        assert len(points) <= 1  # At most 1 point

    def test_identical_atom_positions(self, sasa_ext_available):
        """Test behavior with atoms at identical positions."""
        import sasa_ext

        # Two atoms at same position
        coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # Should handle without crashing
        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=100)
        assert sasa >= 0

    def test_very_close_atoms(self, sasa_ext_available):
        """Test behavior with atoms very close together."""
        import sasa_ext

        # Atoms almost overlapping
        coords = np.array([[0.0, 0.0, 0.0], [1e-6, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=500)
        assert sasa >= 0

        # Should be similar to single atom case due to overlap
        single_sasa, _ = sasa_ext.compute_sasa(coords[:1], radii[:1], n_samples=500)
        # Allow for significant variation due to overlap effects
        assert abs(sasa - single_sasa) / single_sasa < 1.0  # Within 100%

    def test_linear_configuration(self, sasa_ext_available):
        """Test with atoms in perfect linear arrangement."""
        import sasa_ext

        # Five atoms in a line
        coords = np.array([
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [9.0, 0.0, 0.0],
            [12.0, 0.0, 0.0]
        ], dtype=np.float32)
        radii = np.array([1.5, 1.5, 1.5, 1.5, 1.5], dtype=np.float32)

        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=500, seed=42)
        assert sasa > 0
        assert len(points) > 0

        # Middle atoms should have less SASA than end atoms
        # Test by computing individual atom contributions (approximate)
        end_atom_sasa, _ = sasa_ext.compute_sasa(coords[:1], radii[:1], n_samples=500, seed=42)

        # Total should be less than 5x single atom due to neighbor effects
        assert sasa < 4.0 * end_atom_sasa, "Linear arrangement should show neighbor effects"

    def test_planar_configuration(self, sasa_ext_available):
        """Test with atoms in a plane."""
        import sasa_ext

        # Atoms arranged in a plane
        coords = np.array([
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [3.0, 3.0, 0.0]
        ], dtype=np.float32)
        radii = np.array([1.5, 1.5, 1.5, 1.5], dtype=np.float32)

        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=500, seed=42)
        assert sasa > 0
        assert len(points) > 0

    def test_very_large_probe(self, sasa_ext_available):
        """Test with probe radius much larger than atoms."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.0], dtype=np.float32)

        # Probe radius 10x larger than atom
        huge_probe = 10.0
        sasa, points = sasa_ext.compute_sasa(coords, radii, probe_radius=huge_probe, n_samples=500)

        # Should approximate convex hull behavior
        expected_radius = radii[0] + huge_probe
        expected_area = 4 * np.pi * expected_radius**2

        relative_error = abs(sasa - expected_area) / expected_area
        assert relative_error < 0.05, f"Large probe radius error: {relative_error*100:.2f}%"

    def test_memory_stress(self, sasa_ext_available):
        """Test with configurations that stress memory usage."""
        import sasa_ext

        # Many atoms with many samples
        n_atoms = 50
        coords = np.random.rand(n_atoms, 3).astype(np.float32) * 20
        radii = np.full(n_atoms, 1.5, dtype=np.float32)

        # High sample count
        sasa, points = sasa_ext.compute_sasa(coords, radii, n_samples=2000, seed=42)
        assert sasa > 0
        assert len(points) <= n_atoms * 2000  # Upper bound on points

    def test_reproducibility_different_seeds(self, sasa_ext_available):
        """Test that different seeds give different results."""
        import sasa_ext

        # Use fewer samples or more complex geometry where variance is expected
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # Same parameters, different seeds - use fewer samples for more variance
        sasa1, points1 = sasa_ext.compute_sasa(coords, radii, n_samples=50, seed=42)
        sasa2, points2 = sasa_ext.compute_sasa(coords, radii, n_samples=50, seed=12345)

        # Results should be different (with very high probability)
        # Check if SASA values are different or if point counts are different
        sasa_different = abs(sasa1 - sasa2) > 1e-6
        point_count_different = len(points1) != len(points2)

        # If same number of points, check if coordinates are different
        coords_different = False
        if len(points1) == len(points2):
            coords_different = not np.allclose(points1, points2, atol=1e-6)

        assert sasa_different or point_count_different or coords_different, \
            "Different seeds should give different results"

        # But both should be reasonable
        assert sasa1 > 0 and sasa2 > 0
        assert len(points1) > 0 and len(points2) > 0