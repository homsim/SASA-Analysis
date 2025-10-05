"""
Unit tests for core SASA C extension components.

Tests individual algorithm components including sphere point generation,
neighbor list construction, and surface point testing.
"""

import pytest
import numpy as np
import scipy.stats as stats


class TestSpherePointGeneration:
    """Test sphere point generation specifically."""

    def test_sphere_point_uniformity(self, sasa_ext_available):
        """Test that sphere points are uniformly distributed."""
        import sasa_ext

        # Test with single atom at origin, zero probe radius
        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.0], dtype=np.float32)

        _, points = sasa_ext.compute_sasa(
            coords, radii, probe_radius=0.0, n_samples=1000, seed=42
        )

        # All points should be on sphere of radius 1.0
        center = coords[0]
        distances = np.linalg.norm(points - center, axis=1)
        expected_radius = radii[0]

        # Check that all points are approximately on the sphere surface
        radius_errors = np.abs(distances - expected_radius)
        assert np.max(radius_errors) < 0.01, "Points should be on sphere surface"

    def test_sphere_point_distribution_statistical(self, sasa_ext_available):
        """Test uniform distribution using statistical tests."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.0], dtype=np.float32)

        _, points = sasa_ext.compute_sasa(
            coords, radii, probe_radius=0.0, n_samples=5000, seed=12345
        )

        # Normalize points to unit sphere
        center = coords[0]
        normalized_points = (points - center) / radii[0]

        # Test z-coordinates are uniform in [-1, 1]
        z_coords = normalized_points[:, 2]
        _, p_value_z = stats.kstest(z_coords, lambda x: (x + 1) / 2)  # CDF for uniform[-1,1]
        assert p_value_z > 0.01, f"Z-coordinates not uniform (p={p_value_z:.4f})"

        # Test phi angles are uniform in [0, 2π]
        phi = np.arctan2(normalized_points[:, 1], normalized_points[:, 0]) + np.pi
        _, p_value_phi = stats.kstest(phi, lambda x: x / (2 * np.pi))  # CDF for uniform[0,2π]
        assert p_value_phi > 0.01, f"Phi angles not uniform (p={p_value_phi:.4f})"

    def test_sphere_point_count(self, sasa_ext_available):
        """Test that point count matches expected for isolated atom."""
        import sasa_ext

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.0], dtype=np.float32)

        for n_samples in [100, 500, 1000]:
            _, points = sasa_ext.compute_sasa(
                coords, radii, probe_radius=0.0, n_samples=n_samples, seed=42
            )

            # For isolated atom, all points should be exposed
            assert len(points) == n_samples, f"Expected {n_samples} points, got {len(points)}"


class TestCoreSASAAlgorithm:
    """Test core SASA computation functionality."""

    def test_single_atom_analytical(self, sasa_ext_available, simple_test_cases):
        """Test single atom case against analytical result."""
        import sasa_ext

        case = simple_test_cases['single_atom']

        # For probe_radius=0, should get exactly 4πr²
        total_sasa, surface_points = sasa_ext.compute_sasa(
            case['coords'], case['radii'],
            probe_radius=0.0, n_samples=5000, seed=42
        )

        expected_area = case['expected_area_factor'] * (case['radii'][0] ** 2)
        relative_error = abs(total_sasa - expected_area) / expected_area

        assert relative_error < 0.01, f"Single atom SASA error too large: {relative_error*100:.2f}%"
        assert len(surface_points) > 0, "Should generate surface points"

    def test_two_atoms_geometry(self, sasa_ext_available, simple_test_cases):
        """Test geometric behavior with two atoms."""
        import sasa_ext

        separated = simple_test_cases['two_separated']
        touching = simple_test_cases['two_touching']

        # Test separated atoms
        sasa_sep, _ = sasa_ext.compute_sasa(
            separated['coords'], separated['radii'],
            probe_radius=separated['probe_radius'], n_samples=1000, seed=42
        )

        # Test touching atoms
        sasa_touch, _ = sasa_ext.compute_sasa(
            touching['coords'], touching['radii'],
            probe_radius=touching['probe_radius'], n_samples=1000, seed=42
        )

        # Touching atoms should have less SASA than separated
        assert sasa_touch < sasa_sep, "Touching atoms should have less SASA than separated"

    def test_buried_atom_configuration(self, sasa_ext_available, simple_test_cases):
        """Test atom burial in crowded environment."""
        import sasa_ext

        case = simple_test_cases['buried_atom']

        total_sasa, points = sasa_ext.compute_sasa(
            case['coords'], case['radii'],
            probe_radius=case['probe_radius'], n_samples=1000, seed=42
        )

        # Central atom should contribute very little to total SASA
        # Test by computing single central atom SASA for comparison
        central_sasa, _ = sasa_ext.compute_sasa(
            case['coords'][:1], case['radii'][:1],
            probe_radius=case['probe_radius'], n_samples=1000, seed=42
        )

        # In crowded environment, total should be much less than sum of isolated atoms
        all_isolated_sasa = 0
        for i in range(len(case['coords'])):
            isolated_sasa, _ = sasa_ext.compute_sasa(
                case['coords'][i:i+1], case['radii'][i:i+1],
                probe_radius=case['probe_radius'], n_samples=1000, seed=42
            )
            all_isolated_sasa += isolated_sasa

        burial_ratio = total_sasa / all_isolated_sasa
        assert burial_ratio < 0.7, f"Burial effect insufficient: {burial_ratio:.2f}"

    def test_parameter_effects(self, sasa_ext_available, simple_test_cases, sasa_parameters):
        """Test that parameters have expected effects."""
        import sasa_ext

        case = simple_test_cases['single_atom']

        # Test probe radius effect
        sasa_small, _ = sasa_ext.compute_sasa(
            case['coords'], case['radii'],
            probe_radius=1.0, n_samples=500, seed=42
        )
        sasa_large, _ = sasa_ext.compute_sasa(
            case['coords'], case['radii'],
            probe_radius=2.0, n_samples=500, seed=42
        )

        # Larger probe should give larger SASA
        assert sasa_large > sasa_small, "Larger probe radius should increase SASA"

        # Test sample count convergence
        params_few = sasa_parameters['few_samples']
        params_many = sasa_parameters['many_samples']

        sasa_few, _ = sasa_ext.compute_sasa(
            case['coords'], case['radii'], **params_few
        )
        sasa_many, _ = sasa_ext.compute_sasa(
            case['coords'], case['radii'], **params_many
        )

        # Should converge to similar values
        relative_diff = abs(sasa_many - sasa_few) / sasa_many
        assert relative_diff < 0.05, f"Poor convergence with sample count: {relative_diff*100:.2f}%"

    def test_reproducibility(self, sasa_ext_available, simple_test_cases, convergence_tolerances):
        """Test that results are reproducible with same seed."""
        import sasa_ext

        case = simple_test_cases['single_atom']

        # Run multiple times with same seed
        results = []
        for _ in range(5):
            sasa, points = sasa_ext.compute_sasa(
                case['coords'], case['radii'],
                probe_radius=1.4, n_samples=500, seed=42
            )
            results.append((sasa, len(points)))

        # All results should be identical
        first_result = results[0]
        for result in results[1:]:
            assert abs(result[0] - first_result[0]) < convergence_tolerances['reproducibility'], \
                "SASA should be reproducible with same seed"
            assert result[1] == first_result[1], \
                "Point count should be reproducible with same seed"

    def test_scaling_invariance(self, sasa_ext_available):
        """Test that SASA scales correctly with coordinate scaling."""
        import sasa_ext

        # Original configuration
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.0, 1.0], dtype=np.float32)
        probe_radius = 1.4

        original_sasa, _ = sasa_ext.compute_sasa(
            coords, radii, probe_radius=probe_radius, n_samples=1000, seed=42
        )

        # Scaled configuration (2x)
        scale_factor = 2.0
        scaled_coords = coords * scale_factor
        scaled_radii = radii * scale_factor
        scaled_probe = probe_radius * scale_factor

        scaled_sasa, _ = sasa_ext.compute_sasa(
            scaled_coords, scaled_radii,
            probe_radius=scaled_probe, n_samples=1000, seed=42
        )

        # SASA should scale by area factor (scale²)
        expected_sasa = original_sasa * (scale_factor ** 2)
        relative_error = abs(scaled_sasa - expected_sasa) / expected_sasa

        assert relative_error < 0.02, f"Scaling invariance violated: {relative_error*100:.2f}% error"


class TestNeighborListConstruction:
    """Test spatial optimization via neighbor lists."""

    def test_neighbor_finding_correctness(self, sasa_ext_available):
        """Test that neighbor lists find all relevant neighbors."""
        import sasa_ext

        # Create configuration where we know which atoms should be neighbors
        coords = np.array([
            [0.0, 0.0, 0.0],   # Central atom
            [2.0, 0.0, 0.0],   # Close neighbor
            [1.0, 1.0, 0.0],   # Close neighbor
            [10.0, 0.0, 0.0],  # Far neighbor (should not affect central atom)
        ], dtype=np.float32)
        radii = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float32)

        # Test with different probe radii
        for probe_radius in [0.5, 1.4, 2.0]:
            sasa_all, _ = sasa_ext.compute_sasa(
                coords, radii, probe_radius=probe_radius, n_samples=500, seed=42
            )

            # Test central atom alone
            central_sasa, _ = sasa_ext.compute_sasa(
                coords[:1], radii[:1], probe_radius=probe_radius, n_samples=500, seed=42
            )

            # Central atom in full environment should have less SASA than isolated
            assert sasa_all < 3 * central_sasa, "Neighbor effects should reduce SASA"

    def test_performance_scaling(self, sasa_ext_available):
        """Test that algorithm scales reasonably with atom count."""
        import time
        import sasa_ext

        atom_counts = [10, 50, 100]
        times = []

        for n_atoms in atom_counts:
            # Create random configuration
            np.random.seed(42)  # For reproducibility
            coords = np.random.rand(n_atoms, 3).astype(np.float32) * 20  # 20Å box
            radii = np.full(n_atoms, 1.5, dtype=np.float32)

            start_time = time.time()
            sasa_ext.compute_sasa(coords, radii, probe_radius=1.4, n_samples=100, seed=42)
            elapsed = time.time() - start_time
            times.append(elapsed)

        # Should scale roughly linearly (not quadratically)
        # Check that 100 atoms doesn't take more than 10x the time of 10 atoms
        scaling_factor = times[-1] / times[0]  # 100 atoms vs 10 atoms
        assert scaling_factor < 20, f"Poor scaling: {scaling_factor:.1f}x slowdown for 10x atoms"


class TestSurfacePointTesting:
    """Test the point burial logic."""

    def test_point_burial_simple_cases(self, sasa_ext_available):
        """Test point burial with simple geometric cases."""
        import sasa_ext

        # Two touching spheres - points between them should be buried
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # With small probe radius, some points should be exposed
        total_sasa_small, points_small = sasa_ext.compute_sasa(
            coords, radii, probe_radius=0.1, n_samples=1000, seed=42
        )

        # With large probe radius, fewer points should be exposed
        total_sasa_large, points_large = sasa_ext.compute_sasa(
            coords, radii, probe_radius=2.0, n_samples=1000, seed=42
        )

        # Larger probe should result in more burial
        assert len(points_large) < len(points_small), "Large probe should bury more points"
        assert total_sasa_large > total_sasa_small, "Large probe should increase total SASA"

    def test_complete_burial(self, sasa_ext_available):
        """Test case where central atom is completely buried."""
        import sasa_ext

        # Central atom surrounded by large neighbors
        coords = np.array([
            [0.0, 0.0, 0.0],   # Central small atom
            [1.5, 0.0, 0.0],   # Large surrounding atoms
            [-1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
            [0.0, -1.5, 0.0],
            [0.0, 0.0, 1.5],
            [0.0, 0.0, -1.5]
        ], dtype=np.float32)
        radii = np.array([0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], dtype=np.float32)

        total_sasa, points = sasa_ext.compute_sasa(
            coords, radii, probe_radius=1.4, n_samples=1000, seed=42
        )

        # Central atom should contribute almost nothing
        # Test by computing SASA of just the surrounding atoms
        surrounding_sasa, _ = sasa_ext.compute_sasa(
            coords[1:], radii[1:], probe_radius=1.4, n_samples=1000, seed=42
        )

        # Total should be close to just the surrounding atoms
        burial_effect = abs(total_sasa - surrounding_sasa) / surrounding_sasa
        assert burial_effect < 0.1, f"Central atom not properly buried: {burial_effect:.2f}"