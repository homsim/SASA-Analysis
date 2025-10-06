"""
Performance and numerical stability tests.

Tests performance scaling, numerical stability, and memory usage of the SASA implementation.
"""

import pytest
import numpy as np
import time
import gc
import sasa_ext

class TestPerformanceBenchmarks:
    """Test performance characteristics and scaling."""

    def test_performance_scaling_with_atoms(self):
        """Test that algorithm scales reasonably with atom count."""

        atom_counts = [10, 25, 50, 100]
        times = []
        sasa_per_atom = []

        for n_atoms in atom_counts:
            # Create random configuration in a box
            np.random.seed(42)  # For reproducibility
            coords = np.random.rand(n_atoms, 3).astype(np.float32) * 30  # 30Å box
            radii = np.full(n_atoms, 1.5, dtype=np.float32)

            # Time the computation
            start_time = time.time()
            sasa, points = sasa_ext.compute_sasa(
                coords, radii, probe_radius=1.4, n_samples=200, seed=42
            )
            elapsed = time.time() - start_time

            times.append(elapsed)
            sasa_per_atom.append(sasa / n_atoms)

            # Basic sanity checks
            assert sasa > 0, f"SASA should be positive for {n_atoms} atoms"
            assert len(points) > 0, f"Should generate points for {n_atoms} atoms"

        # Performance should scale reasonably (not worse than quadratic)
        if len(times) >= 2:
            # Check scaling between smallest and largest
            time_ratio = times[-1] / times[0]
            atom_ratio = atom_counts[-1] / atom_counts[0]

            # Should be better than O(n^2) scaling
            max_expected_ratio = atom_ratio ** 1.5  # Allow up to O(n^1.5)
            assert time_ratio < max_expected_ratio, \
                f"Poor scaling: {time_ratio:.1f}x time for {atom_ratio:.1f}x atoms"

        # Should complete in reasonable time
        max_time_per_atom = 0.1  # 0.1 seconds per atom is reasonable upper bound
        for i, (n_atoms, elapsed) in enumerate(zip(atom_counts, times)):
            time_per_atom = elapsed / n_atoms
            assert time_per_atom < max_time_per_atom, \
                f"Too slow for {n_atoms} atoms: {time_per_atom:.3f}s per atom"

    def test_performance_scaling_with_samples(self):
        """Test that performance scales linearly with sample count."""

        coords = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        sample_counts = [100, 200, 500, 1000]
        times = []

        for n_samples in sample_counts:
            start_time = time.time()
            sasa, points = sasa_ext.compute_sasa(
                coords, radii, probe_radius=1.4, n_samples=n_samples, seed=42
            )
            elapsed = time.time() - start_time
            times.append(elapsed)

            assert sasa > 0, f"SASA should be positive for {n_samples} samples"

        # Should scale roughly linearly with sample count
        if len(times) >= 2:
            time_ratio = times[-1] / times[0]
            sample_ratio = sample_counts[-1] / sample_counts[0]

            # Allow some overhead, but should be roughly linear
            assert time_ratio < 1.5 * sample_ratio, \
                f"Sample scaling poor: {time_ratio:.1f}x time for {sample_ratio:.1f}x samples"

    def test_memory_usage_reasonable(self):
        """Test that memory usage is reasonable."""

        # Large but manageable test case
        n_atoms = 200
        coords = np.random.rand(n_atoms, 3).astype(np.float32) * 50
        radii = np.full(n_atoms, 1.5, dtype=np.float32)

        # Force garbage collection before test
        gc.collect()

        # Run computation
        sasa, points = sasa_ext.compute_sasa(
            coords, radii, probe_radius=1.4, n_samples=500, seed=42
        )

        # Check results are reasonable
        assert sasa > 0
        assert len(points) <= n_atoms * 500  # Upper bound on points

        # Check point array size is reasonable
        points_memory_mb = points.nbytes / (1024 * 1024)
        max_expected_mb = 50  # 50MB should be more than enough
        assert points_memory_mb < max_expected_mb, \
            f"Surface points use too much memory: {points_memory_mb:.1f}MB"

    def test_large_molecule_performance(self):
        """Test performance on a large molecule."""

        # Simulate a large protein (500 atoms)
        n_atoms = 500

        # Create a more realistic distribution (not completely random)
        # Protein-like shape: roughly spherical with some density variation
        center = np.array([25.0, 25.0, 25.0])
        coords = []

        for i in range(n_atoms):
            # Gaussian distribution around center with some spread
            r = np.random.exponential(scale=15.0)  # Distance from center
            theta = np.random.uniform(0, 2*np.pi)  # Azimuthal angle
            phi = np.random.uniform(0, np.pi)      # Polar angle

            x = center[0] + r * np.sin(phi) * np.cos(theta)
            y = center[1] + r * np.sin(phi) * np.sin(theta)
            z = center[2] + r * np.cos(phi)

            coords.append([x, y, z])

        coords = np.array(coords, dtype=np.float32)

        # Mixed atom types (realistic protein composition)
        atom_types = np.random.choice(['C', 'N', 'O', 'H'], size=n_atoms, p=[0.4, 0.15, 0.25, 0.2])
        radii_map = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'H': 1.20}
        radii = np.array([radii_map[atom] for atom in atom_types], dtype=np.float32)

        # Time the computation
        start_time = time.time()
        sasa, points = sasa_ext.compute_sasa(
            coords, radii, probe_radius=1.4, n_samples=500, seed=42
        )
        elapsed = time.time() - start_time

        # Should complete in reasonable time (adjust threshold as needed)
        max_time = 10.0  # 10 seconds for 500 atoms
        assert elapsed < max_time, f"Too slow for large molecule: {elapsed:.2f}s"

        # Results should be reasonable
        assert sasa > 0, "Large molecule should have positive SASA"
        assert len(points) > 0, "Large molecule should generate surface points"

        # SASA per atom should be reasonable for a compact structure
        sasa_per_atom = sasa / n_atoms
        assert 10 < sasa_per_atom < 200, f"SASA per atom unreasonable: {sasa_per_atom:.1f} Ų"


class TestNumericalStability:
    """Test numerical stability and robustness."""

    def test_repeated_computation_stability(self, convergence_tolerances):
        """Test that repeated computations are stable."""

        coords = np.array([[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # Run the same computation multiple times
        results = []
        for i in range(10):
            sasa, points = sasa_ext.compute_sasa(
                coords, radii, probe_radius=1.4, n_samples=1000, seed=42  # Same seed
            )
            results.append(sasa)

        # All results should be identical (same seed)
        for result in results[1:]:
            assert abs(result - results[0]) < convergence_tolerances['reproducibility'], \
                "Results should be identical with same seed"

        # Now test with different seeds - should have low variance
        different_seed_results = []
        for seed in range(100, 110):  # 10 different seeds
            sasa, _ = sasa_ext.compute_sasa(
                coords, radii, probe_radius=1.4, n_samples=1000, seed=seed
            )
            different_seed_results.append(sasa)

        # Coefficient of variation should be small
        mean_sasa = np.mean(different_seed_results)
        std_sasa = np.std(different_seed_results)
        cv = std_sasa / mean_sasa

        assert cv < 0.02, f"High variance across seeds: CV = {cv:.4f}"

    def test_coordinate_precision_stability(self):
        """Test stability with different coordinate precisions."""

        base_coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # Test with different precision levels
        precisions = [1e-6, 1e-3, 1.0]  # Perturbations
        results = []

        for precision in precisions:
            # Add small perturbations
            perturbed_coords = base_coords + np.random.uniform(
                -precision, precision, base_coords.shape
            )
            coords_f32 = perturbed_coords.astype(np.float32)

            sasa, _ = sasa_ext.compute_sasa(
                coords_f32, radii, probe_radius=1.4, n_samples=1000, seed=42
            )
            results.append(sasa)

        # Results should be similar for small perturbations
        base_result = results[0]
        for i, precision in enumerate(precisions[1:], 1):
            relative_change = abs(results[i] - base_result) / base_result
            expected_stability = precision * 10  # Allow some amplification

            if precision < 0.1:  # Only test for small perturbations
                assert relative_change < expected_stability, \
                    f"Unstable with {precision} precision: {relative_change:.4f} change"

    def test_extreme_configurations_stability(self):
        """Test stability with extreme but valid configurations."""

        # Very close atoms
        close_coords = np.array([[0.0, 0.0, 0.0], [1e-3, 0.0, 0.0]], dtype=np.float32)
        radii_close = np.array([1.0, 1.0], dtype=np.float32)

        sasa_close, points_close = sasa_ext.compute_sasa(
            close_coords, radii_close, probe_radius=1.4, n_samples=500, seed=42
        )
        assert sasa_close >= 0, "SASA should be non-negative for close atoms"

        # Very far atoms
        far_coords = np.array([[0.0, 0.0, 0.0], [1000.0, 0.0, 0.0]], dtype=np.float32)
        radii_far = np.array([1.0, 1.0], dtype=np.float32)

        sasa_far, points_far = sasa_ext.compute_sasa(
            far_coords, radii_far, probe_radius=1.4, n_samples=500, seed=42
        )

        # Should be approximately 2x single atom
        single_sasa, _ = sasa_ext.compute_sasa(
            close_coords[:1], radii_close[:1], probe_radius=1.4, n_samples=500, seed=42
        )
        expected_far = 2 * single_sasa
        relative_error = abs(sasa_far - expected_far) / expected_far
        assert relative_error < 0.05, f"Far atoms don't behave independently: {relative_error:.3f}"

        # Very different radii
        mixed_coords = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float32)
        mixed_radii = np.array([0.1, 5.0], dtype=np.float32)  # Very different sizes

        sasa_mixed, points_mixed = sasa_ext.compute_sasa(
            mixed_coords, mixed_radii, probe_radius=1.4, n_samples=500, seed=42
        )
        assert sasa_mixed >= 0, "SASA should be non-negative for mixed radii"

    def test_convergence_with_sample_size(self):
        """Test that results converge properly with increasing sample size."""

        # Use two atoms to create a system where variance is expected
        coords = np.array([[0.0, 0.0, 0.0], [3.2, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        # Test convergence to analytical result
        sample_counts = [100, 500, 1000, 2000, 5000]
        results = []

        for n_samples in sample_counts:
            # Average over multiple seeds to reduce noise
            sasa_values = []
            for seed in range(10):
                sasa, _ = sasa_ext.compute_sasa(
                    coords, radii, probe_radius=0.0, n_samples=n_samples, seed=seed
                )
                sasa_values.append(sasa)

            mean_sasa = np.mean(sasa_values)
            std_sasa = np.std(sasa_values)
            results.append((mean_sasa, std_sasa))

        # For two atoms, just check that values are reasonable and converging
        # (analytical result is complex for this geometry)
        first_result = results[0][0]
        last_result = results[-1][0]

        for i, (mean_sasa, std_sasa) in enumerate(results):
            # Basic sanity checks
            assert mean_sasa > 0, f"SASA should be positive at {sample_counts[i]} samples"

            # Results should be reasonably close to first measurement
            relative_diff = abs(mean_sasa - first_result) / first_result
            assert relative_diff < 0.1, f"Results varying too much: {relative_diff:.4f}"

        # Standard deviation should generally decrease or stay low
        stds = [result[1] for result in results]

        # Check that we have some reasonable variance in the first few measurements
        # (if std is 0 everywhere, the test setup might be wrong)
        max_std = max(stds)
        if max_std > 0:
            # If we have variance, check that it generally decreases
            for i in range(1, len(stds)):
                if stds[i-1] > 1e-10:  # Avoid division by very small numbers
                    ratio = stds[i] / stds[i-1]
                    # Standard deviation should not increase significantly
                    assert ratio < 2.0, \
                        f"Standard deviation increased too much: {stds[i]:.6f} vs {stds[i-1]:.6f}"


class TestConcurrencyAndThreadSafety:
    """Test behavior under concurrent usage (if applicable)."""

    @pytest.mark.parametrize("n_threads", [5, 7, 8, 15])
    def test_multiple_simultaneous_calls(self, n_threads):
        """Test that multiple simultaneous calls don't interfere."""
        import threading
        import queue

        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        results = queue.Queue()
        errors = queue.Queue()

        def compute_sasa_thread(thread_id):
            try:
                sasa, points = sasa_ext.compute_sasa(
                    coords, radii, probe_radius=1.4, n_samples=500, seed=thread_id
                )
                results.put((thread_id, sasa, len(points)))
            except Exception as e:
                errors.put((thread_id, str(e)))

        # Start multiple threads
        threads = []

        for i in range(n_threads):
            thread = threading.Thread(target=compute_sasa_thread, args=(i,))
            threads.append(thread)
            thread.start()

        # Wait for all to complete
        for thread in threads:
            thread.join()

        # Check that all completed successfully
        assert errors.empty(), f"Errors occurred: {list(errors.queue)}"
        assert results.qsize() == n_threads, f"Not all threads completed: {results.qsize()}"

        # All results should be reasonable
        while not results.empty():
            thread_id, sasa, n_points = results.get()
            assert sasa > 0, f"Thread {thread_id} got invalid SASA: {sasa}"
            assert n_points > 0, f"Thread {thread_id} got no points"
