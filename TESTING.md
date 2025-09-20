
## Testing Strategy for SASA Implementation

### 1. **Fixed Seed Reference Tests**

Since VMD uses a fixed seed (38572111), you can create deterministic reference data:

```python
import numpy as np
from vmd import molecule, measure
import pytest

class TestSASAAgainstVMD:
    @pytest.fixture
    def reference_structures(self):
        """Load test structures and calculate reference SASA values"""
        test_cases = []
        
        # Simple test cases with known geometry
        structures = [
            'single_atom.pdb',      # Trivial case
            'two_atoms.pdb',        # Simple overlap
            'small_protein.pdb',    # Real-world small case
            'buried_residue.pdb',   # Edge case with buried atoms
        ]
        
        for structure_file in structures:
            mol = molecule.load('pdb', structure_file)
            
            # Test with different parameters
            for npts in [100, 500, 1000]:
                for probe_radius in [1.4, 0.0, 2.8]:
                    sasa = measure.sasa(
                        selection=mol.atoms,
                        probe_radius=probe_radius,
                        samples=npts
                    )
                    
                    test_cases.append({
                        'file': structure_file,
                        'npts': npts,
                        'probe_radius': probe_radius,
                        'total_sasa': sasa.sum(),
                        'per_atom_sasa': sasa,
                        'params': (npts, probe_radius)
                    })
        
        return test_cases
    
    def test_total_sasa_matches(self, reference_structures):
        """Test that total SASA matches VMD within tolerance"""
        for ref_case in reference_structures:
            result = your_sasa_module.calculate_sasa(
                ref_case['file'],
                probe_radius=ref_case['probe_radius'],
                n_points=ref_case['npts'],
                seed=38572111  # Match VMD's seed
            )
            
            # Total SASA should match closely
            np.testing.assert_allclose(
                result['total_sasa'],
                ref_case['total_sasa'],
                rtol=1e-3,  # 0.1% tolerance
                err_msg=f"Failed for {ref_case['file']} with params {ref_case['params']}"
            )
```

### 2. **Statistical Convergence Tests**

Test that your Monte Carlo sampling converges properly:

```python
def test_monte_carlo_convergence():
    """Test that SASA converges as n_points increases"""
    test_structure = 'test_protein.pdb'
    probe_radius = 1.4
    
    # Calculate SASA with increasing sample points
    n_points_list = [50, 100, 200, 500, 1000, 2000, 5000]
    sasa_values = []
    
    for n_points in n_points_list:
        # Run multiple times to get statistics
        runs = []
        for seed in range(10):
            sasa = your_sasa_module.calculate_sasa(
                test_structure,
                probe_radius=probe_radius,
                n_points=n_points,
                seed=seed
            )
            runs.append(sasa['total_sasa'])
        
        sasa_values.append({
            'n_points': n_points,
            'mean': np.mean(runs),
            'std': np.std(runs),
            'cv': np.std(runs) / np.mean(runs)  # Coefficient of variation
        })
    
    # Check convergence properties
    for i in range(1, len(sasa_values)):
        # Standard deviation should decrease
        assert sasa_values[i]['std'] <= sasa_values[i-1]['std'] * 1.1
        
        # Coefficient of variation should decrease roughly as 1/sqrt(n)
        if sasa_values[i]['n_points'] >= 500:
            assert sasa_values[i]['cv'] < 0.02  # Should be < 2% for 500+ points
    
    # Mean should stabilize
    high_n_means = [s['mean'] for s in sasa_values if s['n_points'] >= 1000]
    assert np.std(high_n_means) / np.mean(high_n_means) < 0.005  # < 0.5% variation
```

### 3. **Geometric Invariant Tests**

```python
def test_geometric_invariants():
    """Test properties that must hold regardless of sampling"""
    
    # Test 1: Single isolated atom
    def test_isolated_atom():
        # For a single atom with no neighbors
        atom_radius = 1.5
        probe_radius = 1.4
        expected_sasa = 4 * np.pi * (atom_radius + probe_radius) ** 2
        
        result = your_sasa_module.calculate_single_atom_sasa(
            atom_radius, probe_radius, n_points=5000
        )
        
        np.testing.assert_allclose(result, expected_sasa, rtol=0.01)
    
    # Test 2: Fully buried atom
    def test_buried_atom():
        # Create configuration where central atom is completely buried
        result = your_sasa_module.calculate_sasa('buried_atom.pdb')
        
        # Central atom should have zero SASA
        assert abs(result['per_atom_sasa'][central_atom_idx]) < 1e-6
    
    # Test 3: Scaling invariance
    def test_scaling():
        original = your_sasa_module.calculate_sasa('molecule.pdb', probe_radius=1.4)
        
        # Scale all coordinates by factor of 2
        scaled = your_sasa_module.calculate_sasa('molecule_scaled_2x.pdb', probe_radius=2.8)
        
        # SASA should scale by factor of 4 (area scaling)
        np.testing.assert_allclose(
            scaled['total_sasa'],
            original['total_sasa'] * 4,
            rtol=0.02
        )
```

### 4. **Algorithm Component Tests**

Test individual components of your algorithm:

```python
class TestAlgorithmComponents:
    def test_sphere_point_generation(self):
        """Test uniform distribution of points on sphere"""
        n_points = 10000
        points = your_sasa_module._generate_sphere_points(n_points, seed=42)
        
        # Check all points are on unit sphere
        radii = np.linalg.norm(points, axis=1)
        np.testing.assert_allclose(radii, 1.0, atol=1e-10)
        
        # Check uniform distribution using statistical tests
        # Test z-coordinates are uniform in [-1, 1]
        z_coords = points[:, 2]
        _, p_value = stats.kstest(z_coords, 'uniform', args=(-1, 2))
        assert p_value > 0.01  # Should not reject uniformity
        
        # Check phi angles are uniform in [0, 2π]
        phi = np.arctan2(points[:, 1], points[:, 0]) + np.pi
        _, p_value = stats.kstest(phi, 'uniform', args=(0, 2*np.pi))
        assert p_value > 0.01
    
    def test_neighbor_search(self):
        """Test spatial grid optimization"""
        # Create test configuration
        atoms = create_test_atoms()
        cutoff = 2.0 * (max_radius + probe_radius)
        
        # Get neighbors using your implementation
        neighbors = your_sasa_module._find_neighbors(atoms, cutoff)
        
        # Verify against brute force
        for i, atom in enumerate(atoms):
            brute_force_neighbors = []
            for j, other in enumerate(atoms):
                if i != j:
                    dist = np.linalg.norm(atom.pos - other.pos)
                    if dist < cutoff:
                        brute_force_neighbors.append(j)
            
            assert set(neighbors[i]) == set(brute_force_neighbors)
    
    def test_point_burial_check(self):
        """Test the point burial logic"""
        # Setup: atom at origin, one neighbor
        center = np.array([0, 0, 0])
        center_radius = 1.0
        probe_radius = 1.4
        
        neighbor_pos = np.array([3.0, 0, 0])
        neighbor_radius = 1.0
        
        # Test point that should be exposed
        exposed_point = np.array([0, 0, 1]) * (center_radius + probe_radius)
        assert your_sasa_module._is_exposed(
            exposed_point, neighbor_pos, neighbor_radius + probe_radius
        )
        
        # Test point that should be buried
        buried_point = np.array([1, 0, 0]) * (center_radius + probe_radius)
        assert not your_sasa_module._is_exposed(
            buried_point, neighbor_pos, neighbor_radius + probe_radius
        )
```

### 5. **Performance Regression Tests**

```python
def test_performance():
    """Ensure optimizations don't break"""
    import time
    
    test_cases = [
        ('small_protein.pdb', 100, 0.1),   # 100 atoms, should be < 0.1s
        ('medium_protein.pdb', 1000, 1.0),  # 1000 atoms, should be < 1s
        ('large_protein.pdb', 5000, 10.0),  # 5000 atoms, should be < 10s
    ]
    
    for pdb_file, n_atoms, max_time in test_cases:
        start = time.time()
        result = your_sasa_module.calculate_sasa(pdb_file, n_points=500)
        elapsed = time.time() - start
        
        assert elapsed < max_time, f"Too slow for {n_atoms} atoms: {elapsed:.2f}s"
        
        # Also verify O(n) scaling with neighbor lists
        # Should NOT be O(n²)
```

### 6. **Edge Case Tests**

```python
def test_edge_cases():
    """Test boundary conditions and edge cases"""
    
    # Zero probe radius
    result = your_sasa_module.calculate_sasa('protein.pdb', probe_radius=0.0)
    assert result['total_sasa'] > 0
    
    # Very large probe radius
    result = your_sasa_module.calculate_sasa('protein.pdb', probe_radius=10.0)
    # Should approximate convex hull
    
    # Single point sampling (n_points=1)
    result = your_sasa_module.calculate_sasa('protein.pdb', n_points=1)
    # Should still run, though inaccurate
    
    # Overlapping atoms
    result = your_sasa_module.calculate_sasa('overlapping_atoms.pdb')
    # Should handle gracefully
```

### 7. **Integration Test Suite**

```python
@pytest.mark.integration
def test_full_pipeline():
    """Test complete workflow matches VMD"""
    
    # Create a fixture that runs both implementations
    test_structures = glob.glob('test_data/*.pdb')
    
    for structure in test_structures:
        # Run VMD
        mol = molecule.load('pdb', structure)
        vmd_sasa = measure.sasa(mol.atoms, probe_radius=1.4, samples=500)
        
        # Run your implementation
        your_sasa = your_sasa_module.calculate_sasa(
            structure, 
            probe_radius=1.4, 
            n_points=500,
            seed=38572111
        )
        
        # Compare
        correlation = np.corrcoef(vmd_sasa, your_sasa['per_atom_sasa'])[0, 1]
        assert correlation > 0.999, f"Poor correlation for {structure}: {correlation}"
        
        # Check relative error
        rel_error = abs(your_sasa['total_sasa'] - vmd_sasa.sum()) / vmd_sasa.sum()
        assert rel_error < 0.01, f"Large error for {structure}: {rel_error*100:.2f}%"
```

### Test Data Generation Script

```python
def generate_test_data():
    """Generate reference data from VMD for testing"""
    import json
    
    test_configs = []
    for pdb_file in glob.glob('test_structures/*.pdb'):
        mol = molecule.load('pdb', pdb_file)
        
        for n_points in [100, 500, 1000]:
            sasa = measure.sasa(mol.atoms, samples=n_points)
            
            test_configs.append({
                'file': pdb_file,
                'n_points': n_points,
                'seed': 38572111,
                'total_sasa': float(sasa.sum()),
                'per_atom_sasa': sasa.tolist(),
                'n_atoms': len(sasa)
            })
    
    with open('test_reference_data.json', 'w') as f:
        json.dump(test_configs, f, indent=2)
```

This comprehensive testing approach ensures your C extension correctly replicates VMD's SASA calculation while being robust to the non-deterministic nature of Monte Carlo sampling.
