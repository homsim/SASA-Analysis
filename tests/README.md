# SASA Test Suite

This directory contains the comprehensive test suite for the SASA C extension implementation, following the testing strategy outlined in `../TESTING.md`.

## Test Structure

### Core Test Files

- **`test_core_components.py`** - Unit tests for core C extension components
  - Sphere point generation uniformity and distribution
  - Core SASA algorithm validation against analytical results
  - Neighbor list construction and optimization
  - Surface point burial logic
  - Parameter effects and geometric relationships

- **`test_input_validation.py`** - Input validation and error handling
  - Array shape and type validation
  - Edge cases and boundary conditions
  - Error handling for invalid inputs
  - Extreme parameter values

- **`test_vmd_integration.py`** - Integration tests against VMD reference
  - Comparison with VMD reference data (when available)
  - Algorithm fidelity validation
  - Parameter compatibility testing
  - Integration with replacement modules

- **`test_performance.py`** - Performance and numerical stability
  - Performance scaling with atom count and sample size
  - Memory usage validation
  - Numerical stability and reproducibility
  - Convergence testing

- **`test_python_modules.py`** - Python module testing
  - XYZ file parsing
  - VdW radius assignment
  - Integration with sasa_lammps package
  - Error propagation

### Support Files

- **`conftest.py`** - Pytest configuration and shared fixtures
- **`__init__.py`** - Test package initialization

## Running Tests

### Run All Tests
```bash
pytest tests/ -v
```

### Run Specific Test Categories
```bash
# Core algorithm tests only
pytest tests/test_core_components.py -v

# Performance tests only
pytest tests/test_performance.py -v

# Integration tests only
pytest tests/test_vmd_integration.py -v
```

### Run Tests by Marker
```bash
# Skip slow tests
pytest tests/ -m "not slow" -v

# Run only performance tests
pytest tests/ -m "performance" -v
```

### Test with Coverage
```bash
pytest tests/ --cov=sasa_ext --cov=sasa_lammps --cov-report=html
```

## Test Requirements

### Minimal Requirements
- `pytest` >= 6.0.0
- `numpy` >= 1.15.0
- `scipy` (for statistical tests)

### Full Testing Requirements
- SASA C extension compiled and importable
- VMD reference data (optional, for integration tests)
- All package dependencies installed

## Test Categories

### Unit Tests
Test individual components in isolation:
- Sphere point generation algorithms
- Neighbor list construction
- Point burial detection
- Input validation

### Integration Tests
Test complete workflows:
- End-to-end SASA computation
- XYZ file processing pipeline
- VMD compatibility validation
- Package integration

### Performance Tests
Validate performance characteristics:
- Scaling with problem size
- Memory usage patterns
- Numerical stability
- Convergence rates

### Edge Case Tests
Test boundary conditions:
- Extreme parameter values
- Malformed inputs
- Empty or minimal datasets
- Numerical edge cases

## Expected Test Results

### With SASA Extension Available
All tests should pass, demonstrating:
- ✅ Accurate SASA computation
- ✅ Proper input validation
- ✅ Good performance scaling
- ✅ Numerical stability
- ✅ VMD compatibility (if reference data available)

### Without SASA Extension
Tests requiring the C extension will be skipped:
- ⏭️ Core algorithm tests skipped
- ⏭️ Performance tests skipped
- ✅ Python module tests should still pass
- ✅ Input validation structure tests pass

## Troubleshooting

### Common Issues

1. **Extension Not Found**
   ```
   ModuleNotFoundError: No module named 'sasa_ext'
   ```
   Solution: Install package with `pip install -e .`

2. **VMD Reference Tests Skipped**
   ```
   SKIPPED [1] ... VMD reference data not available
   ```
   Solution: Generate reference data first (requires VMD)

3. **Performance Tests Failing**
   ```
   AssertionError: Too slow for large molecule
   ```
   Solution: May indicate compilation issues or slow system

### Debugging Failed Tests

1. **Increase Verbosity**
   ```bash
   pytest tests/test_core_components.py::TestCoreSASAAlgorithm::test_single_atom_analytical -vvv
   ```

2. **Run Single Test**
   ```bash
   pytest tests/test_core_components.py::TestCoreSASAAlgorithm::test_single_atom_analytical
   ```

3. **Check Test Output**
   ```bash
   pytest tests/ --tb=long --show-capture=all
   ```

## Contributing Tests

When adding new tests:

1. Follow existing naming conventions
2. Use appropriate fixtures from `conftest.py`
3. Add proper docstrings explaining test purpose
4. Include both positive and negative test cases
5. Use appropriate markers for categorization

### Test Naming Convention
- `test_<functionality>_<specific_case>`
- Example: `test_sphere_point_generation_uniformity`

### Fixture Usage
```python
def test_my_feature(sasa_ext_available, simple_test_cases, convergence_tolerances):
    """Test description."""
    # Test implementation
```

This test suite ensures the SASA implementation is robust, accurate, and maintains compatibility with existing VMD-based workflows.