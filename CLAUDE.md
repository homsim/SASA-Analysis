- The code in this directory does some scientific calculations by using the so called solvent-accessible surface analysis (SASA). The directory @vmd-python/ is a library that I cloned here for you to see its source code. I use this library by installing in a conda environment. What I want to do eventually is to remove the need for this library. I want to do this because it is very bloated for my purposes. In my code (which you find in @sasa_lammps/ ) I really only use one major class and function of vmd-python, namely the calculation a sasa() on an AtomSel object (atom selection). I would like to implement this as a C extension and use it in python, because I expect it be too CPU-heavy to implement it in python. First, try to understand, how the code from vmd-python works. You can do this by tracking back how I use it in @sasa_lammps/conversion.py in the method _create_sasa_xyz()

## SASA Algorithm Analysis (from vmd-python)

### Core Algorithm Overview
The SASA calculation uses **numerical integration via Monte Carlo sampling** on spheres around each atom.

### Implementation Details (`vmd/vmd_src/src/Measure.C:measure_sasa()`):

1. **Sphere Point Generation** (lines 1431-1441):
   - Generates `npts` random points on unit sphere using uniform sampling
   - Uses spherical coordinates: `z = 2*u1 - 1`, `phi = 2*π*u2`
   - Seeds RNG (38572111) for reproducible results
   - Default: 500 points per atom (converges to ~1% accuracy)

2. **Spatial Optimization** (lines 1408-1421):
   - Grid search finds neighbors within `2.0 * (maxrad + srad)`
   - Creates pairlist per atom to avoid O(n²) all-vs-all checks

3. **Surface Point Testing** (lines 1455-1479):
   - For each atom: scale unit sphere by `(atom_radius + probe_radius)`
   - Test each surface point against all neighbors
   - Point is buried if inside any neighbor's expanded radius
   - Count exposed points for area calculation

4. **Area Calculation** (line 1480):
   - `Surface_area = (4π/npts) * radius² * exposed_points`
   - Sum over all atoms

### Usage in sasa_lammps/conversion.py:59:
```python
_, sasa_points = sel.sasa(srad=srad, samples=samples, points=True)
```
- `srad`: probe radius (typically 1.4Å for water)
- `samples`: number of Monte Carlo points per atom
- `points=True`: returns 3D coordinates of surface points

### Key Parameters:
- **Default samples**: 500 (good accuracy/performance balance)
- **Probe radius**: Added to VdW radii for solvent accessibility
- **Thread support**: Available for multiple selections
- **Random seed**: Fixed for reproducibility

## Minimal C Extension Components for SASA Replacement

### Core Components Needed:

1. **Core Algorithm Function**
```c
int compute_sasa(
    const float *coords,      // atom coordinates [3*n_atoms]
    const float *radii,       // atom radii [n_atoms]
    int n_atoms,
    float probe_radius,       // srad parameter
    int n_samples,           // samples parameter
    float *sasa_area,        // output: total SASA area
    float **surface_points,  // output: surface point coordinates
    int *n_surface_points    // output: number of surface points
);
```

2. **Spatial Optimization (Critical for Performance)**
```c
// Grid-based neighbor finding (from Measure.C:1408-1421)
typedef struct {
    int *neighbors;
    int n_neighbors;
} neighbor_list_t;

neighbor_list_t* build_neighbor_lists(
    const float *coords,
    int n_atoms,
    float cutoff
);
```

3. **Random Sphere Point Generator**
```c
// Generate uniform random points on unit sphere (Measure.C:1431-1441)
void generate_sphere_points(
    float *sphere_points,  // output [3*n_samples]
    int n_samples,
    unsigned int seed      // for reproducibility
);
```

4. **Python C Extension Interface**
```c
// Python module interface
static PyObject* py_compute_sasa(PyObject *self, PyObject *args) {
    // Parse numpy arrays for coords, radii
    // Call compute_sasa()
    // Return (sasa_area, surface_points_array)
}
```

### Implementation Size Estimate:
- **Core SASA algorithm**: ~150 lines
- **Grid neighbor search**: ~100 lines
- **Sphere point generation**: ~50 lines
- **Python bindings**: ~100 lines
- **Total**: ~400 lines of C code

### What You DON'T Need from VMD:
- AtomSel class hierarchy
- Molecule/MoleculeList classes
- Graphics/visualization code
- TCL bindings
- Multi-frame handling
- File I/O systems
- Threading (initially - can add later)

### Simplified Python Interface:
```python
import sasa_ext

# Your current usage:
# _, sasa_points = sel.sasa(srad=srad, samples=samples, points=True)

# Becomes:
coords = atom_coordinates  # numpy array (n_atoms, 3)
radii = atom_radii        # numpy array (n_atoms,)
sasa_area, sasa_points = sasa_ext.compute_sasa(
    coords, radii, srad=srad, samples=samples
)
```

### Dependencies:
- Standard C math library
- Python C API
- NumPy C API (for array handling)

This implementation would be roughly **5-10% the size** of the full VMD codebase while providing exactly the functionality needed.