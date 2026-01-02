# Overview

This package calcukates the solvent-accesible surface analysis (SASA) on a given macro-molecule. It places a given probe-molecule on the surface points to compute a 3D interaction energy surface. The potential energy is calculated with a ReaxFF potential that needs to be provided (see below) and uses [LAMMPS](https://www.lammps.org) as actual implementation for the energy calculations.

# Build and Install

This package uses a custom C extension for SASA calculations and can be installed as a standard Python package using pip.

## Requirements

- Python 3.9 or higher
- A C compiler (gcc, clang)
- Python development headers

## Installation

It is recommended to use a python environment.
```bash
# Create a python environment (optional) or use an existing one
python -m venv ~/venv/sasa
# Activate your Python environment
source ~/venv/sasa/bin/activate
```

From this projects root directory install the package with all dependencies
```
pip install .
```

## Verification

After installation, you can either verify that the SASA extension is loaded correctly (tp not confuse python with the local `sasa_ext` directory, execute this from any other than the source-file directory):
```bash
cd && python -c "import sasa_ext; print('✓ SASA extension loaded successfully')"
```

or/and execute the tests in the `tests` directory. This requires the installation of the test-requirements and the install in editable mode:
```bash
pip install -e ".[test]"
pytest ./tests
```

# Usage

The package really only has one usable method `sasa_lammps.sasa()`:

```python
>>> from sasa_lammps import sasa
```

```
>>> print(sasa.__doc__)
    Run the SASA (solvet accasible surface analysis) on a given macromolecule
    using a given probe molecule.
    The package was designed to start from a gromacs file of the macromolecule.
    For good simulation practices the macromolecule should be pre-equilibrated in water.
    Care must be taken for N-atomic probe molecules: The package does not identify
    a plane or something in the probe molecule. It just makes sure that at every
    interaction site the probe faces the macromolecule with the same orientation.
    However, the orientation itself is purely determined by the configuration
    given in the mol_file.

    Parameters
    ----------
    gro_file : str
        Name of the gromacs file of the macromolecule
    data_file : str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
        Name of the LAMMPS mol file to use as probe of the SAS (solvet accasible surface)
    ff_str : str
        Force field parameters to provide to LAMMPS. See examples directory
        https://docs.lammps.org/pair_style.html
        https://docs.lammps.org/pair_coeff.html
        Care must be taken because currently the 'unit real' in the in.template basically restricts to only use pair_style reaxff.
    dump_str : str
        Dump command to provide to LAMMPS. See examples directory
        https://docs.lammps.org/dump.html
    lammps_exe : str, optional
        Full path to the LAMMPS executable. If not provided, will automatically
        download and use pre-built LAMMPS binaries from https://download.lammps.org/static/
    n_procs : int, optional
        Number of LAMMPS instances to run in parallel (Default: 1)
    srad : float, optional
        Probe radius: Effectively a scaling factor for the vdW radii
        (Default: 1.4, which is the most commonly used because its approx. the
        radius of water)
    samples : int, optional
        Maximum points on the atomic vdW sphere to generate per atom (Default: 100)
    path : str, optional
        Execution path (Default: .)

    Returns
    -------
    None

```

## Example

In the directory `example` you find an example on the protein lysozyme with either H or H2O2 as probe molecules. There are two files lysozyme_part.gro and lysozyme_full.gro in the folder. The first contains only a part of the protein to speed up the test (sampling 2691 SAS points), the latter contains the full protein (sampling 8130 SAS points). 
To succesfully convert the .gro file to a lammps data file you need a library file of all elements included in you system, element_library.txt. This files should contain the ElementType (H,C,O, ...), the gromacs_ParticleType (1, 6, 8, ...), the lammps_ParticleType (1,2,3, ...) and the mass (1.00, 12.0, 16.0 ...). 
You should doubble check first which gromacs_ParticleType number ovito gives the elements. It is mostly related to the atomic number, but for Mg, for example, gromacs_ParticleType is 0. So check that first! 
You can choose any lammps_ParticleType you want for each of the elements, but the numbers have to be consecutive (1,2,3, ...).
For the ReaxFF simulation with lammps the elements have to be listed in the right order in the "pair_coeff" command. Resulting in: ```pair_coeff     * * protein2013.ff H C O```. If you use the protein2013.ff force field and H has the lammps_ParticleType 1, C is 2 and O is 3. 
The following files need to be provided for LAMMPS:
- lysozyme_part.gro (.gro file)
- h.mol or h2o2.mol (.mol file)
- element_library.txt
- protein2013.ff

```python
from sasa_lammps import sasa

# Import the gromacs file
gro_file = "lysozyme_part.gro"
# Choose a name for the lammps data file
data_file = "data.lysozyme_part"
# Import the molecule file (example also contains h2o2.mol)
mol_file = "h.mol"
# Path to your lammps executable (optional - will auto-download if not provided)
# lammps_exe =  "/opt/lammps-23Jun2022/src/lmp_mpi"
# Specify the force filed parameters, in this case reaxFF parameters*
ff_str = """
pair_style      reaxff NULL safezone 1.6 mincap 100 minhbonds 150
pair_coeff      * * protein2013.ff H C N O S 
fix             QEq all qeq/reax 1 0.0 10.0 1e-6 reaxff
"""
# Specify the output, if needed. Caution: This file gets very lage.
dump_str = """ """
"""
dump            traj all custom 1 traj.lmp id mol type element x y z vx vy vz q 
dump_modify     traj append yes element H C N O S 
"""
# Run sasa (lammps_exe is now optional)
sasa(gro_file, data_file, mol_file, ff_str, dump_str)
# Or with custom LAMMPS executable:
# sasa(gro_file, data_file, mol_file, ff_str, dump_str, lammps_exe)
```
