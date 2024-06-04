# Build and Install

The package requires the `vmd-python` package, which is only distributed in the `conda-forge` channel. For this reason, this package has to be installed in a conda environment as its dependencies cannot be installed from PyPI alone. 
This package will be build locally and then imported in a new conda environment that is created from the `env.yml` file. This is only supposed to be a temporary solution.

## Environment

First a conda environment has to be created from a `env.yml` file. The name of the resulting environment will be `sasa`, but can be changed by modifying the first line in the `env.yml` file:

```bash
conda env create -f env.yml
```

```bash
conda activate <env-name>    # <env-name>='sasa' or however you named the env
```

## Build

For building the package the created environment should now have all the
depedencies.
Make sure that atomatic uploading to the Anaconda repository is turned of (this is the default). Otherwise execute:

```bash
conda config --set anaconda_upload false
```

In the `SASA-Analysis` directory execute:

```bash
conda build build_recipe/  
```

This builds the package locally in your `~/anaconda3/conda-bld`. 

## Install

It can then be installed as a package via

```bash
conda install --use-local sasa_lammps_multi 
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
    lammps_exe : str
        Full path to the LAMMPS executable
    n_procs : int
        Number of LAMMPS instances to run in parallel (Default: 1)
    srad : float
        Probe radius: Effectively a scaling factor for the vdW radii
        (Default: 1.4, which is the most commonly used because its approx. the
        radius of water)
    samples : int
        Maximum points on the atomic vdW sphere to generate per atom (Default: 100)
    path : str
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
# Path to you lammps executable
lammps_exe =  "/opt/lammps-23Jun2022/src/lmp_mpi"
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
# Run sasa
sasa(gro_file, data_file, mol_file, ff_str, dump_str, lammps_exe)
```
