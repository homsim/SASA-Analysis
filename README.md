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
    Run the SASA analysis on a given macromolecule using a given probe molecule.
    Care must be taken for N-atomic probe molecules: The package does not identify
    a plane or something in the probe molecule. It just makes sure that at every
    interaction site the probe faces the macromolecule with the same orientation.
    However, the orientation itself is purely determined by the configuration
    given in the mol_file.

    Parameters
    ----------
    data_file : str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
        Name of the LAMMPS mol file to use as probe of the SAS
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

In the directory `example` you find an example on the CviUPO with H2O2 which samples 12848 points on the SAS. The following files need to be provided for LAMMPS:
- data.Cvi_nowater
- h2o2.mol
- protein2013.ff

```python
from sasa_lammps import sasa

data_file = "data.Cvi_nowater"
mol = "h2o2.mol"
lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi" 
ff_str = """
pair_style      reaxff NULL checkqeq no safezone 1.6 mincap 100 minhbonds 150
pair_coeff      * * protein2013.ff H C N O S X X Cl  
"""
dump_str = """
dump            traj all custom 1 traj.lmp id mol type element x y z vx vy vz q 
dump_modify     traj append yes element H C N O S Mg Fe Cl
"""


sasa(data_file, mol_file, ff_str, dump_str, lammps_exe, n_procs = 2)
```
