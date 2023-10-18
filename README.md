# ToDo

## General

- [X] Create function(s) to write SASA positions
- [x] Write proper `README.md`
- [X] Implement as importable package -> How to partition the methods properly?
- [X] ~~Write/Finalize `setup.py` and/or `requirements.txt`~~
- [X] ~~(Maybe add entry point in `setup.py`)~~
- [X] -> Have to use `conda` because `vmd-python` is not available in PyPi, so a .yml file will have to do to install dependencies
- [ ] Figure out the import...
- [ ] Create a class hirarchy to avoid the countless repeating function arguments

## 1-atomic probe

- [X] Check routine to see if all the files are present
- [X] Write general method for executing LAMMPS

## N_atomic probe

- [X] Rotation of the molecule in order to get defined interactions

## Analysis

- [X] Proper output file format with probe position, energy, residue


# Build and Install

The package requires the `vmd-python` package, which is only distributed in the `conda-forge` channel. For this reason, this package has to be installed in a conda environment as its dependencies cannot be installed from PyPI alone. 
This package will be build locally and then imported in a new conda environment that is created from the `env.yml` file. This is only supposed to be a temporary solution.

## Environment

First a conda environment has to be created from a `env.yml` file. The name of the resulting environment will be `sasa`, but can be changed by modifying the first line in the `env.yml` file:

```
conda env create -f env.yml
conda activate <env-name>    # <env-name>='sasa' or however you named the env
```

## Build

For building the package the created environment should now have all the depedencies.
Make sure that atomatic uploading to the Anaconda repository is turned of (this is the default). Otherwise execute:

```
conda config --set anaconda_upload false
```

In the `SASA-Analysis` directory execute:

```
conda build build_recipe/  
```

This builds the package locally in your `~/anaconda3/conda-bld`. 

## Install

It can then be installed as a package via

```
conda install --use-local sasa_lammps 
```

# Usage

The package really only has one usable method `sasa_lammps.sasa()`:

```
>>> from sasa_lammps import sasa
>>> print(sasa.__doc__)
    Run the SASA analysis on a given macromolecule using a given probe molecule.
    Care must be taken for N-atomic probe molecules: The script does not identify
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
    lammps_exe : str
        Full path to the LAMMPS executable
    params_file : str
        Name of the file with informations on the force field. The template LAMMPS
        input requires:
            pair_style <...>
            pair_coeff <...>
        See the examples/ff_params.dat for an example. As of now the 'units real' command in the in.template kind of restricts to use reaxff pair_style
        (Default: ff_params.dat)
    n_procs : int
        Number of MPI processes to start LAMMPS with (Default: 1)
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
- ff_params.dat

```
from sasa_lammps import sasa

data_file = "data.Cvi_nowater"
mol = "h2o2.mol"
lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi" 

sasa(data_file, mol, lammps_exe)
```
