# ToDo

## General

- [X] Create function(s) to write SASA positions
- [x] Write proper `README.md`
- [X] Implement as importable package -> How to partition the methods properly?
- [X] ~~Write/Finalize `setup.py` and/or `requirements.txt`~~
- [X] ~~(Maybe add entry point in `setup.py`)~~
- [X] -> Have to use `conda` because `vmd-python` is not available in PyPi, so a .yml file will have to do to install dependencies
- [ ] Figure out the import...

## 1-atomic probe

- [X] Check routine to see if all the files are present
- [X] Write general method for executing LAMMPS

## N_atomic probe

- [X] Rotation of the molecule in order to get defined interactions

## Analysis

- [X] Proper output file format with probe position, energy, residue


# Build and Install

The package requires the `vmd-python` package, which is only distributed in the `conda-forge` channel. For this reason, this package has to be installed in a conda environment as its dependencies cannot be installed from PyPI alone. 
The package can be build and then installed locally. 

## Build

For building the package one can use a dedicated conda environment, which has all the depedencies (missing will be installed upon building).

```
conda create --name sasa 
conda activate sasa
```

Make sure that atomatic uploading to the Anaconda repository is turned of (this is the default). Otherwise execute:

```
conda config --set anaconda_upload false
```

In the `SASA-Analysis` directory execute (will take some time):

```
conda build build_recipe/ -c https://conda.ovito.org -c conda-forge 
```

This builds the package locally in your `~/anaconda3/conda-bld`. 

## Install

It can then be installed as a package via

```
conda install --strict-channel-priority --use-local sasaanalysis -c https://conda.ovito.org -c conda-forge
```

Again, this can take some time, because it needs to solve a lot of dependencies.
The order of the `-c` options is important. `https://conda.ovito.org` has to be stated first in order to install `ovito` from this channel and not from the `conda-forge`. For some reason installing `ovito` from the conda-forge leads to it not being found upon importing. If that still does not work (`conda-forge` listed as source for `ovito` in `conda list`) consider first installing `ovito` via 

```
conda install --strict-channel-priority ovito -c https://conda.ovito.org
```

and only after this is done install `sasaanalysis` locally. 
Sadly, as far as I know there is no option to state the channels individually in the `meta.yaml`.

# Usage

...
