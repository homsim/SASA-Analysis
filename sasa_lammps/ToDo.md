# ToDo

## General

- [X] Create function(s) to write SASA positions
- [x] Write proper `README.md`
- [X] Implement as importable package -> How to partition the methods properly?
- [X] ~~Write/Finalize `setup.py` and/or `requirements.txt`~~
- [X] ~~(Maybe add entry point in `setup.py`)~~
- [X] -> Have to use `conda` because `vmd-python` is not available in PyPi, so a .yml file will have to do to install dependencies
- [X] Create a class hirarchy to avoid the countless repeating function arguments
- [X] Enable 'real' multiprocessing by running multiple LAMMPS instances in parallel
- [X] Strings for the force field and dump commands
- [ ] Function argument to initialize the calculation with out running LAMMPS
- [X] Add a proper termination criterion to kill all running processes
- [ ] Redirect error outputs of LAMMPS (`capture_output`)
- [X] Add star-import for `sasa()`
- [ ] Run simple benchmark on the MP performance

## 1-atomic probe

- [X] Check routine to see if all the files are present
- [X] Write general method for executing LAMMPS

## N_atomic probe

- [X] Rotation of the molecule in order to get defined interactions

## Analysis

- [X] Proper output file format with probe position, energy, residue

