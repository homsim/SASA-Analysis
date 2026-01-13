import os
import shutil
from pathlib import Path

from sasa_lammps.constants import (
    ETOT,
    TRAJ,
    THERMOLOG,
    SPEC,
    IN_PRE,
    IN_TEMPLATE
)
from sasa_lammps.resource_loader import resources


def _check_files(path: str) -> None:
    """
    Check if all the relevant files are present to exec LAMMPS
    and copy LAMMPS input file templates.
    """
    fs = os.listdir(path)

    to_rm = [ETOT, TRAJ, THERMOLOG, SPEC]
    for f in to_rm:
        if f in fs:
            print(f"{f} file already exists. Will be overwritten...")
            os.remove(Path(path) / f)

    # copy the LAMMPS input template to the working dir
    shutil.copy(resources / IN_PRE, Path(path) / IN_PRE)
    shutil.copy(resources / IN_TEMPLATE, Path(path) / IN_TEMPLATE)

    return 0


def _count_atoms_in_macromol(data_file: str) -> int:
    """Count number of atoms in macro molecule file"""
    with open(data_file, "r") as f:
        for line in f:
            if "atoms" in line:
                atom_number = int(line.split()[0])
                break

    return atom_number

def _count_atoms_in_mol(mol_file: str) -> int:
    """Count number of atoms in molecule file"""
    with open(mol_file, "r") as f:
        for line in f:
            if "atoms" in line:
                N = int(line.split()[0])
                break

    return N


def _write_params_file(string: str, file_name: str) -> None:
    """Write string to a file. Needed to create the include files for LAMMPS to read"""
    with open(file_name, "w") as f:
        f.write(f"{string}\n")

    return 0


def _read_last_two(path: str, file_name: str) -> list[float, float]:
    """Read last to lines of file, convert to float and return the read values"""
    with open(Path(path) / file_name, "r") as f:
        lines = f.readlines()[-2:]
    l1 = float(lines[0])
    l2 = float(lines[1])

    return l1, l2
