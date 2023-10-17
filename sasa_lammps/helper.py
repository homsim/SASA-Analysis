import os
import sasa_lammps
import shutil


def check_files(path):
    """
    Check if all the relevant files are present to exec LAMMPS
    and copy LAMMPS input file templates.
    """
    fs = os.listdir(path)

    to_rm = ["etot", "traj.lmp", "thermolog1", "spec.xyz"]
    for f in to_rm:
        if f in fs:
            print(f"{f} file already exists. Will be overwritten...")
            os.remove(os.path.join(path, f))

    # copy the LAMMPS input template to the working dir
    module_dir = os.path.dirname(sasa_lammps.__file__)
    shutil.copy(os.path.join(module_dir, "in.pre"), os.path.join(path, "in.pre"))
    shutil.copy(
        os.path.join(module_dir, "in.template"), os.path.join(path, "in.template")
    )

    return 0


def count_atoms_in_mol(mol_file):
    """Count number of atoms in molecule file"""
    with open(mol_file, "r") as f:
        for line in f:
            if "atoms" in line:
                N = int(line.split()[0])
                break

    return N
