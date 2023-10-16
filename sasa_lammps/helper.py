import os


def check_files(path):
    """Check if all the relevant files are present to exec LAMMPS"""
    fs = os.listdir(path)

    to_rm = ["in.SASA", "etot", "traj.lmp", "thermolog1", "spec.xyz"]
    for f in to_rm:
        if f in fs:
            print(f"{f} file already exists. Will be overwritten...")
            os.remove(os.path.join(path, f))

    """
    prefixes = ["in.", "data."]
    suffixes = [".mol", ".ff"]

    for p in prefixes:
        if [f for f in fs if f.startswith(p)]:
            continue
        else:
            print(f"No {p}* file found. Aborting...")
            return False
    for s in suffixes:
        if [f for f in fs if f.endswith(s)]:
            continue
        else:
            print(f"No *{s} file found. Aborting...")
            return False
    """

    return True


def count_atoms_in_mol(mol):
    """Count number of atoms in molecule file"""
    with open(mol, "r") as f:
        for line in f:
            if "atoms" in line:
                N = int(line.split()[0])
                break

    return N
