"""
Package to execute instances of LAMMPS to perform probe analysis of the
solvent-accessible-surface-area (SASA).
"""

import os
import subprocess
import numpy as np
from vmd import atomsel, molecule
from ovito.io import import_file, export_file


def convert_data_file(path, data_file):
    """Use the Ovito API to convert LAMMPS data file to xyz"""
    pipeline = import_file(f"{path}/{data_file}")

    ### temporary solution
    xyz_file = f"{data_file.split('.')[-1]}.xyz"

    export_file(
        pipeline,
        os.path.join(path, xyz_file),
        "xyz",
        columns=[
            "Particle Type",
            "Position.X",
            "Position.Y",
            "Position.Z",
        ],
    )

    return xyz_file


def create_sasa_xyz(path, xyz_file, srad=1.4, samples=100):
    """
    Use an unofficial (?) VMD API to create van der Waals surface points

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    xyz_file : str
        Name of the xyz-file to use
    srad : float
        Probe radius: Effectively a scaling factor for the vdW radii
        (Default: 1.4, which is the most commonly used because its approx. the
        radius of water)
    samples : int
        Maximum points on the atomic vdW sphere to generate per atom
        (Default: 100)

    Returns
    -------
    None

    """
    export_file = "sasa.xyz"
    mol = molecule.load("xyz", os.path.join(path, xyz_file))
    sel = atomsel()

    _, sasa_points = sel.sasa(srad=srad, samples=samples, points=True)

    export_points = np.array(sasa_points, dtype=str)
    export_points = np.insert(export_points, 0, "He", axis=1)

    header = f"{len(export_points)}\n "
    np.savetxt(
        os.path.join(path, export_file),
        export_points,
        header=header,
        comments="",
        fmt="%s",
    )

    return 0


def count_atoms_in_mol(mol):
    """Count number of atoms in molecule file"""
    with open(mol, "r") as f:
        for line in f:
            if "atoms" in line:
                N = int(line.split()[0])
                break

    return N


def check_files(path):
    """Check if all the relevant files are present to exec LAMMPS"""
    fs = os.listdir(path)
    prefixes = ["in.", "data."]
    suffixes = [".mol", ".ff"]

    if "in.SASA" in fs:
        print("in.SASA file already exists. Will be overwritten...")

    if "etot" in fs:
        print("etot file already exists. Will be overwritten...")
        os.remove(os.path.join(path, "etot"))

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

    return True


def run_lmp(lammps_exe, in_file, iterat, n_procs=1, capture_output=True):
    """
    Run LAMMPS by running a subprocess. May not be the most elegant way,
    because it cannot handle LAMMPS errors and is dependent on OS etc...
    Also as of now assumes LAMMPS to be build in MPI mode.

    Parameters
    ----------
    exe : str
        Absolute path to the LAMMPS executable
    in_file : str
        Name of the LAMMPS input file
    iterat : int
        Number of iteration
    n_procs : int
        Number of processors to use in MPI execution (Default: 1)
    capture_output : bool
        Whether to print the LAMMPS output to stdout or not (Default: True)

    Returns
    -------
    None

    """
    cmd = f"mpirun -np {n_procs} {lammps_exe} -in {in_file}"
    # somehow the vmd-python module messes up $TMPDIR which is why I need
    # to reset the env variables in the subprocess
    try:
        subprocess.run(
            [cmd],
            shell=True,
            env=os.environ,
            check=True,
            capture_output=not capture_output,
        )
    except subprocess.CalledProcessError as e:
        print(e)
    finally:
        print(f"Finished iteration {iterat}...")

    return 0


def exec_lammps_one(path, data_file, in_file, mol_file, lammps_exe, n_procs):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a 1-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom

    Returns
    -------
    None

    """
    if check_files(path):
        pass
    else:
        return 0

    sasa_positions = np.genfromtxt(
        os.path.join(path, "sasa.xyz"), skip_header=2, usecols=(1, 2, 3)
    )

    for i, pos in enumerate(sasa_positions):
        with open(os.path.join(path, "in.SASA"), "w") as infile:
            with open(os.path.join(path, in_file), "r") as template:
                for line in template:
                    if "create_atoms" in line:
                        # insert positions of probe
                        infile.write(
                            f"create_atoms \t 0 single {pos[0]} {pos[1]} {pos[2]} mol h 39802 \n"
                        )
                    elif "molecule" in line:
                        # insert correct molecule file
                        infile.write(line.replace("h.mol", mol_file))
                    elif "read_data" in line:
                        # insert correct data file
                        infile.write(line.replace("data.Cvi_nowater", data_file))
                    else:
                        infile.write(line)

        run_lmp(lammps_exe, "in.SASA", i, n_procs=n_procs)

    return 0


def exec_lammps_more(path, sasa_xyz, in_file, mol_file, num_atoms):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a N-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    sasa_xyz : str
        Name of the SASA coordinates xyz file to use
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom
    num_atoms : int
        Number of atoms in the probe molecule

    Returns
    -------
    None

    """

    return 0


def sasa(data_file, mol, lammps_exe, n_procs=1, srad=1.4, samples=100, path="."):
    """
    Run the SASA analysis on a given macromolecule using a given probe molecule.

    Parameters
    ----------
    data_file : str
        Name of the LAMMPS data file of the macromolecule
    mol : str
        Name of the LAMMPS mol file to use as probe of the SAS
    lammps_exe : str
        Full path to the LAMMPS executable
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

    """
    # convert data file
    xyz_file = convert_data_file(path, data_file)
    # create sasa position file
    create_sasa_xyz(path, xyz_file, srad=srad, samples=samples)

    # count atoms in mol
    num_atoms = count_atoms_in_mol(os.path.join(path, mol))

    # execute
    if num_atoms == 1:
        exec_lammps_one(path, data_file, "in.append", mol, lammps_exe, n_procs=n_procs)
    else:
        exec_lammps_more(path, data_file, "in.append", mol, num_atoms)


def main():
    """
    main function mainly for testing
    """
    # declare files etc

    data_file = "data.Cvi_nowater"
    mol = "h.mol"
    lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi"

    sasa(data_file, mol, lammps_exe)

    return 0


if __name__ == "__main__":
    main()
