import os
import subprocess
import numpy as np

from sasa_lammps.helper import check_files
from sasa_lammps.conversion import rotate_probe


def exec_lammps_one(path, data_file, in_file, mol_file, lammps_exe, n_procs):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a 1-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom
    lammps_exe : str
        Full path to the LAMMPS executable
    n_procs : int
        Number of MPI processes to start LAMMPS with (Default: 1)


    Returns
    -------
    None

    """
    # check if all the files are there...
    # basically not necessary because of exception handling in run_lmp()
    if check_files(path):
        pass
    else:
        return 0

    # create the sasa positions
    sasa_positions = np.genfromtxt(
        os.path.join(path, "sasa.xyz"), skip_header=2, usecols=(1, 2, 3)
    )

    for i, pos in enumerate(sasa_positions):
        with open(os.path.join(path, "in.SASA"), "w") as infile:
            with open(os.path.join(path, in_file), "r") as template:
                # overwrite probe mol positions etc
                for line in template:
                    if "create_atoms" in line:
                        # insert positions of probe
                        infile.write(
                            f"create_atoms \t 0 single {pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f} mol h 39802 \n"
                        )
                    elif "molecule" in line:
                        # insert correct molecule file
                        infile.write(line.replace("h.mol", mol_file))
                    elif "read_data" in line:
                        # insert correct data file
                        infile.write(line.replace("data.Cvi_nowater", data_file))
                    else:
                        infile.write(line)

        # execute LAMMPS
        run_lmp(lammps_exe, "in.SASA", i, n_procs=n_procs)

    return 0


def exec_lammps_more(
    path, data_file, in_file, mol_file, lammps_exe, n_procs, neighbors
):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a N-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom
    lammps_exe : str
        Full path to the LAMMPS executable
    n_procs : int
        Number of MPI processes to start LAMMPS with (Default: 1)

    Returns
    -------
    None

    """
    # check if all the files are there...
    # basically not necessary because of exception handling in run_lmp()
    if check_files(path):
        pass
    else:
        return 0

    # create the sasa positions
    sasa_positions = np.genfromtxt(
        os.path.join(path, "sasa.xyz"), skip_header=2, usecols=(1, 2, 3)
    )

    # rotate the probe molecule
    rot = rotate_probe(path, data_file, sasa_positions, neighbors)

    for i, pos in enumerate(sasa_positions):
        with open(os.path.join(path, "in.SASA"), "w") as infile:
            with open(os.path.join(path, in_file), "r") as template:
                # overwrite probe mol positions etc
                for line in template:
                    if "create_atoms" in line:
                        # insert positions of probe and rotate
                        infile.write(
                            f"""create_atoms \t 0 single {pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f} mol h 39802 &
                            \t rotate {rot[i][0]:.3f} {rot[i][1]:.3f} {rot[i][2]:.3f} {rot[i][3]:.3f}
                            """
                        )
                    elif "molecule" in line:
                        # insert correct molecule file
                        infile.write(line.replace("h.mol", mol_file))
                    elif "read_data" in line:
                        # insert correct data file
                        infile.write(line.replace("data.Cvi_nowater", data_file))
                    else:
                        infile.write(line)

        # execute LAMMPS
        run_lmp(lammps_exe, "in.SASA", i, n_procs=n_procs)
        """
        # instead of overwriting the lammps in file, I could just call lammps
        # with -var , i.e. parse the variables to lammps directly. this just
        # requires a modified lammps in file...
        run_args = [lammps_exe,
                    "in.SASA",
                    i,
                    n_procs=n_procs,
                    sasa_pos,
                    sasa_rot_vec,
                    sasa_rot_ang
                    ]
        run_lmp(*run_args)
        """

    return 0


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

    ### temporary placeholders
    data_file = "data.bla"
    mol_file = "bla.mol"
    sasa_pos = [0.0, 0.0, 0.0]
    sasa_rot_vec = [0.0, 0.0, 0.0]
    sasa_rot_ang = 0.0

    new_cmd = f"""
    mpirun -np {n_procs} {lammps_exe} -in {in_file} \
        -var DataFile {data_file} -var MolFile {mol_file} \
        -var sasaX {sasa_pos[0]:.3f} -var sasaY {sasa_pos[1]:.3f} \
        -var sasaZ {sasa_pos[2]:.3f} -var rotAng {sasa_rot_ang:.3f} \
        -var rotVecX {sasa_rot_vec[0]:.3f} -var rotVecY {sasa_rot_vec[1]:.3f} \
        -var rotVecZ {sasa_rot_vec[2]:.3f} 
    """
    ###

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
    except subprocess.CalledProcessError as exc:
        print(exc)
        raise
    finally:
        finish_str = f" Finished iteration {iterat}... "
        print("{s:{c}^{n}}".format(s="", n=60, c="#"))
        print("{s:{c}^{n}}".format(s=finish_str, n=60, c="#"))
        print("{s:{c}^{n}}".format(s="", n=60, c="#"))

    return 0
