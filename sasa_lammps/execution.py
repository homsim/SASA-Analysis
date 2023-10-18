import os
import subprocess
import numpy as np

from sasa_lammps.helper import check_files, count_atoms_in_mol, write_params_file
from sasa_lammps.conversion import rotate_probe


def exec_lammps_iterations(
    path, data_file, mol_file, ff_str, dump_str, lammps_exe, n_procs, neighbors
):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a N-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
        Name of the molecule file of the probe atom
    ff_str : str
        Force field parameters to provide to LAMMPS
    dump_str : str
        Dump command to provide to LAMMPS
    lammps_exe : str
        Full path to the LAMMPS executable
    n_procs : int
        Number of MPI processes to start LAMMPS with (Default: 1)
    neighbors : dict
        Dictionary of neighbor list informations.
        Output of the conversion.neighbor_finder() method

    Returns
    -------
    None

    """

    # remove existing files and copy input templates...
    check_files(path)

    # write ff_str and dump_str to files for LAMMPS to read in
    write_params_file(ff_str, "ff_params.dat")
    write_params_file(dump_str, "dump_com.dat")

    # get the energies for the isolated macro- and probe molecule, respectively
    e_mol, e_prob = pre_calc(
        path, lammps_exe, data_file, mol_file, n_procs
    )

    # create the sasa positions
    sasa_positions = np.genfromtxt(
        os.path.join(path, "sasa.xyz"), skip_header=2, usecols=(1, 2, 3)
    )
    n_probes = len(sasa_positions)

    # create final output file header and write to spec.xyz
    header = f"{n_probes}\natom\tx\ty\tz\tres\tetot [eV]\teint [eV]\n"
    with open(os.path.join(path, "spec.xyz"), "w") as f:
        f.write(header)

    # rotate the probe molecule for n-atomic probes (n > 1)
    if count_atoms_in_mol(os.path.join(path, mol_file)) > 1:
        rotations = rotate_probe(path, data_file, sasa_positions, neighbors)
    else:
        rotations = np.zeros((n_probes, 4))
        # add some direction otherwise LAMMPS raises an error because of the zero vector
        rotations[:, 1] += 1.000

    iterators = [sasa_positions, neighbors["res"], rotations]
    for i, (pos, res, rot) in enumerate(zip(*iterators)):
        # execute LAMMPS
        run_args = [
            lammps_exe,
            "in.template",
            i,
            n_probes,
            data_file,
            mol_file,
            pos,
            rot,
            n_procs,
        ]
        run_lmp(*run_args)

        # get the current total energy
        with open(os.path.join(path, "etot"), "r") as f:
            for line in f:
                pass
            etot = float(line)
        eint = etot - (e_mol + e_prob)

        # convert kcal/mole to eV. This is correct as long as LAMMPS
        # uses the units real command
        etot = etot * 0.04336
        eint = eint * 0.04336

        # append the final output file
        with open(os.path.join(path, "spec.xyz"), "a") as f:
            f.write("He\t")  # "He" is only a dummy
            f.write(f"{pos[0]:.3f}\t")
            f.write(f"{pos[1]:.3f}\t")
            f.write(f"{pos[2]:.3f}\t")
            f.write(f"{res}\t")
            f.write(f"{etot:.3f}\t")
            f.write(f"{eint:.3f}\n")

    return 0


def pre_calc(path, lammps_exe, data_file, mol_file, n_procs):
    """
    Do two pre-runs in LAMMPS: One of the isolated macromolecule and one of the isolated probe molecule

    Parameters
    ----------
    path : str
        Execution path
    lammps_exe : str
        Absolute path to the LAMMPS executable
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
        Name of the molecule file of the probe atom
    n_procs : int
        Number of MPI processes to start LAMMPS with (Default: 1)


    Returns
    -------
    e_mol : float
        Energy of the macro molecule in kcal/mole
    e_prob : float
        Energy of the probe molecule in kcal/mole

    """

    run_lmp(
        lammps_exe,
        "in.pre",
        0,
        0,
        data_file,
        mol_file,
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        n_procs,
    )

    # the unpacking here is crap. This must be more elegant somehow...
    with open(os.path.join(path, "etot"), "r") as f:
        e_mol, e_prob = zip(f.readlines()[1:])

    return float(e_mol[0]), float(e_prob[0])


def run_lmp(
    lammps_exe,
    in_file,
    iterat,
    max_iterat,
    data_file,
    mol_file,
    pos,
    rot,
    n_procs,
):
    """
    Run LAMMPS by running a subprocess. May not be the most elegant way,
    because it cannot handle LAMMPS errors and is dependent on OS etc...
    Also as of now assumes LAMMPS to be build in MPI mode.

    Parameters
    ----------
    lammps_exe : str
        Absolute path to the LAMMPS executable
    in_file : str
        Name of the LAMMPS input file
    iterat : int
        Number of iteration
    max_iterat : int
        Max Number of iterations
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
        Name of the molecule file of the probe atom
    pos : list
        x, y, z position list of the SAS positions
    rot : list
        List of rotation data:
        rot[0]: Rotation angle
        rot[1]: X-component of rotation vector
        rot[2]: Y-component of rotation vector
        rot[3]: Z-component of rotation vector
    n_procs : int
        Number of processors to use in MPI execution

    Returns
    -------
    None

    """

    capture_output = True  # Whether to capute LAMMPS output to stdout or not

    cmd = f"""
    mpirun -np {n_procs} {lammps_exe} -in {in_file} \
        -var DataFile {data_file} -var MolFile {mol_file} \
        -var sasaX {pos[0]:.3f} -var sasaY {pos[1]:.3f} \
        -var sasaZ {pos[2]:.3f} -var rotAng {rot[0]:.3f} \
        -var rotVecX {rot[1]:.3f} -var rotVecY {rot[2]:.3f} \
        -var rotVecZ {rot[3]:.3f}
    """

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
        lead_spaces = " " * (len(str(max_iterat)) - len(str(iterat)))
        finish_str = f" Finished iteration |{lead_spaces}{iterat}/{max_iterat}| "
        print("{s:{c}^{n}}".format(s="", n=70, c="#"))
        print("{s:{c}^{n}}".format(s=finish_str, n=70, c="#"))
        print("{s:{c}^{n}}".format(s="", n=70, c="#"))

    return 0
