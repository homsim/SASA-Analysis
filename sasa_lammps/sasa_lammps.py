"""
Package to execute instances of LAMMPS to perform probe analysis of the
solvent-accessible-surface-area (SASA).
"""

import os

from sasa_lammps.conversion import convert_data_file
from sasa_lammps.execution import exec_lammps_one, exec_lammps_more
from sasa_lammps.helper import count_atoms_in_mol
from sasa_lammps.conversion import create_sasa_xyz, neighbor_finder
from sasa_lammps.out_analysis import out_analysis


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
    sasa_positions = create_sasa_xyz(path, xyz_file, srad=srad, samples=samples)

    # count atoms in mol
    num_atoms = count_atoms_in_mol(os.path.join(path, mol))

    # build neigbor list
    neighbors = neighbor_finder(path, data_file, sasa_positions)

    # execute
    if num_atoms == 1:
        exec_lammps_one(path, data_file, "in.append", mol, lammps_exe, n_procs)
    else:
        exec_lammps_more(
            path, data_file, "in.append", mol, lammps_exe, n_procs, neighbors
        )

    out_analysis(path, sasa_positions, neighbors)

    return 0


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
