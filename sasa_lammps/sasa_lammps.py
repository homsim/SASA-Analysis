"""
Package to execute instances of LAMMPS to perform probe analysis of the
solvent-accessible-surface-area (SASA).

...Instead of calling LAMMPS repeatedly one could also write a loop in the 
LAMMPS input instead. This would probably avoid the huge overhead of the LAMMPS
initialization, which probably takes the most amount of computational time right now. 
"""

import os

from sasa_lammps.conversion import convert_data_file
from sasa_lammps.execution import exec_lammps_iterations
from sasa_lammps.conversion import create_sasa_xyz, neighbor_finder


def sasa(data_file, mol_file, lammps_exe, n_procs=1, srad=1.4, samples=100, path="."):
    """
    Run the SASA analysis on a given macromolecule using a given probe molecule.

    Parameters
    ----------
    data_file : str
        Name of the LAMMPS data file of the macromolecule
    mol_file : str
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
    sasa_positions = create_sasa_xyz(path, xyz_file, srad, samples)

    # build neigbor list
    neighbors = neighbor_finder(path, data_file, sasa_positions)

    # execute
    exec_lammps_iterations(
        path, data_file, mol_file, lammps_exe, n_procs, neighbors
    )

    return 0


def main():
    """main function mainly for testing"""
    
    # declare files etc
    data_file = "data.Cvi_nowater"
    mol = "h.mol"
    lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi"

    sasa(data_file, mol, lammps_exe)

    return 0


if __name__ == "__main__":
    main()
