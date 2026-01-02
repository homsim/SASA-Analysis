import os
import numpy as np
from ovito.io import import_file, export_file
from ovito.data import NearestNeighborFinder

from sasa_lammps.constants import SASAXYZ
from sasa_lammps.sasa_core import _create_sasa_xyz as _create_sasa_xyz_impl


def _convert_data_file(path, data_file):
    """Use the Ovito API to convert LAMMPS data file to xyz"""
    pipeline = import_file(os.path.join(path, data_file))

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


def _create_sasa_xyz(path, xyz_file, srad, samples):
    """
    Create van der Waals surface points using efficient C extension.

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    xyz_file : str
        Name of the xyz-file to use
    srad : float
        Probe radius: Effectively a scaling factor for the vdW radii
    samples : int
        Maximum points on the atomic vdW sphere to generate per atom

    Returns
    -------
    sasa_points : numpy.ndarray
        (N, 3) Array of coordinates on the SAS.
        N is loosely determined by 'samples' argument.

    """


    ## is this wrapper really needed????
    return _create_sasa_xyz_impl(path, xyz_file, srad, samples)


def _neighbor_finder(path, data_file, sasa_positions):
    """
    Compute informations on the nearest neighbors of every SAS point, i.e. which
    atom of the macromolecule is closest to the SAS point.

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    sasa_positions : numpy.array
        Coordinates of the points on the SAS

    Returns
    -------
    neighbors: dict{"id", "pos", "res", "dist"}
        Dictionary of informations on nearest neighbors:
        neighbors["id"]: Particle ID of the nearest neighbor
        neighbors["pos"]: Position of the nearest neighbor
        neighbors["res"]: Molecule ID of the nearest neighbor, which can be used for residue identification
        neighbors["dist"]: Distance of the nearest neighbor

    """

    pipeline = import_file(os.path.join(path, data_file))
    data = pipeline.compute()

    finder = NearestNeighborFinder(1, data)  # 1st nearest neighbor
    neighbors = {"id": [], "pos": [], "res": [], "dist": []}
    for pos in sasa_positions:
        for neigh in finder.find_at(pos):
            neighbors["id"].append(data.particles["Particle Identifier"][neigh.index])
            neighbors["pos"].append(data.particles["Position"][neigh.index])
            neighbors["res"].append(data.particles["Molecule Identifier"][neigh.index])
            neighbors["dist"].append(neigh.distance)

    return neighbors


def _rotate_probe(path, data_file, sasa_positions, neighbors):
    """
    Rotate the probe molecule on the SAS. Finds the nearest SAS point for each atom
    of the macromolecule. The rotation is then done with respect to the center-of-mass
    of the macromolecule.

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    data_file: str
        Name of the LAMMPS data file of the macromolecule
    sasa_positions : numpy.array
        Coordinates of the points on the SAS
    neighbors : dict
        Dictionary of neighbor list informations.
        Output of the neighbor_finder() method

    Returns
    -------
    rot_export: numpy.ndarray
        Array of shape (N, 4) where N is the number of points on the SAS.
        rot_export[N, 0]: rotation angle
        rot_export[N, 1:]: x, y, z coordinates of the respective rotation vector

    """

    pipeline = import_file(os.path.join(path, data_file))

    # calculate the center-of-mass of the macromolecule
    data = pipeline.compute()
    data_positions = data.particles.positions
    data_masses = data.particles.masses
    com = np.sum([m * p for m, p in zip(data_masses, data_positions)], axis=0) / np.sum(
        data_masses
    )

    # get the indices of the particle container according to the IDs
    ordering = np.argsort(data.particles.identifiers)
    indices = ordering[
        np.searchsorted(data.particles.identifiers, neighbors["id"], sorter=ordering)
    ]

    # calculate rotation angles and vectors to parse them to LAMMPS
    rot_export = []
    for index, sasa_pos in zip(indices, sasa_positions):
        curr_id = neighbors["id"][index]
        com_vec = data_positions[curr_id] - com
        sas_vec = sasa_pos - data_positions[curr_id]

        angle = np.rad2deg(
            np.arccos(
                np.dot(com_vec, sas_vec)
                / (np.linalg.norm(com_vec) * np.linalg.norm(sas_vec))
            )
        )

        vector = np.cross(com_vec, sas_vec)
        rot_export.append(np.insert(vector, 0, angle))

    # writing the file is not necessary
    """
    header = f"{np.shape(rot_export)[0]}\n "
    np.savetxt(
        os.path.join(path, "rotations.dat"),
        rot_export,
        header=header,
        comments="",
        fmt="%s",
    )
    """

    return rot_export
