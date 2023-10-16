import os
import numpy as np

from sasa_lammps.conversion import neighbor_finder


def out_analysis(
    path, sasa_positions, neighbors, e_mol=-405549.69, e_prob=-263.49651
):
    """

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    data_file : str
        Name of the LAMMPS data file of the macromolecule
    e_mol: float
        Energy of the isolated macromolecule in kcal/mole
    e_prob: float
        Energy of the isolated probe molecule in kcal/mole

    Returns
    -------
    None

    """
    data_out = {
        "atom": ["He"] * len(sasa_positions),
        "x": [],
        "y": [],
        "z": [],
        "res": [],
        "etot": [],
        "eint": [],
    }

    data_out["x"] = sasa_positions[:, 0]
    data_out["y"] = sasa_positions[:, 1]
    data_out["z"] = sasa_positions[:, 2]

    data_out["etot"] = np.genfromtxt(os.path.join(path, "etot"))
    data_out["eint"] = data_out["etot"] - (e_mol + e_prob)

    data_out["res"] = neighbors["res"]

    header = f"{len(data_out['res'])}\natom\tx\ty\tz\tres\tetot [kcal/mole]\teint [kcal/mole]"
    np.savetxt(
        os.path.join(path, "spec.xyz"),
        np.array([data_out[d] for d in data_out.keys()]).T,
        fmt="%s",
        delimiter="\t",
        comments="",
        header=header,
    )

    return 0
