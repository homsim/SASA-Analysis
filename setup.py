"""SASA-Analysis

Repeatedly executes LAMMPS instances to probe the solvent-accessible surface-area
of a given macromolecule in order to create an energy landscape of that molecule.
"""


DOCSTRING = __doc__.split("\n")

from setuptools import setup, find_packages

requirements = [
]

setup(
    name="sasa_lammps_multi",
    version=0.1,
    install_requires=requirements,
    include_package_data=True,
    packages=find_packages(),
)
