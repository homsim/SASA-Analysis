"""SASA-Analysis

Repeatedly executes LAMMPS instances to probe the solvent-accessible surface-area
of a given macromolecule in order to create an energy landscape of that molecule.
"""


DOCSTRING = __doc__split('\n')

from distutils.core import setup

requirements = [
    'numpy',
#    'ovito',
#    'vmd-python'
]

setup(
    name = 'sasaanalysis',
    version = 0.1,
    install_requires = requirements,
)
