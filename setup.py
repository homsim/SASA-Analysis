"""SASA-Analysis

Perform solvent-accessible surface area (SASA) analysis of a given 
macromolecule by probing the SAS with a probe molecule. 
Uses LAMMPS singlepoint calculations for that.
"""

DOCLINES = __docs__.split('\n')

from distutils.core import setup

setup(name = 'sasaanalysis',
      version = 0.1,
      description = DOCLINES,
      author = '',
      author_email = '',
      url = 'https://github.com/hpoggemann/SASA-Analysis',
      packages = ['src'],
      install_requires = [
          'vmd-python',
          'ovito',
          'numpy',
      ],
      tests_require = ['pytest'],
     )
