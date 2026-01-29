.. SASA-Lammps documentation master file, created by
   sphinx-quickstart on Thu Jan 29 17:29:15 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SASA-Lammps documentation
=========================

**SASA-Lammps** calculates the solvent-accesible surface analysis (SASA) on a given macro-molecule (protein). It places a given probe-molecule on the surface points to compute a 3D interaction energy surface. The potential energy is calculated with a ReaxFF potential that needs to be provided (see below) and uses `LAMMPS <https://www.lammps.org>`__ as actual implementation for the energy calculations.

This package was written for the publication `J. Phys. Chem. B 2025, 129, 44, 11374–11386 by Poggemann et al. <https://pubs.acs.org/doi/10.1021/acs.jpcb.5c03518>`__. The original version of the package relied on the VMD molecular visualization program for the calculation of the solvent-accesible surface, this dependency was removed in a major updated and repaced by a custom implementation. 


.. toctree::
   installation
   api
   :maxdepth: 2
   :caption: Contents:

