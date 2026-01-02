"""
Fast SASA (Solvent Accessible Surface Area) computation using Monte Carlo sampling.

This module provides a C extension implementation of the SASA algorithm that replaces
the VMD dependency in the sasa_lammps package.
"""

from .sasa_ext import compute_sasa

__all__ = ['compute_sasa']