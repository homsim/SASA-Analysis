#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

from scipy import stats
from scipy.signal import savgol_filter
import argparse as ag
import subprocess



parser = ag.ArgumentParser(prog='Run SASA for several in files')
parser.add_argument('-f','--file',
                    type=str,
                    default='h.mol',
                    help='mol file for current SASA analysis'
                    )
args = parser.parse_args()

mol_file=args.file

pos=np.loadtxt('sasa.xyz',skiprows=2, usecols=(1,2,3))

# wirte a lot of new input files for each new position
for i in range(len(pos)):
    with open('in.SASA', 'w') as nfile:
        with open('../in.append', 'r') as file:
            for line in file:
                if 'create_atoms' in line:
                        nfile.write('create_atoms \t 0 single %s %s %s mol h 39802 \n'
                                    %(pos[i][0],
                                      pos[i][1],
                                      pos[i][2], 
                                     ))
                elif 'molecule' in line:
                    line.split()
                    nfile.write(line.replace('h.mol', mol_file))
                else:
                    nfile.write(line)
        nfile.close()
        
    for files in os.walk('.', topdown=True):
        if 'in.' and 'data.' and '.mol' and '.ff':
            subprocess.run(['mpirun -np 1 /home/hanna/Programs/lammps-23Jun2022/src/lmp_mpi < in.SASA'], shell=True)

