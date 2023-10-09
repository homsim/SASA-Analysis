"""
Package to execute instances of LAMMPS to perform probe analysis of the
solvent-accessible-surface-area (SASA).
"""

import os
from vmd import atomsel, molecule
from ovito.io import import_file, export_file
import numpy as np
import subprocess

## ignore warnings about ovito being installed via PyPi
#import warnings
#warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

LAMMPS_EXE = '/home/hanna/Programs/lammps-23Jun2022/src/lmp_mpi'

def convert_data_file(path, data_file):
    """Use the Ovito API to convert LAMMPS data file to xyz"""
    pipeline = import_file(f'{path}/{data_file}')

    ### temporary solution
    xyz_file = f"{data_file.split('.')[-1]}.xyz"

    export_file(pipeline, f'{path}/{xyz_file}', 'xyz',
               columns = ['Particle Type', 
                          'Position.X',
                          'Position.Y',
                          'Position.Z',
                         ])
    
    return xyz_file


def create_sasa_xyz(path, xyz_file, srad = 1.4, samples = 100, export_file = 'sasa.xyz'):
    """
    Use an unofficial (?) VMD API to create van der Waals surface points

    Parameters
    ----------
    path : str
        Path to xyz_file and where to export files to
    xyz_file : str
        Name of the xyz-file to use
    srad : float
        Probe radius: Effectively a scaling factor for the vdW radii
        (Default: 1.4, which is the most commonly used because its approx. the
        radius of water)
    samples : int
        Maximum points on the atomic vdW sphere to generate per atom
        (Default: 100)
    export_file : str
        Name of the exported SASA coordinates xyz file
        (Default: sasa.xyz)

    Returns 
    -------
    None

    """
    mol = molecule.load('xyz', f'{path}/{xyz_file}')
    sel = atomsel()

    _, sasa_points = sel.sasa(srad = srad, 
                              samples = samples,
                              points = True)

    export_points = np.array(sasa_points, dtype=str)
    export_points = np.insert(export_points, 0, 'He', axis=1)
    header = f'{len(export_points)}\n '
    np.savetxt(f'{path}/{export_file}', export_points, 
               header=header, comments='', fmt='%s')

    return 0


def count_atoms_in_mol(mol):
    """Count number of atoms in molecule file"""
    with open(mol, 'r') as f:
        for line in f:
            if 'atoms' in line:
                N = int(line.split()[0])
                break

    return N


def check_files(path):
    """Check if all the relevant files are present to exec LAMMPS"""
    fs = os.listdir(path)
    prefixes = ['in.', 'data.']
    suffixes = ['.mol', '.ff']

    if 'in.SASA' in fs:
        print('in.SASA files already exists. Will be overwritten...')

    for p in prefixes:
        if [f for f in fs if f.startswith(p)]: 
            continue
        else:
            print(f'No {p}* file found. Aborting...')
            return False
    for s in suffixes:
        if [f for f in fs if f.endswith(s)]:
            continue
        else:
            print(f'No *{s} file found. Aborting...')
            return False

    return True


def run_lammps(exe, in_file, n_procs = 1):
    """
    Run LAMMPS by running a subprocess. May not be the most elegant way,
    because it cannot handle LAMMPS errors is dependent on OS etc... 
    Also as of now assumes LAMMPS to be build in MPI mode.

    Parameters
    ----------
    exe : str
        Absolute path to the LAMMPS executable
    in_file : str
        Name of the LAMMPS input file
    n_procs : int
        Number of processors to use in MPI execution (Default: 1)

    Returns
    -------
    None

    """
    cmd = f'mpirun -np {n_procs} {exe} -in {in_file}'
    print('Exec: ' + cmd)
    #subprocess.run([cmd], shell = True)

    return 0


def exec_lammps_one(path, sasa_xyz, data_file, in_file, mol_file):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a 1-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    sasa_xyz : str
        Name of the SASA coordinates xyz file to use
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom

    Returns 
    -------
    None

    """
    if check_files(path):
        pass
    else:
        return 0

    sasa_positions = np.genfromtxt(f'{path}/{sasa_xyz}', 
                                   skip_header = 2,
                                   usecols = (1, 2, 3))

    for pos in sasa_positions:
        with open(f'{path}/in.SASA', 'w') as infile:
            with open(f'{path}/{in_file}', 'r') as template:
                for line in template:
                    if 'create_atoms' in line:
                        # insert positions of probe
                        infile.write('create_atoms \t 0 single %s %s %s mol h 39802 \n'
                                    %(pos[0],
                                      pos[1],
                                      pos[2], 
                                     ))
                    elif 'molecule' in line:
                        # insert correct molecule file
                        line.split()
                        infile.write(line.replace('h.mol', mol_file))
                    elif 'read_data' in line:
                        # insert correct data file
                        line.split()
                        infile.write(line.replace('data.Cvi_nowater', 
                                                  data_file))
                    else:
                        infile.write(line)
        
        run_lammps(LAMMPS_EXE, 'in.SASA')
        ###
        break
        ###

    return 0


def exec_lammps_more(path, sasa_xyz, in_file, mol_file, N):
    """
    Execute LAMMPS singlepoints on SASA coordinates using a N-atomic probe

    Parameters
    ----------
    path : str
        Execution path
    sasa_xyz : str
        Name of the SASA coordinates xyz file to use
    in_file : str
        Name of the LAMMPS input file to use as template
    mol_file : str
        Name of the molecule file of the probe atom

    Returns 
    -------
    None

    """

    return 0


def main():
    # declare files etc
    path = '../testing'
    data_file = 'data.Cvi_nowater'
    mol = 'h.mol'

    # convert data file
    xyz_file = convert_data_file(path, data_file)
    # create sasa position file
    create_sasa_xyz(path, xyz_file)

    # count atoms in mol
    N = count_atoms_in_mol(f'{path}/{mol}')
    
    # execute
    if N == 1:
        exec_lammps_one(path, 
                        xyz_file, 
                        data_file, 
                        'in.append', 
                        mol)
    else:
        exec_lammps_more(path,
                        xyz_file, 
                        data_file, 
                        'in.append', 
                        mol,
                        N)

    return 0


if __name__ == '__main__':
    main()
