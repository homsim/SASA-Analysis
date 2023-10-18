from sasa_lammps import sasa

data_file = "data.Cvi_nowater"
mol = "h2o2.mol"
lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi" 

sasa(data_file, mol, lammps_exe)