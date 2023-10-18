from sasa_lammps import sasa

data_file = "data.Cvi_nowater"
mol = "h2o2.mol"
lammps_exe = "/opt/lammps-23Jun2022/src/lmp_mpi" 

ff_str = """
pair_style      reaxff NULL checkqeq no safezone 1.6 mincap 100 minhbonds 150
pair_coeff      * * protein2013.ff H C N O S X X Cl  
"""
dump_str = """
dump            traj all custom 1 traj.lmp id mol type element x y z vx vy vz q 
dump_modify     traj append yes element H C N O S Mg Fe Cl
"""


sasa(data_file, mol_file, ff_str, dump_str, lammps_exe)
