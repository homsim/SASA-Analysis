#!/usr/bin/env python
# coding: utf-8

# In[58]:


import numpy as np
import os
import pandas as pd
from helperfunctions import *

from ovito.io import *
from ovito.modifiers import *
from ovito.vis import Viewport
from ovito.data import *
from ovito.pipeline import *


# In[59]:


e_mol=-562980.36365695# # Enyzme energy in kcal/mol
e_species =0 # species in kcal/mol
path= r'/home/hanna/simulations/MM/lammps/AaeUPO/04_plasma/SASA/'
folder='O_ion/'
file='sasa.xyz'


# In[60]:


data={'atom':[], 'x':[], 'y':[], 'z':[],'res':[], 'etot':[], 'e_int':[]}
# load species and coords into dict
with open(path+file) as f:
    for i,l in enumerate(f):
        if i>1:
            data['atom'].append('H')
            data['x'].append((l.split()[1]))
            data['y'].append((l.split()[2]))
            data['z'].append((l.split()[3]))
# load final energies into dict
d=np.loadtxt(path+folder+'etot')
for i in range(len(d)):
    data['etot'].append(d[i])
# calc. interaction energy    
for i in range(len(data['etot'])):
    data['e_int'].append(data['etot'][i]-(e_mol+e_species))


# In[61]:


# Create NearestNeighbors list
# load positions from sasa file
pos = {'coords':[]}
with open(path+file) as f:
    for i,l in enumerate(f):
        if i>1:
            floats = [float(x) for x in l.split()[1:]]
            coords=np.asarray(floats)
            pos['coords'].append(coords)

neighbor ={'id':[],'pos':[],'res':[],'dist':[]}
# Load input simulation file.
pipeline = import_file(path+'data.Aae_nowater',atom_style = 'full')
d = pipeline.compute()
# Initialize neighbor finder object.
# Visit the 12 nearest neighbors of each particle.
N = 1
finder = NearestNeighborFinder(N, d)
# Visit particles closest to some spatial point (x,y,z):
for i in range(len(pos['coords'])):
    for neigh in finder.find_at(pos['coords'][i]):
        neighbor['id'].append(d.particles['Particle Identifier'][neigh.index])
        neighbor['pos'].append(d.particles['Position'][neigh.index])
        neighbor['res'].append(d.particles['Molecule Identifier'][neigh.index])
        neighbor['dist'].append(neigh.distance)
        
with open(path+'NearestNeighbors.txt', 'w') as fp:
    fp.write('id_surf_atom \t Residue \n')
    for i in range(len(neighbor['id'])):
        fp.write('%s ' %neighbor['id'][i])
        fp.write("\t %s\n" %neighbor['res'][i])

with open(path+'NearestNeighbors.txt')as f2:
    for i,l in enumerate(f2):
        if i>0:
            data['res'].append(l.split()[1])


# In[62]:


df = pd.DataFrame(data)
df


# In[63]:


# save output file with residue numbers
n_atm=len(data['etot'])
#df.to_csv(path+'H_ion/Species.txt',index=None sep='\t', na_rep="none")
np.savetxt(path+folder+'spec.xyz', df, delimiter='  ',fmt='%s', header='%s \n atom  x  y  z  res  etot[kcal/mol]'%n_atm)


# In[64]:


# save output file without residue numbers
ignore_keys = {'res'}
new = {k:v for k,v in data.items() if k not in ignore_keys}
df2 = pd.DataFrame(new)
np.savetxt(path+folder+'energy.xyz', df2, delimiter='  ',fmt='%s', header='atom  x  y  z  etot[kcal/mol]')


# In[ ]:




