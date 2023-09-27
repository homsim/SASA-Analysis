#!/usr/bin/env python
# coding: utf-8

# In[17]:


import os
import numpy as np

from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *


# In[21]:


path=r'/home/hanna/simulations/MM/lammps/GapA/04_plasma/SASA/'
folder='OH_ion/'
file='sasa.xyz'
mol_file='oh_minus.mol'


# In[22]:


# calc the vactor standing on mol vector and FeS as roation vector
pos = {'coords':[]}
with open(path+file) as f:
    for i,l in enumerate(f):
        if i>1:
            floats = [float(x) for x in l.split()[1:]]
            coords=np.asarray(floats)
            pos['coords'].append(coords)

neighbor ={'id':[],'pos':[],'res':[],'dist':[]}
# Load input simulation file.
pipeline = import_file(path+'data.GapA_nowater',atom_style = 'full')
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


# In[23]:


# load positions from sasa file    
pos=np.loadtxt(path+file,skiprows=2, usecols=(1,2,3))
# 2 atomic plasma species mols have to be rotated to be orientet right to the surface
# calc the angle between the mol vector and the FeS
#get coordinates from mol file
mol=[]
COUNT_FLAG = False
with open(path+folder+mol_file) as f:
    for i, l in enumerate(f):
        if COUNT_FLAG:
            if "Types" in l:
                COUNT_FLAG = False
                continue
            # ...
            if l=="\n": continue
            floats = [float(x) for x in l.split()]
            mol.append(np.asarray(floats[1:]))
            #...
        if "Coords" in l:
            COUNT_FLAG = True

#calc vectors and rotation angle
rot = {'rotation_angle':[], 'rotation_vector':[]}
if len(mol)>=2:
    for i in range(len(pos)):
        pos_O=(pos[i]+mol[0])
        pos_H=(pos[i]+mol[1])
        OH=pos_H-pos_O
        NN=neighbor['pos'][i]
        N=np.cross(OH,NN)
        rot['rotation_vector'].append('%s %s %s'%(round(N[0],3),round(N[1],3),round(N[2],3)))
        S=np.dot(OH,NN)
        lx=np.sqrt((OH*OH).sum())
        ly=np.sqrt((NN*NN).sum())
        angle=(np.arccos(S/(lx*ly)))* 360 / 2 / np.pi
        if angle <= 90:
            rot['rotation_angle'].append(round(90-angle,3))
        if angle > 90:
            rot['rotation_angle'].append(-round(angle-90,3))

with open(path+folder+'rotation.txt', 'w') as fp:
    fp.write('rotation_angle rotation_vector \n')
    for i in range(len(pos)):
        fp.write('%s %s %s ' %(pos[i][0],pos[i][1],pos[i][2]))
        fp.write('%s ' %rot['rotation_angle'][i])
        fp.write("%s\n" %rot['rotation_vector'][i])


# In[ ]:




