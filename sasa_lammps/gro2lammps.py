from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import os
import pandas as pd

'''
To succesfully run gro2lammps you need a library file of all elements included in you system.
This files needs the ElementType (H,C,O, ...), the gromacs_ParticleType (1, 6, 8, ...), 
the lammps_ParticleType (1,2,3, ...) and the mass (1.00, 12.0, 16.0 ...). 
You should doubble check first which gromacs_ParticleType number ovito gives the elements.
It is mostly related to the atomic number, but for Mg, for example, gromacs_ParticleType is 0.
So check that first. 
You can choose any lammps_ParticleType you want for each of the elements,
but the numbers have to be consecutive (1,2,3, ...).
And for a follow up ReaxFF simulation with lammps the elements have to be listed in the right order
in the "pair_coeff" command. Resulting in: "pair_coeff     * * protein2013.ff H C O" 
if you use the protein2013.ff force field and H has the lammps_ParticleType 1, C is 2 and O is 3. 
'''

class gro2lammps:
    def __init__(self, path, element_library):
        #_init: load data from libary file and write it to dictionary
        self.path = path
        self.di = pd.read_csv(os.path.join(path, element_library), sep='\s+')
    
    def __delete_solvent(self, infile, pipeline):
        # Find the place to cut
        with open(infile, "r") as rf:
            for i,line in enumerate(rf):
                if 'SOL' in line or 'SOD' in line:
                    atmnr = i - 2
                    break
        # UPO_del_Solv - Expression selection:
        pipeline.modifiers.append(ExpressionSelectionModifier(
            expression = 'ParticleIdentifier>%s'%str(atmnr)))
        #  UPO_del_Solv - Delete selected
        pipeline.modifiers.append(DeleteSelectedModifier())

    def __change_ParticleTypes(self, frame, data):
        types = data.particles_.particle_types_
        for gro_PT, lammps_PT in zip(self.di['gromacs_ParticleType'], self.di['lammps_ParticleType']):
            types[types == gro_PT ] = lammps_PT
    
    def __change_ParticleIDs(self, frame, data):
        for i,item in enumerate(self.di['gromacs_ParticleType']):
            data.particles_.particle_types_.type_by_id_(item).id = i+1
    
    def __change_Masses(self, frame, data):
        for lammps_PT, et, mass in zip(self.di['lammps_ParticleType'], self.di['ElementType'], self.di['mass']):
            data.particles_.particle_types_.type_by_id_(lammps_PT).name = et
            data.particles_.particle_types_.type_by_id_(lammps_PT).mass = mass
    
    def convert(self,infile, outfile):
        # Data import:
        pipeline = import_file(os.path.join(self.path,infile))

        # Delete solvent and ions
        self.__delete_solvent(infile, pipeline)

        # Change Particle IDs
        pipeline.modifiers.append(self.__change_ParticleTypes)

        # Change Particle IDs
        pipeline.modifiers.append(self.__change_ParticleIDs)

        #Change Masses and Names
        pipeline.modifiers.append(self.__change_Masses)

        # Store residue informations in the Molecule Identifier 
        pipeline.modifiers.append(ComputePropertyModifier(
            expressions = ('ResidueIdentifier',), 
            output_property = 'Molecule Identifier'))
        data = pipeline.compute()

        export_file(pipeline, os.path.join(self.path, outfile), "lammps/data",  atom_style="full")
