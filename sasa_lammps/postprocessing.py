# general
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
# ovito specifics
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *

'''
class SASA_postprossing:
    def __init__(self, path):
        self.path = path
'''
        
### Single Atom analysis ###
def neighbor_analysis(path, SASA_outfile, gro_file):
    # Load positions from the spec file
    pos=np.loadtxt(os.path.join(path,SASA_outfile), skiprows=2, usecols=(1,2,3))
    # create dict to store the nighbour informations
    neighbor ={'ID':[], 'ParticleType': [], 'ResidueType':[], 'AtomName':[]}
    # Load the gro file to find the nearest nighbours of the SASA atoms
    pipeline = import_file(os.path.join(path, gro_file))
    # Delete the solvent and Ions
    with open(gro_file, "r") as rf:
        for i,line in enumerate(rf):
            if 'SOL' in line or 'SOD' in line:
                atmnr = i - 2
                break
    # UPO_del_Solv - Expression selection:
    pipeline.modifiers.append(ExpressionSelectionModifier(
        expression = 'ParticleIdentifier>%s'%str(atmnr)))
    #  UPO_del_Solv - Delete selected
    pipeline.modifiers.append(DeleteSelectedModifier())
    d = pipeline.compute()
    # Initialize neighbor finder object and visit the 1 nearest neighbors of each particle.
    finder = NearestNeighborFinder(1, d)
    # Visit particles closest to some spatial point (x,y,z):
    for i in range(len(pos)):
        for neigh in finder.find_at(pos[i]):
            #print(neigh.index)
            neighbor['ID'].append(d.particles['Particle Identifier'][neigh.index])
            neighbor['ParticleType'].append(d.particles['Particle Type'] \
                    .type_by_id(d.particles['Particle Type'][neigh.index]).name)
            neighbor['ResidueType'].append(d.particles['Residue Type'] \
                    .type_by_id(d.particles['Residue Type'][neigh.index]).name)
            neighbor['AtomName'].append(d.particles['Atom Name'] \
                    .type_by_id(d.particles['Atom Name'][neigh.index]).name)
            #neighbor['dist'].append(neigh.distance)
    neighbor = pd.DataFrame(neighbor)
    return neighbor

def atom_analysis(path, SASA_outfile, neighbor):
    # load SASA results
    spec = pd.read_csv(os.path.join(path, SASA_outfile), sep='\s+', skiprows=1).drop(['atom'], axis=1)
    # insert neighbor properties into df with energy results, delete dublicates and sort by lowest energy
    for key in neighbor:
        spec.insert(0, key, neighbor[key])
        # sort all values by lowest energy and define cutoff for strongest interactions (30 values)
    spec=spec.sort_values('eint/eV', ascending=True).drop_duplicates(subset=['ID']).reset_index(drop=True)
    strongest_interaction = spec.iloc[:30]
    # save strongest interactions to file
    strongest_interaction.to_csv('./atom_analysis_total.txt', sep='\t', index=False)
    # count the amount of most attacked particle types and devide by their total number of all attacked (not in general in the protein!)
    total_count = spec['ParticleType'].value_counts()
    s_count =strongest_interaction['ParticleType'].value_counts()
    result = {'labels': [],'total':[], 'percent': [],}
    for item in s_count.index:
        result['labels'].append(str(item))
        result['total'].append(s_count[item])
        result['percent'].append(float(s_count[item]/total_count[item]))
    pd.DataFrame(result).to_csv(os.path.join(path, 'atom_analysis_percent.txt'), 
                                sep='\t', index=False)
    return result

def atom_analysis_plot(path, neighbor, result):
    # set up color list
    colors={'C': 'grey', 'H': 'whitesmoke', 'N': 'b', 'O': 'r', 'S':'y', 'Fe':'organge', 'Mg': 'm'}
    # set up plot
    fig, ax = plt.subplots(1)
    fig.set_size_inches(4, 5)
    ## plot bars ##
    ax.bar(result['labels'],result['total'],edgecolor='black',
            color=[colors[key] for key in result['labels']],)
    #layout
    ax.set_ylabel('Total interaction',fontsize=18)
    ax.set_xlabel('Particle type',fontsize=18)
    ## increase line thickness box and ticks      
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params(labelsize=16, width=1.5, length=6, bottom=True, left=True, right=True,)
    ## minor ticks
    #ax.minorticks_on()
    ax.tick_params(which='minor',width=1, length=6, right=True,)
    #save
    plt.savefig(os.path.join(path, 'attack_vs_atomtype.png'), dpi = 300, bbox_inches='tight')
    plt.show()

### Residue specific analysis ###
def residuelist(path, gro_file):
    # read .gro file
    gro = pd.read_csv(os.path.join(path, gro_file), sep='\s+', skiprows=2, 
                names=['ResType', 'AtomType', 'AtomNumber', 'x', 'y', 'z', 'vx', 'vy', 'vz',] )
    # find the place to cut the dataset (when the waters or ions start) and ignore SOL and counter Ions
    for i,item in enumerate(gro['ResType']):
        if 'SOL' in item or 'SOD' in item:
            index = i
            break
        else:
            pass
    gro = gro[gro.index < index]
    residuelist = pd.DataFrame([re.split(r'(\d+)', s) for s in gro['ResType']]).drop(columns = 0).drop_duplicates(subset=1)
    np.savetxt(os.path.join(path, 'residuelist.txt'), residuelist,fmt='%s', delimiter='  ', 
                header='ResNum ResType', comments='')
    return residuelist

def residue_analysis(path, SASA_outfile, residuelist):
    df = pd.read_csv(os.path.join(path, SASA_outfile), sep='\s+', skiprows=1)
    res = pd.read_csv(os.path.join(path, residuelist), sep='\s+',)
    # sort by lowest energy and identify and remove dublicates of the same residue
    sort= df.sort_values('eint/eV', ascending=True).reset_index(drop=True)
    # only take the 15 lowest
    emin=sort.drop_duplicates(subset=['res']).iloc[0:15]
    for item in emin['res']:
        loc= res.isin([item]).any(axis=1).idxmax()
        resn = res.iloc[loc]['ResType']
        emin = emin.replace(item, resn)
    emin.to_csv(os.path.join(path, 'minimum_energy_values.txt'), sep='\t', index=False)
    # get all residues and count, count residues in emin
    res_names=res['ResType'].reset_index(drop=True)
    res_names = res_names.value_counts()
    count = emin['res'].value_counts()
    
    # get normalized results, attack amount/amount of residue
    result = {'labels': [], 'total':[], 'percent': [],}
    for item in count.index:
        result['percent'].append(float(count[item]/res_names[item]))
        result['total'].append(count[item])
        result['labels'].append(str(item))
    # export attack amount vs residue to file
    pd.DataFrame(result).to_csv('./residue_interaction_probability.txt', sep='\t', index=False)
    return result 
    
    '''
    drop outliners if there are some
    emin = emin.drop([14])
    emin
    have to include that in a way that makes sense somehow
    '''

def residue_analysis_plot(path, result):
    def residue_analysis_total_plot(path, result):
        result = pd.DataFrame(result).sort_values('labels', ascending=True)
    # create color list with fixed color for each residue
    colors={'ALA': '#d6a090',
        'ARG': '#fe3b1e',
        'ASN': '#a12c32',
        'ASP': '#fa2f7a',
        'CYM': '#fb9fda',
        'CYS': '#e61cf7',
        'GLN': '#992f7c',
        'GLU': '#11963b',
        'GLY': '#051155',
        'HEM': '#4f02ec',
        'HISD': '#2d69cb',
        'HISE': '#00a6ee',
        'HISH': '#6febff',
        'ILE': '#08a29a',
        'LEU': '#2a666a',
        'LYS': '#063619',
        'MET': '#000000',
        'MG': '#4a4957',
        'PHE': '#8e7ba4',
        'PRO': '#b7c0ff',
        'SER': '#ffffff',
        'THR': '#acbe9c',
        'TRP': '#827c70',
        'TYR': '#5a3b1c',
        'VAL': '#ae6507'}
    #plot attack amount vs residue
    fig, ax = plt.subplots(1)
    fig.set_size_inches(10, 5)
    ## plot bars ##
    ax.bar(result['labels'],result['total'],edgecolor='black', 
            color=[colors[key] for key in result['labels']])
    #layout
    ax.set_ylabel('Interaction probability / %',fontsize=18)
    ax.set_xlabel('Residue',fontsize=18)
    ## increase line thickness box and ticks      
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params(labelsize=16, width=1.5, length=6, bottom=True, left=True, right=True,)
    ## minor ticks
    #ax.minorticks_on()
    ax.tick_params(which='minor',width=1, length=6, right=True,)
    #save
    plt.savefig(os.path.join(path, 'total_residue_attack.png'), dpi = 300, bbox_inches='tight')
    #plt.show()
