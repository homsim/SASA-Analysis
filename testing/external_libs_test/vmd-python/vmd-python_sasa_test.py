from vmd import atomsel, molecule
import numpy as np

mol = molecule.load('xyz', 'Cvi_nowater.xyz')
sel = atomsel()

_, sasa_points = sel.sasa(srad = 1.4, 
                          samples = 100,
                          points = True)

export_points = np.array(sasa_points, dtype=str)
export_points = np.insert(export_points, 0, 'He', axis=1)
header = f'{len(export_points)}\n '
np.savetxt('./sasa.txt', export_points, header=header, comments='', fmt='%s')

