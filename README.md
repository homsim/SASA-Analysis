# SASA-Analysis
Analyses interaction energyies between protein/macromolecules and probe atom/molecule at the solvent accessible surface area of the molecule
Interaction energies are determined by using ReaxFF in the LAMMPS simulation suite.
## Step 1
- prepare .xyz of the protein/macro molecule
- file should contain only the molecule in vaccum (no water or other solvent molecules)

## Step 2
- Generate SASA coordinates using the vmd.tcl file
- run file with: ``` vmd -dispdev text -eofexit < vmd.tcl ```
    - therefore .xyz strcutre of the protein/enzyme has to be in the same folder and names in the .tcl file

## Step 3
- one atomic probe molecule:
    - run lammps_SASA-1atom.py in respective folder (h-ion in the example)
    - right .mol file, .data file and forcefield file has to be in that folder too
- more atomic probe molecules:
    - run SASA_rotationfile.py first to get right orientation
    - run lammps_SASA-run.py in respective folder (h202 in the example)
    - right .mol file, .data file, rotation.xyz and forcefield file has to be in that folder too

  ## Step 4
  - Interaction energies of the single point calculations are stored in etot file
  - Use SASA_OutputAnalysis.py to generate an xyz file with:
    -  the coordinates of the probe molecules
    -  the interaction energies
    -  the respective molecule residues which are closest do the probe molecule
  -  For running SASA_OutputAnalysis.py you need this files:
    -  etot
    -  sasa.xyz
    -  NearestNeighbors.txt
    -  data. file
