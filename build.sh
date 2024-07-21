#####
# This scipt installs the dependencies for the sasa_lammps package.
# I do like this because I have not found a way to install dependencies in the build recipe
# with a strict-channel-priority on the ovito-channel
#####

# install the dependencies
conda install --strict-channel-priority -c https://conda.ovito.org -c conda-forge ovito matplotlib numpy tqdm scipy pandas vmd-python -y
# build the package locally
conda build build_recipe/ -y
# install the package into the current env
conda install --use-local sasa_lammps -y