#!/bin/bash

# Job name:
#SBATCH --job-name=YOUR_JOB_NAME
#
#### EDIT NEXT LINE TO YOUR PROJECT ACCOUNT ####
#SBATCH --account=nnXXXXk
#
# Wall time limit (DD-HH:MM:SS):
#SBATCH --time=01-00:00:00
#
# Number of processors and memory
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=2
#
#### EDIT NEXT LINES TO YOUR USERNAME ####
# Set output log location
#SBATCH --output=/cluster/home/USERNAME/Oasis/%x-%j.txt

# Set username on cluster
USERNAME=USERNAME

## Set up job environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

# Move simulation files
mv probe Artery.py ICA_values Probe.py Womersley.py $SCRATCH

# Move post-processing files
mv compute_flow_and_simulation_metrics.py compute_hemodynamic_indices.py postprocessing_common.py visualize_probes.py $SCRATCH

# Reset the modules to the system default
module --quiet purge

#### EDIT NEXT LINES TO LOAD FENICS ####
source PATH_TO_FENICS_CONFIGURATION.conf

# Install Oasis
export OASIS_PREFIX=/cluster/home/$USERNAME/Oasis
export PATH=$OASIS_PREFIX/bin:$PATH
export PYTHONPATH=$OASIS_PREFIX/lib/python3.6/site-packages:$PYTHONPATH
cd $OASIS_PREFIX
pip3 install --prefix=$OASIS_PREFIX .

module list

## Move results after simulation
cleanup "cp -r $SCRATCH/results_artery $USERWORK/"
cd $SCRATCH

# Define mesh path
mesh_path=/cluster/home/$USERNAME/Oasis/oasis/mesh/artery.xml.gz

## Run Oasis
srun oasis NSfracStep problem=Artery mesh_path=$mesh_path

## Run post-processing
python compute_hemodynamic_indices.py --case $SCRATCH/results_artery/artery/data/1/Solutions
python compute_flow_and_simulation_metrics.py --case $SCRATCH/results_artery/artery/data/1/Solutions
python visualize_probes.py --case $SCRATCH/results_artery/artery/data/1/Probes
