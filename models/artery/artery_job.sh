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
#SBATCH --output=/cluster/home/USERNAME/%x-%j.txt

## Set up job environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

# Reset the modules to the system default
module --quiet purge

# load the Anaconda3
module load Anaconda3/2019.03

# Set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

# Source the conda environment setup, Miniconda or Anaconda
# source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

# Deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

# Activate the VaMPY environment with Oasis, FEniCS and morphMan installed
conda activate PATH_TO_VAMPY_ENVIRONMENT

# Set username on cluster
USERNAME=USERNAME

# pip install VaMPy
export VAMPY_PREFIX=/cluster/home/$USERNAME/VaMPy
export PATH=$VAMPY_PREFIX/bin:$PATH
export PYTHONPATH=$VAMPY_PREFIX/lib/python3.6/site-packages:$PYTHONPATH
cd $VAMPY_PREFIX
pip3 install --prefix=$VAMPY_PREFIX .

# Move VaMPy directory to SCRATCH
cp -r /cluster/home/$USERNAME/VaMPy $SCRATCH

## Move results after simulation
cleanup "cp -r $SCRATCH/VaMPy/src/vampy/simulation/results_artery $USERWORK/results_artery"

# Move to simulations folder for CFD simulation
cd $SCRATCH/VaMPy/src/vampy/simulation

# Define mesh path
mesh_path=/cluster/home/$USERNAME/VaMPy/src/vampy/simulation/artery.xml.gz

## Run Oasis
mpirun -np $SLURM_NTASKS oasis NSfracStep problem=Artery mesh_path=$mesh_path

## Run post-processing
vampy-convert  --case $SCRATCH/results_artery/artery/data/1/Solutions
vampy-hemo --case $SCRATCH/results_artery/artery/data/1/Solutions
vampy-metrics --case $SCRATCH/results_artery/artery/data/1/Solutions
vampy-probes --case $SCRATCH/results_artery/artery/data/1/Probes
