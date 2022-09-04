#!/bin/bash

# Job name:
#SBATCH --job-name=YourJobname
#
# Project:
#SBATCH --account=nnXXXXk
#
# Wall time limit (DD-HH:MM:SS):
#SBATCH --time=01-00:00:00
#
# Number of processors and memory
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=2
#
# Output file:
#SBATCH --output=/cluster/home/YourUsername/%x-%j.txt

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load SomeProgram/SomeVersion # E.g. FEniCS
module list


## Move results after simulation
cleanup "cp -r $SCRATCH/results_artery $USERWORK/results_artery"
cd $SCRATCH


mesh_path=/cluster/home/YourUserame/Case_test_artery/artery.xml.gz

## Run Oasis
srun oasis NSfracStep problem=Artery mesh_path=$mesh_path

## Run post-processing
srun python automatedPostProcessing/compute_hemodynamic_indices.py --case simulation/results_artery/artery/data/1/Solutions
srun python automatedPostProcessing/compute_flow_and_simulation_metrics.py --case simulation/results_artery/artery/data/1/Solutions




