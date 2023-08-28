#!/bin/bash
# Job name:
#SBATCH --job-name=Debugging_PP_N_C16_160K_cycles_5to10_4cores_2G_2ndtrial

#
# Max running time (DD-HH:MM:SS)
#SBATCH --time=00-01:30:00 
#
# Project:
#SBATCH --account=nn9249k
#SBATCH --output=/cluster/home/gadursn/VaMPy/Output_log_files/%x-%j.txt
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
source /cluster/home/gadursn/oasis_fenics-2019.1.0.conf
module list

## Move to working directory
WORKDIR=$HOME/VaMPy/automatedPostProcessing
cd $WORKDIR

#Flush output
export PYTHONUNBUFFERED=TRUE

## Run FEniCS
#srun python3 compute_hemodynamic_indices_phaseaveraged.py --case /cluster/work/users/gadursn/N_C16_10Melements/C_16/data/1/PostProc --nu 0.0033 --rheology_model Newtonian --rho 1060 --dt 0.1 --velocity-degree 1 --probe-frequency 100 --T 1000 --save-frequency 5 --start-cycle 2 --sample-step 2

srun python3 compute_hemodynamic_indices_phaseaveraged.py --case /cluster/work/users/gadursn/NN_C16_160Kelements_Carreau_Model/C_16/data/1/PostProc --nu 0.0033 --rheology_model Newtonian --rho 1060 --dt 0.1 --velocity-degree 1 --probe-frequency 100 --T 1000 --save-frequency 5 --start-cycle 2 --sample-step 2
