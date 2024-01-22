#!/bin/bash
#SBATCH --constraint hpc-simulation
#SBATCH --job-name=MRG0080
#SBATCH --nodes=1           
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4 
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1  
#SBATCH --time=3-00:00:00     # Time limit days-hrs:min:sec
#SBATCH --output=%x_%j.log  
#SBATCH --exclusive

echo "Job start"

# Run FEniCS / OasisMove
echo "-- Task setup --"
echo "nTasks"
echo $SLURM_NTASKS
echo "Tasks/node"
echo $SLURM_TASKS_PER_NODE
echo "NTasks/node"
echo $SLURM_NTASKS_PER_NODE
echo "----------------"

# Run FEniCS / OasisMove
date

CASE=0080
CONDITION=AF
echo "WORKING ON CASE ${CASE} CONDITION ${CONDITION} "

echo "COMBINING ENERGY"
/usr/mpi/gcc/openmpi-4.1.2a1/bin/mpirun -np $SLURM_NTASKS --oversubscribe pvbatch /app/VaMPy/scripts/combine_and_convert_energy_to_vtu.py --case $CASE --condition $CONDITION

echo "MERGING ENERGY"
/usr/mpi/gcc/openmpi-4.1.2a1/bin/mpirun -np $SLURM_NTASKS --oversubscribe pvbatch /app/VaMPy/scripts/merge_energy_vtu.py --case $CASE --condition $CONDITION

echo "CONVERTING HEMO"
/usr/mpi/gcc/openmpi-4.1.2a1/bin/mpirun -np $SLURM_NTASKS --oversubscribe pvbatch /app/VaMPy/scripts/combine_and_convert_hemo_to_vtp.py --case $CASE --condition $CONDITION

echo "CONVERTING BRT"
/usr/mpi/gcc/openmpi-4.1.2a1/bin/mpirun -np $SLURM_NTASKS --oversubscribe pvbatch /app/VaMPy/scripts/convert_brt_to_vtu.py --case $CASE --condition $CONDITION


echo "Job done"
date


