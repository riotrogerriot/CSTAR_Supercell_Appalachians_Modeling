#!/bin/bash

# ===== SLURM OPTIONS =====
#SBATCH --job-name=cm1_rriggin
#SBATCH --partition=Orion
#SBATCH --mem=64gb
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=60:00:00
#SBATCH --output=%x_%j.out
# ===== END SLURM OPTIONS =====

#######################################################################
# Main Program - OpenMPI (64-bit / GNU / INFINIBAND)
#######################################################################

#module load openmpi/4.0.3
#module load hdf5/1.10.5-mpi
module load netcdf/4.7.2-mpi

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "======================================================"
echo ""

cd $SLURM_SUBMIT_DIR
srun --mpi=pmix_v3 $SLURM_SUBMIT_DIR/cm1.exe

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

