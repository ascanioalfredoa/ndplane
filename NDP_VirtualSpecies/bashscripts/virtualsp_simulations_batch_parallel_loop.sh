#!/bin/bash
#Comment -  to be submitted by: sbatch slurm_job.txt
#SBATCH --time=80:00:00
#SBATCH --partition=batch

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --output=slurm_out/array_%A_%a.out
#SBATCH --error=slurm_out/array_%A_%a.err
#SBATCH --array=1-8
#SBATCH --exclusive

#SBATCH --mail-type=BEGIN,END

#SBATCH --job-name=NDP_vsp_sim
#Comment - batch job setup complete

#Comment - load a program require modules

module load R-4.1.2
module load gdal-2.3.1-intel
module load mpfr-2.4.2
module unload gcc-4.9.4
module unload gcc-9.2.0
module load gcc-6.3.0
module load gcc-9.2.0
module load gcc-6.3.0
module load geos-3.6.2
module load intel-19
module load gmp-4.3.2
module load gcc-9.2.0

#Comment - Set working directory
cd $HOME/Documents/NDP_virtualspecies_sim

#Comment - Loop

#for i in {1..8}
#do
    echo "Node $SLURM_ARRAY_TASK_ID is running"
    #echo export INDEX=$SLURM_ARRAY_TASK_ID
    R CMD BATCH --no-restore --no-save rscripts/virtualsp_simulations_batch_parallel_loop.R slurm_out/output$SLURM_ARRAY_TASK_ID
    #R CMD BATCH --no-restore --no-save rscripts/test.R slurm_out/output$SLURM_ARRAY_TASK_ID
    #srun --nodes=1 --ntasks=1 R CMD BATCH --no-restore --no-save rscripts/test.R slurm_out/output$SLURM_ARRAY_TASK_ID
#done

#Comment Deletes temporary folders created by me
ls -l | grep "Rtmp" | grep "ascaniaa" | cut -c 61-70 | xargs rm -rf
