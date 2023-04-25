#!/bin/bash
#Comment -  to be submitted by: sbatch slurm_job.txt
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=batch

#SBATCH --mail-type=BEGIN,END

#SBATCH --job-name=batchpractice
#Comment - batch job setup complete

#Comment - load a program require modules

module load R-4.1.2
module load gdal-2.3.1-intel
module load mpfr-2.4.2
module load gcc-6.3.0
module load geos-3.6.2
module load intel-19
module load gmp-4.3.2
module load gcc-9.2.0

#Comment - Set working directory
cd $HOME/Documents/NDP_Simulation_batch


R CMD BATCH --no-restore --no-save rscripts/set_parameter_space.R output/Beta_Parameter_Space_log.txt
