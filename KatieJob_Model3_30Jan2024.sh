#!/bin/bash					
#SBATCH --partition=remi			
#SBATCH --job-name=KatieJob			
#SBATCH --output=KatieJob_%j.out			
#SBATCH --error=KatieJob_%j.err			
#SBATCH --mail-type=ALL				
#SBATCH --mail-user=katie.tseng@wsu.edu		
#SBATCH --nodes=1				
#SBATCH --ntasks=1				
#SBATCH --cpus-per-task=5			
#SBATCH --mem-per-cpu=4G			

module load r/4.1.0				
Rscript --vanilla Kamiak_BRT-Model3_30Jan2024.R 	
echo "Completed job $SLURM_JOBID on node(s) $SLURM_JOB_NODELIST" 

























