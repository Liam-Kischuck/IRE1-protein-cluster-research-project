#!/bin/bash -l
#SBATCH --job-name=matlab_test
#SBATCH --account= # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=2-12:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=4G             # adjust this according to the memory requirement per node you need



module load matlab/2021a.5

matlab -nodisplay -singleCompThread -r "KMC_Transitions_0"
