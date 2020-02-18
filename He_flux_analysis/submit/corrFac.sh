#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --ntasks=1
#SBATCH --job-name=CF10corr
#SBATCH --mem=2G
#SBATCH --output=logs/log-%j.out
#SBATCH --error=logs/log-%j.err

# The above lines are options for the SLURM batch system
# Submit this script in interactive mode: bash run.sh
# or submit it to the cluster:  sbatch run.sh

echo $(date) - This is $(hostname), executing task

EXECUTABLE=../corrFac_09
./$EXECUTABLE 

echo $(date +"%d%m%y_%H%M%S") - All done.
