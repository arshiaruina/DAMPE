#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --ntasks=1
#SBATCH --job-name=eta
#SBATCH --mem=2G
#SBATCH --output=logs/log-%j.out
#SBATCH --error=logs/log-%j.err

# The above lines are options for the SLURM batch system
# Submit this script in interactive mode: bash run.sh
# or submit it to the cluster:  sbatch run.sh

echo $(date) - This is $(hostname), executing task

echo "Sourcing CentOS7"
source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh

#echo "Sourcing SL6"
#source /cvmfs/dampe.cern.ch/rhel6-64/etc/setup.sh

hadd -f /beegfs/users/ruina/out/20181019/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/out/20181019/*.root;

echo $(date) - Merge done.
