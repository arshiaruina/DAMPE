#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --ntasks=1
#SBATCH --job-name=demo
#SBATCH --mem=2G
#SBATCH --output=log-%j.out

# The above lines are options for the SLURM batch system
# Submit this script in interactive mode: bash run.sh
# or submit it to the cluster:  sbatch run.sh

echo $(date) - This is $(hostname), executing task

echo "Sourcing CentOS7"
source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh

#echo "Sourcing SL6"
#source /cvmfs/dampe.cern.ch/rhel6-64/etc/setup.sh

dampe_init

python demo.py

echo $(date) - All done.
