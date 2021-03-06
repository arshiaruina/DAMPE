#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --ntasks=1
#SBATCH --job-name=vaeq12
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

dampe_init

export LD_LIBRARY_PATH=/cvmfs/dampe.cern.ch/centos7/opt/DMPSW/latest/lib:${LD_LIBRARY_PATH}

EXECUTABLE=$1
SPLIT_FILE=$2
while read LINE; do
	echo $LINE
	./$EXECUTABLE $LINE
done < $SPLIT_FILE

echo $(date +"%d%m%y_%H%M%S") - All done.
