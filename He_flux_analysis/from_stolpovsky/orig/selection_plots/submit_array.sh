#!/bin/bash
#SBATCH --job-name="DampeSelection"
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --partition=centos7
#SBATCH --export=ALL
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Extract Nlines from the file list to a temporary file
FILETOANALYSE=`pwd`/$(basename -- $FILELIST .txt)_$SLURM_ARRAY_TASK_ID.txt
sed -n $(($SLURM_ARRAY_TASK_ID * $Nlines + 1)),$(($SLURM_ARRAY_TASK_ID * $Nlines + $Nlines))p $FILELIST > $FILETOANALYSE

# user area - change according to your needs!
MYCOMMAND=`pwd`"/exe"
OUTPUTDIR=`pwd`"/Results"
OUTPUTFILE="$OUTPUTDIR/job_$(basename -- $FILETOANALYSE .txt).root"
#CHANNELS=`pwd`"/bad_chan.txt"

DAMPECOMMAND="$MYCOMMAND $FILETOANALYSE $OUTPUTFILE"
# $CHANNELS

$DAMPECOMMAND

rm $FILETOANALYSE
