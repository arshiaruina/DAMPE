#!/bin/bash
#SBATCH --job-name="DampeVA"
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --partition=centos7
#SBATCH --export=ALL
#SBATCH --output=job_logs/job_%j.out
#SBATCH --error=job_logs/job_%j.err

# Extract Nlines from the file list to a temporary file
FILETOANALYSE=`pwd`/$(basename -- $FILELIST .txt)_$SLURM_ARRAY_TASK_ID.txt
sed -n $(($SLURM_ARRAY_TASK_ID * $Nlines + 1)),$(($SLURM_ARRAY_TASK_ID * $Nlines + $Nlines))p $FILELIST > $FILETOANALYSE

# user area - change according to your needs!
MYCOMMAND=`pwd`"/eta_corr"
OUTPUTDIR=`pwd`"/Results"
mkdir -p $OUTPUTDIR
OUTPUTFILE="$OUTPUTDIR/job_$(basename -- $FILETOANALYSE .txt).root"
#CHANNELS=`pwd`"/bad_chan.txt"

DAMPECOMMAND="$MYCOMMAND $FILETOANALYSE $OUTPUTFILE"

echo "$DAMPECOMMAND"
$DAMPECOMMAND

rm $FILETOANALYSE
