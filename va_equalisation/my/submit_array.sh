#!/bin/bash
#SBATCH --job-name="DampeVA"
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --partition=centos7
#SBATCH --export=ALL
#SBATCH --output=job_logs/202001_corr/job_%j.out
#SBATCH --error=job_logs/202001_corr/job_%j.err

# Extract Nlines from the file list to a temporary file
FILETOANALYSE=`pwd`/$(basename -- $FILELIST .txt)_$SLURM_ARRAY_TASK_ID.txt
sed -n $(($SLURM_ARRAY_TASK_ID * $Nlines + 1)),$(($SLURM_ARRAY_TASK_ID * $Nlines + $Nlines))p $FILELIST > $FILETOANALYSE

# user area - change according to your needs!
MYCOMMAND=`pwd`"/fill_va_hist"
OUTPUTDIR=`pwd`"/Results/202001_corr"
mkdir -p $OUTPUTDIR
OUTPUTFILE="$OUTPUTDIR/job_$(basename -- $FILETOANALYSE .txt).root"
CORRFILE="/atlas/users/ruina/DAMPE/va_equalisation/ruina/my/Results/corrFac/201601_201602_020920.root"
#CORRFILE="/atlas/users/ruina/DAMPE/va_equalisation/ruina/my/Results/corrFac/withEnergyAdc_201801_201802_250820.root"
#CHANNELS=`pwd`"/bad_chan.txt"

#DAMPECOMMAND="$MYCOMMAND $FILETOANALYSE $OUTPUTFILE 0" # No correction applied. Fix it later
DAMPECOMMAND="$MYCOMMAND $FILETOANALYSE $OUTPUTFILE $CORRFILE"
# $CHANNELS

$DAMPECOMMAND

rm $FILETOANALYSE
