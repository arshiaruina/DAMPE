#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --ntasks=1
#SBATCH --job-name=merge
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

hadd -f /beegfs/users/ruina/VAequalisation/out/periodApply/20180901_20180909/corrected/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodApply/20180901_20180909/corrected/unmerged/*.root;
hadd -f /beegfs/users/ruina/VAequalisation/out/periodApply/20181201_20181209/corrected/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodApply/20181201_20181209/corrected/unmerged/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodApply/20181001_20181009/corrected/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodApply/20181001_20181009/corrected/unmerged/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodApply/20180901_20180909/not_corrected/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodApply/20180901_20180909/not_corrected/unmerged/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodApply/20181201_20181209/not_corrected/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodApply/20181201_20181209/not_corrected/unmerged/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodB/data_corrected/20181101_20181109/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodB/data_corrected/20181101_20181109/unmerged/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/20181001_20181009/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/20181001_20181009/*.root;
#hadd -f /beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/20181001_20181009/merged/$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/20181001_20181009/*.root;
#for i in /beegfs/users/ruina/VAequalisation/out/periodA/data_selection_cuts/20181001_20181009/*.root;do rm $i;touch $i;done

#hadd -f /beegfs/users/ruina/VAequalisation/out/20181001_20181009/test/AppCorrFac/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/20181001_20181009/test/AppCorrFac/*.root;
#for i in /beegfs/users/ruina/VAequalisation/out/20181001_20181009/test/AppCorrFac/*.root;do rm $i;touch $i;done

#hadd -f /beegfs/users/ruina/VAequalisation/out/20181001_20181009_1/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/20181001_20181009_1/*.root;
#for i in /beegfs/users/ruina/VAequalisation/out/20181001_20181009_1/*.root;do rm $i;touch $i;done
#hadd -f /beegfs/users/ruina/VAequalisation/out/20181001_20181009/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/VAequalisation/out/20181001_20181009/*.root;
#for i in /beegfs/users/ruina/VAequalisation/out/20181001_20181009/*.root;do rm $i;touch $i;done
#hadd -f /beegfs/users/ruina/out/20181019/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/out/20181019/*.root;
#for i in /beegfs/users/ruina/out/20181019/*.root;do rm $i;touch $i;done
#hadd -f /beegfs/users/ruina/out/201810/merged/merged_$(date +"%d%m%y_%H%M%S").root /beegfs/users/ruina/out/201810/*.root
#for i in /beegfs/users/ruina/out/201810/*.root;do rm $i;touch $i;done
echo $(date) - Merge done.
