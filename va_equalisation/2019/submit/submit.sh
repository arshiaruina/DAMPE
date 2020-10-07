#!/bin/bash
filename='../resources/20181019.txt'
while read line; do
	echo $line
	sbatch run.sh $line  
#	./eta.exe $line
done < $filename
#for i in /beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_*/*.root; do
#	echo $i
#	sbatch run.sh $i  
#done
