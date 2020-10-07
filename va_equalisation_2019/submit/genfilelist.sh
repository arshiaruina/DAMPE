#!/bin/bash

mkdir -p DataList

y=20$1
m=$(printf "%02d" $2)

# Number of days in the month
arr=("01" "03" "05" "07" "08" "10" "12")
if [[ ${arr[*]} =~ $m ]]; then
    ndays=31
elif [[ $m =~ "02" ]]; then
    if [[$1 -eq 16]]; then
        ndays=29
    else
        ndays=28
    fi
else
    ndays=30
fi

# for each day in the month
for i in `seq -w 1 $ndays`
do
    mkdir -p DataList/$y/$m
    f=DataList/$y/$m/data_vaeq_$y$m$i.txt
    find /beegfs/dampe/prod/FM/FlightData/2A/$y$m$i/DAMPE*/DAMPE*.root > $f
    split -l$((`wc -l < $f`/10)) $f ${f}_split -da 9
done
