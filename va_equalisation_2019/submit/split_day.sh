#!/bin/bash

mkdir -p ../resources/20181019_split/

DAY_FILE="../resources/20181019.txt"
LINES=$(( `wc -l < $DAY_FILE`/10 ))
LOCATION="../resources/20181019_split/"

split -l $LINES $DAY_FILE $LOCATION

PARTS=`ls $LOCATION | wc -l`

echo "$DAY_FILE splitted into $PARTS parts and stored in $LOCATION"
