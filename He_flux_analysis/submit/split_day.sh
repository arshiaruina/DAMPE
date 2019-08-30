#!/bin/bash

mkdir -p ../resources/20181019_split/

MONTH_FILE="../resources/20181019.txt"
LINES=$(( `wc -l < $MONTH_FILE`/10 ))
LOCATION="../resources/20181019_split/"

split -l $LINES $MONTH_FILE $LOCATION

PARTS=`ls $LOCATION | wc -l`

echo "$MONTH_FILE splitted into $PARTS parts and stored in $LOCATION"
