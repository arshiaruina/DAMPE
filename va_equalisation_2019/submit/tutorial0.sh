#!/bin/bash
# a shbang

ourfile=$0
echo $ourfile

NAME=$1
GREETING="Hi there"
HAT_TIP="*tip of the hat*"
HEAD_SHAKE="*slow head shake*"

if [ "$NAME" = "Dave" ]; then #single square brackets works for me unlike Misha's double brackets
	echo $GREETING
elif [ "$NAME" = "Steve" ]; then
	echo $HAT_TIP
else
	echo $HEAD_SHAKE
fi
