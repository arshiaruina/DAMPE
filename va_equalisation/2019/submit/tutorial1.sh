#!/bin/bash
#testing and control flow with if, [ and [[, and/or

NUM_REQD_ARGS=2

#Do we have at least two arguments?
if [ $# -lt $NUM_REQD_ARGS ]; then
	echo "Not enough arguments. Call this script with
	./{$0} <name> <number>"
fi

## helpers
# && and
# || or
echo "hi" || echo "this won't happen"
$(ls nonexistentfile) || echo "this happens if the first thing fails"
echo $(pwd) && echo "this also happens"

# Not null (-n) or zero length (-z)
notnully="this is not null"
nully=""

if [ -n "$notnully" ]; then
	echo "This is not null"
fi
if [ -z "$nully" ]; then
	echo "This is null"
fi
