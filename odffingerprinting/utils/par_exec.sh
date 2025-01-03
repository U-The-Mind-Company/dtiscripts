#! /bin/bash
# Steven Baete, January 2018, NYU SOM CBI

INPUT=$1
PROC=$2

while read -r LINE; do
	while [ $(pgrep -f $PROC | wc -l) -ge 8 ]; do sleep 0.000001 ; done
	#echo $(pgrep -f $PROC | wc -w)
	#echo "$LINE"
	$LINE
done < "$INPUT"
