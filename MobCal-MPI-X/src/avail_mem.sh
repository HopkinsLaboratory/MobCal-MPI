#!/bin/bash

# Determine the total amount of system memory

FREECMDOUTPUT=(`free --mega | sed "s/ /_/g"`)
PROCESSORS=`cat /proc/cpuinfo | grep "processor" | wc -l`
MEM=0

for LINE in ${FREECMDOUTPUT[@]}; do
    if [[ $LINE =~ Mem:.* ]]; then
	MEM=(`echo $LINE | sed "s/_/ /g"`)
	python3 -c "print(int(${MEM[3]}/(2 * $PROCESSORS)))"
    fi
done

