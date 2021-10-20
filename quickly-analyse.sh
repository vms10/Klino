#!/bin/bash

currentDir=`pwd`

for expNo in `seq 0 99`

do
	cd $currentDir; 
	cd $expNo; 
	qsub -cwd -b y ../$1
	sleep 10
done

