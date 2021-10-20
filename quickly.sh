#!/bin/bash

currentDir=`pwd`

for expNo in `seq 0 19`

do
	cd $currentDir; 
	mkdir $expNo; 
	cd $expNo;
	qsub -cwd -b y ../$1
	sleep 30
done

