#!/bin/sh
while read p;
do
	a=($p)
	module load biology
	module load samtools/1.8
	module load bedtools/2.27.1
	sbatch -J ${a[0]} -o ${a[0]}.o -e ${a[0]}.e -p akundaje,owners,normal --time=2000 ./log.sh ${a[0]} ${a[1]} ${a[2]}
done < $1
