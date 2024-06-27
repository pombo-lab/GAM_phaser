#!/bin/bash
#$ -N quality_trimming
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=80G

export PATH="/home/akukalev/.guix-profile/bin${PATH:+:}$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"

read1="/path/to/read1/*_1.fastq.gz"
read2="/path/to/read2/*_2.fastq.gz"

set -- $read2
for i in $read1
do
	echo $i $1
	cutadapt -q 20,15 --pair-filter=any -o $i.trimmed.fastq.gz -p $1.trimmed.fastq.gz $i $1
    shift
done

