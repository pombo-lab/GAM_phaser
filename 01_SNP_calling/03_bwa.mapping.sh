#!/bin/bash
#$ -N bwa.mem
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=16G

export PATH="/home/akukalev/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Projects/SNP_calling/software/bwa.kit:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"

BWA_mm10_index=/path/to/bwa_indexes/mm10.genome.fa
read1="/path/to/read1/*_1.fastq.gz.trimmed.fastq.gz"
read2="/path/to/read2/*_2.fastq.gz.trimmed.fastq.gz"

set -- $read2
for i in $read1
do
	echo $i $1
	bwa mem -M -t 8 $BWA_mm10_index $i $1 > $i.sam
    shift
done

