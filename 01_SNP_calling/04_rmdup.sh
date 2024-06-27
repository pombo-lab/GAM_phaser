#!/bin/bash
#$ -N cast.rmdup
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=32G

export PATH="/home/akukalev/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Projects/SNP_calling/software/bwa.kit:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"
FILES="/path/to/sam/*.sam"

for f in $FILES
do
	echo $f
	samtools view -bhS $f > $f.bam && samtools sort $f.bam > $f.sorted.bam && samtools rmdup -S $f.sorted.bam $f.rmdup.bam
done



