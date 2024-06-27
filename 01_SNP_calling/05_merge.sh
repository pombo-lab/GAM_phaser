#!/bin/bash
#$ -N cast.merge
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=64G

export PATH="/home/akukalev/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Projects/SNP_calling/software/bwa.kit:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"

samtools merge cast.merged.rmdup.bam /path/to/rmdup/*.rmdup.bam





