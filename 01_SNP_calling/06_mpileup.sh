#!/bin/bash
#$ -N cast.mpileup
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=64G

export PATH="/home/akukalev/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Projects/SNP_calling/software/bwa.kit:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"

bcftools mpileup -Ou -f /path/to/fasta/genome/genome.fasta /path/to/rmdup/merged.rmdup.sorted.bam | bcftools call -mv -Ob -o raw.calls.bcf
 



