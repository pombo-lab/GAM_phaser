#!/bin/bash
#$ -N bwa_index
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=64G

export PATH="/data/pombo/Sasha/Projects/SNP_calling/software/bwa.kit:$PATH"

bwa index mm10.genome.fa
