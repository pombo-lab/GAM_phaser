#!/bin/bash
#$ -N SNP_split
#$ -S /bin/bash
#$ -l m_mem_free=64G
#$ -l data,os=any

FILES=/path/to/mapped/rmdupbam/files/*.rmdup.bam

for f in $FILES
do
	perl SNPsplit --snp_file all_S129_SvJae_SNPs_CAST_reference.based_on_mm10.txt $f
done
