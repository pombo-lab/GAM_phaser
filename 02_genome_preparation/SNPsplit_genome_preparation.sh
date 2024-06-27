#!/bin/bash
#$ -N SNP_genome
#$ -S /bin/bash
#$ -cwd
#$ -l m_mem_free=32G
#$ -l data,os=any

perl SNPsplit_genome_preparation --dual_hybrid --strain CAST --strain2 S129_SvJae --reference_genome /path/to/mm10/genome/fasta --vcf_file cast.S129.mm10.filtered.merged.final.vcf
