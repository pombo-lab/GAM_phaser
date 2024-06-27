#!/bin/bash
#$ -N GAM_phasing
#$ -S /bin/bash
#$ -cwd
#$ -l m_mem_free=32G
#$ -l data,os=any

#=====================================================================================================================
# Script written by Alexander Kukalev 
# This script is part of the study: Irastorza-Azcarate I,  Kukalev A,  Kempfer R et al  
# "Extensive folding variability between homologous chromosomes in mammalian cells", bioRxiv, 2024
# The script perform window calling and creates phased GAM segregation tables from the phased and nonphased coverage tables
#=====================================================================================================================

export PATH="$HOME/.guix-profile/bin:$HOME/.guix-profile/sbin${PATH:+:}$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"
export PYTHONPATH="$HOME/.guix-profile/lib/python3.10/site-packages${PYTHONPATH:+:}$PYTHONPATH"

nonphased_coverage=$1
genome1_coverage=$2
genome2_coverage=$3
dataset_name=$4

python3.10 GAM_phasing.Segregation_tables.py $nonphased_coverage $genome1_coverage $genome2_coverage $dataset_name



