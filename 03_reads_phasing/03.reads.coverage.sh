#!/bin/bash
#$ -N mm10.coverage
#$ -S /bin/bash
#$ -cwd
#$ -l m_mem_free=32G
#$ -l data,os=any

#=====================================================================================================================
# Script written by Alexander Kukalev 
# This script is part of the study: Irastorza-Azcarate I,  Kukalev A,  Kempfer R et al  
# "Extensive folding variability between homologous chromosomes in mammalian cells", bioRxiv, 2024
# The script compute reads and coverage GAM tables for window phasing at multiple resolutions
#=====================================================================================================================

export PATH="$HOME/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Common_files:/data/pombo/Sasha/Common_files/fastqscreen:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"
export PYTHONPATH="$HOME/.guix-profile/lib/python3.8/site-packages${PYTHONPATH:+:}$PYTHONPATH"
CHROM_SIZES=/path/to/mm10.chrom.sizes
WINDOWS="30000 40000 50000 100000 200000 250000 500000 "

FOLDER=$1
NAME=$2

#=========================================================================================

BAM_FILES="${FOLDER}*.bam"

for WINDOWSIZE in $WINDOWS:
do
	bedtools makewindows -g $CHROM_SIZES -w $WINDOWSIZE > ${NAME}.at.${WINDOWSIZE}.windows.bed
	cp ${NAME}.at.${WINDOWSIZE}.windows.bed ${NAME}.at.${WINDOWSIZE}.coverage.windows.bed
	
	# Generate coverage tables
	for f in $BAM_FILES
	do
		bedtools coverage -a ${NAME}.at.${WINDOWSIZE}.windows.bed -b $f | awk '{print $5}' > ${NAME}.at.${WINDOWSIZE}.coverage.bed
		paste ${NAME}.at.${WINDOWSIZE}.coverage.windows.bed ${NAME}.at.${WINDOWSIZE}.coverage.bed > ${NAME}.at.${WINDOWSIZE}.merged.bed && mv ${NAME}.at.${WINDOWSIZE}.merged.bed ${NAME}.at.${WINDOWSIZE}.coverage.windows.bed
		echo "Window coverage was calculated for" $f
	done

	echo -e "chrom\tstart\tstop" $BAM_FILES > ${NAME}.at.${WINDOWSIZE}.list_of_files.txt
	tr ' ' \\t < ${NAME}.at.${WINDOWSIZE}.list_of_files.txt > ${NAME}.at.${WINDOWSIZE}.list_of_files.tab.txt && mv ${NAME}.at.${WINDOWSIZE}.list_of_files.tab.txt ${NAME}.at.${WINDOWSIZE}.list_of_files.txt
	cat ${NAME}.at.${WINDOWSIZE}.list_of_files.txt ${NAME}.at.${WINDOWSIZE}.coverage.windows.bed > ${NAME}.at.${WINDOWSIZE}.windows.header.bed && mv ${NAME}.at.${WINDOWSIZE}.windows.header.bed ${NAME}.coverage.at.${WINDOWSIZE}.table
	rm ${NAME}.at.${WINDOWSIZE}.list_of_files.txt ${NAME}.at.${WINDOWSIZE}.coverage.bed ${NAME}.at.${WINDOWSIZE}.coverage.windows.bed ${NAME}.at.${WINDOWSIZE}.windows.bed
done





