#!/bin/bash
#$ -N Nmasked_mapping
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -l m_mem_free=12G
#$ -l data,os=any

#=====================================================================================================================
# Script written by Alexander Kukalev 
# This script is part of the study: Irastorza-Azcarate I,  Kukalev A,  Kempfer R et al  
# "Extensive folding variability between homologous chromosomes in mammalian cells", bioRxiv, 2024
# The script map reads to N-masked F123 mouse genome, filter and remove PCR duplicates from raw GAM fastq-files
# The script is also able to produce contamination stats and generate bigwig files
# The script takes 3 input variables:
# a) path to the folder with raw fastq-files (required)
# b) true if user wants to get contamination stats (default true)
# c) true if user wants to generate bigwig files for every sample (default false)
#=====================================================================================================================

FOLDER=$1
fastqscreen=${2:-true}
bigwigs=${3:-false}

export PATH="$HOME/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Common_files:/data/pombo/Sasha/Common_files/fastqscreen:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"
export PYTHONPATH="$HOME/.guix-profile/lib/python3.8/site-packages${PYTHONPATH:+:}$PYTHONPATH"
BOWTIE_INDEX=/path/to/F123_N_masked_bowtie2_indexes/F123_mm10_N-masked
CHROM_SIZES=/path/to/mm10.chrom.sizes

outdir=fastqcreen_reports

if $fastqscreen; then
	mkdir $outdir
fi

FASTQ_FILES="${FOLDER}*.fastq.gz"

for f in $FASTQ_FILES
do
	echo
	echo "Processing" $f
	bowtie2 -p 12 -x $BOWTIE_INDEX $f -S $f.sam && samtools view -q 20 -F 4 -bS $f.sam > $f.bam
	echo "Mapping successfuly done for" $f
	samtools sort $f.bam -o $f.sorted.bam && samtools index $f.sorted.bam && samtools rmdup -s $f.sorted.bam $f.rmdup.bam && samtools index $f.rmdup.bam
	echo "PCR duplicates removed for" $f.rmdup.bam
	if $bigwigs; then
		echo
		bamToBed -i $f.rmdup.bam > $f.bed
		sort -k 1,1 $f.bed -o $f.sorted.bed
		genomeCoverageBed -bg -i $f.sorted.bed -g $CHROM_SIZES > $f.bedgraph
		wigToBigWig $f.bedgraph $CHROM_SIZES $f.bw
		echo "Bigwig file was generated" $f.bw
		rm $f.bed $f.sorted.bed $f.bedgraph
	fi
	if $fastqscreen; then
		fastq_screen --force --quiet --threads 6 --outdir $outdir $f
		echo 'Fastqscreen was done for ' $f
	fi
	rm $f.sam $f.bam
done
