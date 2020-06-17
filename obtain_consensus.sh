#!/bin/bash


set -e
set -o pipefail

if [ $# -ne 1 ]; then
        echo ""
        echo "    Usage: ${0} /path/to/sample_directory"
        echo ""
        echo "    ex.) ${0} /home/my_name/my_sample/"
        echo ""
        exit 1
fi

##
# Set PATH
##


# Programs required:
#  fastp, bwa, ivar, bcftools,samtools ,snpEff,water(EMBOSS) 
#  you can install those softwares with conda.
#  e.g. 
#   $ conda create -n python=3.7 covid19 fastp bwa ivar bcftools samtools emboss
#   $ conda activate covid19
# Before running this script, please make snpEff database for covid19.

#Don't forget to set path for SnpEff directory and ploidy.txt


#base_path=/usr/bin

#fastp_path=${base_path}/fastp
#bwa_path=${base_path}/bwa
#samtools_path=${base_path}/samtools-1.9
#bcftools_path=${base_path}/bcftools-1.9
#$ivar_path=${base_path}/ivar
#PATH=${bwa_path}:${samtools_path}:${bcftools_path}:$PATH
#export ${PATH}

snpEff_path=/opt/covid19/snpEFF


##
#
# Prepare
#
##

sample_directory=${1}
if [ ! -d ${sample_directory} ]; then
	echo "Error: No such directory: ${sample_directory}"
	exit 1
fi

cd ${sample_directory}

sample_name=$(basename ${sample_directory})

 
_sample_fastq_1_file=$(ls -1 ${sample_directory}/*1.fastq.gz | head -1 )
if [ ! -f ${_sample_fastq_1_file} ]; then
	echo "Error: No such file: ${sample_fastq_1_file}"
	exit 1
fi

_sample_fastq_2_file=$(ls -1 ${sample_directory}/*2.fastq.gz | head -1 )
if [ ! -f ${_sample_fastq_2_file} ]; then
	echo "Error: No such file: ${sample_fastq_2_file}"
	exit 1
fi


#if $_sample_fastq_1_file is blank,  exit
if [ -z ${_sample_fastq_1_file} ] && [ -z ${_sample_fastq_2_file} ]; then
  echo "Error: no valid fastq file: ${sample_fastq_1_file}"
  exit 1
fi


#########################################################
#
# Start Analysis
#
########################################################

# Set parameters

#Reference 

# Reference : MN908947.3.fasta
# https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3

ref_data_direcotry=/opt/covid19/reference_data
ref_genome=${ref_data_direcotry}/MN908947.fasta
bedfile=${ref_data_direcotry}/nCoV-2019.bed
gff_genome=${ref_data_direcotry}/covid19.gff3
ploidy_file=${ref_data_directory}/ploidy.txt
cpu=1

#iVar Params

#Length of illumina reads to keep after primer trimming
illuminaKeepLen=20
# Sliding window quality threshold for keeping reads after primer trimming (illumina)
illuminaQualThreshold=20
# Mpileup depth for ivar
mpileupDepth=100000
# iVar frequency threshold for consensus variant (ivar consensus: -t)
ivarFreqThreshold=0.75
# Minimum coverage depth to call variant (ivar consensus: -m; ivar variants -m) 
ivarMinDepth=20
# iVar frequency threshold to call variant (ivar variants: -t )
ivarMinFreqThreshold=0.25
# iVar minimum mapQ to call variant (ivar variants: -q)
ivarMinVariantQuality=20

#
# Step 1: read trimming 
#

#fastp

sample_fastq_1_file=${_sample_fastq_1_file%.fastq.gz}_val_1.fq.gz
sample_fastq_2_file=${_sample_fastq_2_file%.fastq.gz}_val_2.fq.gz

fastp -i ${_sample_fastq_1_file} -I ${_sample_fastq_2_file} \
      -o ${sample_fastq_1_file} -O ${sample_fastq_2_file} \
      -l 20 --low_complexity_filter --thread 2

#
# Step 2:  Alignment with a reference genome file 
#

if [ ! -f ${ref_genome} ]; then
	echo "Error: No such reference genome file: ${ref_genome}"
	exit 1
fi

# create index file for the reference genome, if not.
if [ ! -f ${ref_genome}.bwt ]; then
	bwa index ${ref_genome}
fi


# Alinment by BWA
bwa mem -t ${cpus} ${ref_genome} ${sample_fastq_1_file} ${sample_fastq_2_file} \
  | samtools sort  -O BAM - > ${sample_name}.sorted.bam
samtools index ${sample_name}.sorted.bam

if [ ! -f ${ref_genome}.fai ]; then
	samtools faidx ${ref_genome}
	exit
fi


#
# Step 3: iVar::Trim Primer Sequences
#

samtools view -F4 -o ${sample_name}.mapped.bam ${sample_name}.sorted.bam
samtools index ${sample_name}.mapped.bam
ivar trim -i ${sample_name}.mapped.bam -b ${bedfile} -m ${illuminaKeepLen} -q ${illuminaQualThreshold} -p ivar.out
samtools sort -o ${sample_name}.mapped.primertrimmed.sorted.bam ivar.out.bam
samtools index ${sample_name}.mapped.primertrimmed.sorted.bam

#
# Step 4: iVar:: Call Varianttions
#
samtools mpileup -A -d 0 --reference ${ref_genome} -B -Q 0 ${sample_name}.mapped.primertrimmed.sorted.bam |\
  ivar variants -r ${ref_genome} -m ${ivarMinDepth} -p ${sample_name}_ivar -q ${ivarMinVariantQuality} -t ${ivarMinFreqThreshold}

#
# Step 5: iVar::Make Consensus
#
samtools mpileup -A -d ${mpileupDepth} -Q0 ${sample_name}.mapped.primertrimmed.sorted.bam | \
  ivar consensus -t ${ivarFreqThreshold} -m ${ivarMinDepth} -n N -p ${sample_name}.primertrimmed.ivar.consensus

#
# Step 6:
# bcftools::  Call Variations
#

bcftools mpileup -Ou -A -B -Q 0 -d 0  -f ${ref_genome} ${sample_name}.mapped.primertrimmed.sorted.bam \
 | bcftools call -mv -Ou --ploidy-file ${ploidy_file} \
 | bcftools filter -s LowQual -e '%QUAL<20' -Oz -o ${sample_name}_bcftools_calls.vcf.gz

bcftools index ${sample_name}_bcftools_calls.vcf.gz

# normalize indels
bcftools norm -f ${ref_genome} ${sample_name}_bcftools_calls.vcf.gz -Ob -o ${sample_name}_bcftools_calls.norm.bcf

bcftools index ${sample_name}_bcftools_calls.norm.bcf

# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 ${sample_name}_bcftools_calls.norm.bcf -Ob -o ${sample_name}_bcftools_calls.norm.flt-indels.bcf

bcftools index ${sample_name}_bcftools_calls.norm.flt-indels.bcf

# extract PASS variations
bcftools view -f PASS ${sample_name}_bcftools_calls.norm.flt-indels.bcf > ${sample_name}_bcftools_calls.norm.flt-indels.vcf


# Step 7: 
# Annotation ( by snpEff)
#

java -Xmx4g -jar ${snpEff_path}/snpEff.jar  \
	-no-downstream -no-upstream -no-utr -classic -formatEff covid19 \
	${sample_name}_bcftools_calls.norm.flt-indels.vcf > ${sample_name}_bcftools_calls.norm.flt-indels.eff.vcf

#
# Step 8: pairwise alignment ( by EMBOSS Water)
#

# Reference genome fasta file  v.s. ivar consensus fasta file

water -outfile stdout -gapopen 10 -gapextend 0.5 \
	 ${ref_genome} \
	 ${sample_name}.primertrimmed.ivar.consensus.fa > ${sample_name}_water_ref_vs_ivar.txt

exit 0

