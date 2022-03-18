#!/bin/bash 

bam=$1
ref=${2:-../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta}

vcf=${bam}.vcf
normvcf=${bam}.normvcf

# compute read pileup
eval "$(conda shell.bash hook)"
conda activate base 
bcftools mpileup -a INFO/AD --max-idepth 1000000 --max-depth 1000000 --fasta-ref ${ref} ${bam} \
| bcftools call -o ${vcf} -Ov --ploidy 1 --keep-alts --variants-only --multiallelic-caller 
bcftools norm -a --atom-overlaps . ${vcf} > ${normvcf}

cov=${bam}.cov
bedtools genomecov -d -ibam ${bam} > $cov

conda deactivate
