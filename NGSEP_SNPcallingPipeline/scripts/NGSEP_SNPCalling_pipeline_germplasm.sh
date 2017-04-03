#!/bin/sh

### This PIPELINE is for POST PROCESSING the alignment BAM files and SNPS calling using the NGSEP algorithm
### Prior to run this Pipeline you will need to:
### 1) Demultiplex and quality filter the reads
### 2) Align the reads to reference genome
### 3) Filter and sort the reads
### The bam files should be in the format SampleID_sort.bam
### This is the script for analyzing a germplasm collection.

NGSEP=~/Downloads/NGSEP/NGSToolsApp_2.1.5.jar
REF=~/pvulgaris/assembly/bwa_index/Pvulgaris_218.fa
ALN=aln ## Aligned with BWA MEM. The reads are in the format SampleID_sort.bam


### Call SNPs for each alignment file in the folder
### This part will create a bunch of vcf files in the current directory with the name like NGSEP_SAMPLEID.vcf
### Do not use any separator in your file name
for i in $ALN/*sort.bam
do
	newID=$(echo $i|cut -f 2 -d '/'|cut -f 1 -d '_')
	echo 'Start first Calling for' $newID 
	java -jar $NGSEP FindVariants -h 0.0001 -noRep -noRD -noRP -maxBaseQS 30 -minQuality 40 -maxAlnsPerStartPos 100 -sampleId $newID -ignoreXS $REF $i NGSEP_$newID 
done

## Transfer vcf files in 
mkdir mapping
mv NGSEP*vcf mapping

### Find HQ sites
### First generate SEQNAME FILE

grep '>' $REF|cut -f 2 -d '>' > seqNames.txt 

### Merge the variants 
java -jar $NGSEP MergeVariants seqNames.txt variants_list.vcf mapping/NGSEP*vcf

### Now you need to call again the SNPs in the positions identified in the Merge Variants
### This part will create a bunch of files called NGSEP_HQ_SAMPLEID.vcf with SNPs called in the variants_list.vcf 
for i in $ALN/*sort.bam
do
	newID=$(echo $i|cut -f 2 -d '/'|cut -f 1 -d '_')
	echo 'Start second Calling for' $newID
	java -jar $NGSEP FindVariants -knownVariants variants_list.vcf -h 0.0001 -noRep -noRD -noRP -maxBaseQS 30 -maxAlnsPerStartPos 100 -sampleId $newID -ignoreXS $REF $i NGSEP_HQ_$newID  
done

mkdir genotyping
mv NGSEP_HQ*vcf genotyping

### Final mergeVCF files in a single VCF
java -jar $NGSEP MergeVCF seqNames.txt genotyping/NGSEP_HQ*vcf > finalSNPs_all.vcf


