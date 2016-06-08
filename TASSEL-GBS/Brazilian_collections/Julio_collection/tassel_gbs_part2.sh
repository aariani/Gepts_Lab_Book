#!/bin/sh

##### second part of tassel script
##### set TASSEL variable
TASSEL=~/Downloads/tassel4_mod/tassel4-src/run_pipeline.pl

#### convert bam files in TOPM

mkdir topm
$TASSEL -fork1 -SAMConverterPlugin -i aln/GBSTags_uniq.sam -o topm/GBSmastertags.topm -endPlugin -runfork1
echo SAM conversion Done
#### make single tags by taxa file
mkdir tbt
$TASSEL -Xmx12g -fork1 -FastqToTBTPlugin -i fastq -k barcodes.csv -e cviaii -o tbt -y -c 5 -t mergedTagCounts/GBSTags.cnt -endPlugin -runfork1 | tee fastToTBT.log
echo TBT done
### merge tags by taxa files
mkdir mergedTBT
$TASSEL -Xmx12g -fork1 -MergeTagsByTaxaFilesPlugin -i tbt -o mergedTBT/GBSWild.tbt.byte -endPlugin -runfork1 | tee MetgeTBT.log

### SNP calling
### I obtained few SNPs ~6K on 85 individuals. I will save it in oldSNPcall folder and redo the last part by relaxing some parameters
### I will merge the SNP call and imputation step together in the script snpCall_Impute.sh
