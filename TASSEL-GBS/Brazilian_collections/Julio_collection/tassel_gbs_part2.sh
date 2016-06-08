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

### After this see snpCall_Impute.sh
### I recalled and Imputed the SNPS with that script

mkdir hapmap
cd hapmap
mkdir raw
mkdir mergedSNP
mkdir filt
cd ..

$TASSEL -Xmx12g -fork1 -DiscoverySNPCallerPlugin -i mergedTBT/GBSWild.tbt.byte -y -m topm/GBSmastertags.topm -mUpd topm/GBSMasterTagsWithVariants.topm -o hapmap/raw/GBSGenos_chr+.hmp.txt -mnF 0.8 -mnLCov 0.5 -mnMAF 0.01 -ref ~/pvulgaris/assembly/tassel/tasselPvul.fa -sC 1 -eC 11 -endPlugin -runfork1 | tee DiscoverySNPs.log

$TASSEL -Xmx12g -fork1 -MergeDuplicateSNPsPlugin -hmp hapmap/raw/GBSGenos_chr+.hmp.txt -o hapmap/mergedSNP/GBSGenos_mergedSNPs_chr+.hmp.txt -sC 1 -eC 11 -endPlugin -runfork1 |tee mergeSNPs.log

$TASSEL -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp hapmap/mergedSNP/GBSGenos_mergedSNPs_chr+.hmp.txt -o hapmap/filt/myGBSGenos_mergedSNPsFilt_chr+.hmp.txt -mnTCov 0.1 -mnSCov 0.1 -mnMAF 0.05 -mnF 0.8 -sC 1 -eC 11 -endPlugin -runfork1 | tee filteredSNPs.log
