#!/bin/sh

### Script for final SNP imputation using the TASSEL GBS algorithm
##### set variables
TASSEL4=~/Downloads/tassel4_mod/tassel4-src/run_pipeline.pl
TASSEL5=~/TASSEL5/run_pipeline.pl

### Repeat the SNP calling step for relaxing some parameter
mkdir hapmap
cd hapmap
mkdir raw
mkdir mergedSNP
mkdir filt
cd ..

### Modified nmLCov to 0.1 (default) instead of 0.5, and removed mnMAF param (MAF default is 0.0)
$TASSEL4 -Xmx12g -fork1 -DiscoverySNPCallerPlugin -i mergedTBT/GBSWild.tbt.byte -y -m topm/GBSmastertags.topm -mUpd topm/GBSMasterTagsWithVariants.topm -o hapmap/raw/GBSGenos_chr+.hmp.txt -mnF 0.8 -mnLCov 0.1 -ref ~/pvulgaris/assembly/tassel/tasselPvul.fa -sC 1 -eC 11 -endPlugin -runfork1 | tee DiscoverySNPs.log

$TASSEL4 -Xmx12g -fork1 -MergeDuplicateSNPsPlugin -hmp hapmap/raw/GBSGenos_chr+.hmp.txt -o hapmap/mergedSNP/GBSGenos_mergedSNPs_chr+.hmp.txt -sC 1 -eC 11 -endPlugin -runfork1 |tee mergeSNPs.log

###  Removed mnMAF param (MAF default is 0.0)
$TASSEL4 -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp hapmap/mergedSNP/GBSGenos_mergedSNPs_chr+.hmp.txt -o hapmap/filt/myGBSGenos_mergedSNPsFilt_chr+.hmp.txt -mnTCov 0.1 -mnSCov 0.1 -mnF 0.8 -sC 1 -eC 11 -endPlugin -runfork1 | tee filteredSNPs.log


### IMPUTATION step
### Need to be in the hapmap/filt folder
cd hapmap/filt

### Concatenate and sort files in the different chromosomes
grep '#' myGBSGenos_mergedSNPsFilt_chr1.hmp.txt > filtSNPs_all.hmp.txt
grep -v '#' * |cut -f2 -d ':'|sort -k3,3n -k4,4n >> filtSNPs_all.hmp.txt


### Imputation part1
$TASSEL5 -FILLINFindHaplotypesPlugin -hmp filtSNPs_all.hmp.txt -o imputeSNP

### Imputation part2
$TASSEL5 -FILLINImputationPlugin -hmp filtSNPs_all.hmp.txt -d imputeSNP -o finalSNPs.hmp.txt

### Final SNP filtering with TASSEL4
### Filter my MAF and Missingness only after the IMPUTATION
$TASSEL4  -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp finalSNPs.hmp.txt -o final_LowMiss_SNPs.hmp.txt -mnSCov 0.8 -mnMAF 0.05 -sC 1 -eC 11 -endPlugin -runfork1

## Convert to vcf
$TASSEL5 -fork1 -h final_LowMiss_SNPs.hmp.txt -export -exportType VCF -runfork1
