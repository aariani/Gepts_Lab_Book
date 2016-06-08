#!/bin/sh

## SET up Variables
TASSEL5=<PATH_TO_TASSEL5>
### Final step of SNPs imputation
cd hapmap/filter
#Concatenate the files:
grep '#' myGBSGenos_mergedSNPsFilt_chr1.hmp.txt > filtSNPs_all.hmp.txt
grep -v '#' * |cut -f2 -d ':'|sort -k3,3n -k4,4n >> filtSNPs_all.hmp.txt

#Impute
#Part1)
~/TASSEL5/run_pipeline.pl -FILLINFindHaplotypesPlugin -hmp filtSNPs_all.hmp.txt -o imputeSNP

#Part2)
~/TASSEL5/run_pipeline.pl -FILLINImputationPlugin -hmp filtSNPs_all.hmp.txt -d imputeSNP -o finalSNPs_bras.hmp.txt



Last Filtering:

Filtering with tassel
~/Downloads/tassel4_mod/tassel4-src/run_pipeline.pl  -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp finalSNPs_bras.hmp.txt -o final_LowMiss_SNPs.hmp.txt -mnSCov 0.8 -mnMAF 0.05 -sC 1 -eC 11 -endPlugin -runfork1


Convert to VCF
~/TASSEL5/run_pipeline.pl -fork1 -h final_LowMiss_SNPs.hmp.txt -export -exportType VCF -runfork1
