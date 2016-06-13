### Script for calculating ancestry coefficient in the different genes
### I will run the analysis and then I will reload the project later
library(LEA)
## 1 st convert the vcf to a geno file format
message('Starting vcf file conversion')
vcf2geno(<VCFFile of INTEREST>)
message('File conversion DONE!\n\nStarting sNMF analysis')

## 2nd start the sNMF analysis using the new file created and cluster from 2:15
project=NULL
project=snmf('SNPs_4_analysis.recode.geno', K=2:15, entropy=T, rep=10, project='new', iterations=10000)
