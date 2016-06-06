#! /usr/bin/env python
#########################################
#########################################
##### Script for create fasta file with the snps position in a multisample vcf file
##### need vcftools
##### Select only SNPs and indels located in CDS sequences (subjected to lower evolutionary pressure)
##### Since bean is an high selfing species we treat called heterozygotes as missing data, because the
##### probably are errors. Keep only genotypes with a quality >=10
##### only for bi-allelic snps
##### 
from __future__ import print_function
import sys, os
from commands import getoutput

vcf=sys.argv[1]#### INPUT VCF
fasta=sys.argv[2]#### OUTPUT fasta
gff=sys.argv[3]


cds=getoutput('grep CDS %s' % gff).split('\n')
print('#cds_coordinates', file=open('cds.bed', 'a'))
for i in cds:
	i=i.split('\t')
	print(i[0], int(i[3]) -1, int(i[4]), sep='\t', file=open('cds.bed', 'a')) #### modify the position tarlucco!!

os.system('vcftools --vcf %s --bed cds.bed --recode --out cds_snp' % vcf)


samplename=getoutput('grep CHROM %s' % vcf).strip().split('\t')[9:]   
#print(samplename)
fastaseq={k:'' for k in samplename}
for line in open('cds_snp.recode.vcf'):
	if '#' not in line:
		line=line.strip().split('\t')
		ref=line[3]
		alt=line[4]
		n=max(len(ref), len(alt))
		if len(ref)!= len(alt):
			if len(ref)==n:
				add_alt='-'*(len(ref)-len(alt))
				add_ref=''
			elif len(alt)==n:
				add_ref='-'*(len(alt)-len(ref))
				add_alt=''
		else:
			add_alt=''
			add_ref=''
		geno=line[9:]
		geno=[i.split(':') for i in geno]
		for i in range(len(samplename)):
			sample=samplename[i]
			snp=list(set(geno[i][0].split('/')))
			if int(geno[i][3])< 10: ##### low quality genotype, like --minGQ 10 in vcftools
				fastaseq[sample]+=n*'N' ##### N is a missing data
			elif len(snp) > 1: ##### heterozygous sites	
				fastaseq[sample]+=n*'N'	##### N is missing data
			else:######## add gap in high quality variations and consider the gaps valids
				if snp[0]=='0':
					fastaseq[sample]+=ref+add_ref
				elif snp[0]=='1':
					fastaseq[sample]+=alt+add_alt
				elif snp[0]=='.':
					fastaseq[sample]+=n*'N'

os.system('rm cds*')
for i in fastaseq.keys():##### conversion of names with genotypes here
#	n=i.split('_')[1]
	print('>%s\n%s' % (i, fastaseq[i]), file=open(fasta, 'a'))	
					
	
