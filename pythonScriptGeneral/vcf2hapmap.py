##### convert vcf file to hapmap format
##### GAPIT freak out if there's only one snp per chr, so I will need to merge the scaffold together using a walk
#####
from __future__ import print_function
import sys, random

def modify_alleles(r,v):
	while True:
		newR=random.choice(r)
		newV=random.choice(v)
		if newR!=newV:
			break
	return newR, newV
def convert_geno(genolist, r, v):
	dgeno={'0/0':r*2, '0/1': r+v, '1/1':v*2, './.':'NN'}
	g=[dgeno[i] for i in genolist]
	return g

vcf=sys.argv[1]
hmp=vcf.split('.')[0]+'.hmp.txt'

chrom=[]
chrCount=0
scaff_start=0
scaff_walk=0
scaff_name=''

header='rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\t'

for line in open(vcf):
	if 'CHROM' in line:
		line=line.strip().split('\t')
		geno=line[9:]
		geno='\t'.join(geno)
		print(header + geno, file=open(hmp, 'a'))
	elif '#' not in line:
		if 'Chr' in line:
			line=line.strip().split('\t')
			if line[0] not in chrom:
				chrom.append(line[0])
				chrCount+=1
			snpid='_'.join([line[0], line[1]])
			pos=line[1]
			ref, alt=line[3], line[4]
			print(snpid, ref, alt, sep='\t', file=open('varInfo.txt', 'a'))
			if len(ref) > 1 or len(alt)> 1:
				ref, alt=modify_alleles(ref, alt)
			alleles='/'.join([ref, alt])
			geno=[i.split(':')[0] for i in line[9:]]
			geno=convert_geno(geno, ref, alt)
			geno='\t'.join(geno)
			print(snpid, alleles, chrCount, pos, '\t'.join(['NA']*7), geno, sep='\t', file=open(hmp, 'a'))
#### part for scaffold walking
		if 'scaffold' in line:
                        line=line.strip().split('\t')
                        if scaff_name!=line[0]:
                                scaff_name=line[0]
                                scaff_start+=scaff_walk
				scaff_walk=int(line[1])-1
                	else:
				scaff_walk=int(line[1])-1       
			snpid='_'.join([line[0], line[1]])
			pos=scaff_start+ int(line[1])
                        ref, alt=line[3], line[4]
                        print(snpid, ref, alt, sep='\t', file=open('varInfo.txt', 'a'))
			if len(ref) > 1 or len(alt)> 1:
                                ref, alt=modify_alleles(ref, alt)
                       	alleles='/'.join([ref, alt])
                       	geno=[i.split(':')[0] for i in line[9:]]
                        geno=convert_geno(geno, ref, alt)
                        geno='\t'.join(geno)
                        print(snpid, alleles, 12, pos, '\t'.join(['NA']*7), geno, sep='\t', file=open(hmp, 'a'))
        
		


