#! /usr/bin/env python
###### script for analyzing read number and extract only the good genotypes
###### Reads are already cleaned, just select the samples that has a read count >= x percentile of reads distribution
###### it's the second step of GBS pipeline
###### need a file containing the number of reads for each demultiplexed sample
###### its a wrapper for samtools and bwa so you need this program installed in your computer

from __future__ import print_function
import os
import numpy as np
from commands import getoutput
from optparse import OptionParser

parser=OptionParser(description='Script for evaluating read count and align genotypes that has enough read counts after cleaning steps')
parser.add_option('-c', '--read_count_file', dest='read_count', default='stats_file2align.txt', help='file with read counts for each sample after preprocessing step (default stats_file2align.txt)')
parser.add_option('-l', '--lower_limit', dest='lowLimit', type='int', default='10', help='lower limit of percentile for read count distibution, sample under this limit will be discarded from the alignment (default 10th percentile)')
parser.add_option('-r', '--ref', dest='ref', help='the path for your reference genome file')
parser.add_option('-I', '--index_genome', dest='bwa_index', action='store_true', default=False, help='Do you want to index your genome?? (Just do it once cause it is time consuming') 
(opt, args)=parser.parse_args()

read_count=opt.read_count
lowLimit=opt.lowLimit
ref=opt.ref
BWAindex=opt.bwa_index
#### report error if no input args
if 'None' in str(opt):
	parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

##### make count dictionary
count_dict={}
for line in open(read_count):
	if 'blank' not in line.lower():
		line=line.strip().split('\t')
		count_dict[line[0]]=int(line[1])

l=[count_dict[i] for i in count_dict.keys()]
l=np.percentile(l, lowLimit)

valid_geno=[]
for i in count_dict.keys():
	if count_dict[i]>=l:
		valid_geno.append(i)
for i in valid_geno:
	print(i, count_dict[i], sep='\t', file=open('stats_aligned_genotypes.txt', 'a'))

print(valid_geno)
if BWAindex:
        print('Starting genome indexing')
        os.system('bwa index -a bwtsw %s' % ref)
        os.sytem('samtools faidx %s' % ref)

for i in valid_geno:
#####opt for bwa mem
        sam_out=i.split('.')[0]+'_aln.sam'
        cmd='bwa mem -t 4 %s %s > %s' % (ref, i, sam_out)
        print(cmd)
        os.system(cmd)

print('remove multimapping reads')
aln=getoutput('ls *sam').split('\n')
for i in aln:
        name=i.split('.')[0] + '_sort'
        cmd='samtools view -bShq 10 %s | samtools sort - %s' % (i, name)
        print(cmd)
        os.system(cmd)

#os.system('rm *sam')

bam_files=getoutput('ls *_sort.bam').split('\n')
for i in bam_files:
        cmd='samtools index %s' % i
        print(cmd)
        os.system(cmd)

###snp calling
print('start snp calling')
os.system('samtools mpileup -uDgf %s *_sort.bam|bcftools view -vcg -d 0.3 - > all_snp_raw.vcf' % ref)


	
