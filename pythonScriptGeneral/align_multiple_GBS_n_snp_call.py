#! /usr/bin/env python
#### script for aligning reads of a GBS experiment using multiple libraries.
#### does not take into account duplicate genotypes in different libraries
#### Put all the file to align in a single directory 
#### all the input is piped for reducing the memory input of the analysis
from __future__ import print_function
import os
from commands import getoutput
from optparse import OptionParser

parser=OptionParser(description='Script for aligning reads of a GBS experiment using bwa mem algorithm with default parameters. Implemented for reducing memory input of the analysis')
parser.add_option('-i', '--input', dest='reads', help='your folder with the preprocessed reads ready for being aligned')
parser.add_option('-o', '--output', dest='aln_folder', help='your output folder')
parser.add_option('-r', '--ref', dest='ref', help='the path for your reference genome file')
parser.add_option('-I', '--index_genome', dest='bwa_index', action='store_true', default=False, help='Do you want to index your genome?? (Just do it once cause it is time consuming')
parser.add_option('-q', '--min_quality', dest='minQ', type='int', default='10', help='Minimum read mapping quality (default 10)')
(opt, args)=parser.parse_args()

clean_reads=opt.reads
out_folder=opt.aln_folder
ref=opt.ref
minQ=opt.minQ
BWAindex=opt.bwa_index

#### report error if no input args
if 'None' in str(opt):
        parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

#####
if BWAindex:
        print('Starting genome indexing')
        os.system('bwa index -a bwtsw %s' % ref)
        os.sytem('samtools faidx %s' % ref)


all_reads=getoutput('ls %s*' % clean_reads).split()
print(all_reads)
for i in all_reads:
	print('Start processing %s' % i)
	bam_out='hq_aln_'+i.split('/')[-1].split('.')[0]
        cmd='bwa mem -t 4 %s %s | samtools view -bShq %s -|samtools sort - %s' % (ref, i, minQ, bam_out)
        print(cmd)
        os.system(cmd)

print('Start indexing the bam files')
bam_files=getoutput('ls *bam').split()
for i in bam_files:
	os.system('samtools index %s' % i)


###snp calling
#print('start snp calling')
#os.system('samtools mpileup -uDgf %s *.bam|bcftools view -vcg -d 0.3 - > all_snp_raw.vcf' % ref)

os.system('mkdir %s' % out_folder)
os.system('mv *.bam* %s' % out_folder)
	
 


