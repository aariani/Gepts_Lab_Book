#! /usr/bin/env python
### Script for aligning reads
from __future__ import print_function
import argparse, os
from commands import getoutput

parser=argparse.ArgumentParser(prog='BWAAln', description='Align fastq file to reference genome')
parser.add_argument('-i', '--input', dest='fastqreads', help='The folder with your fastq reads')
parser.add_argument('-r', '--referenceGenome', dest='ref', help='yout reference genome')
parser.add_argument('-o', '--output', dest='outfolder', help='The output folder')

arg=parser.parse_args()

## Check if there is some parameter missing
if 'None' in str(arg):
        parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

## sample params
fastq=arg.fastqreads.split('/')[0]
ref=arg.ref
outfolder=arg.outfolder

reads=getoutput('ls %s/*' % fastq).split()

for i in reads:
	print('Start aligning sample %s' % i)
	sam_out=i.split('.')[0].split('/')[1]+'_aln.sam'
	cmd='bwa mem -t 4 %s %s > %s' % (ref, i, sam_out)
	print(cmd)
	os.system(cmd)

print('remove multimapping reads')
aln=getoutput('ls *sam').split()
for i in aln:
        name=i.split('.')[0] + '_sort'
        cmd='samtools view -bShq 10 %s | samtools sort - %s' % (i, name)
        print(cmd)
        os.system(cmd)

os.system('rm *sam')

bam_files=getoutput('ls *_sort.bam').split('\n')
for i in bam_files:
        cmd='samtools index %s' % i
        print(cmd)
        os.system(cmd)
                          
os.system('mkdir %s' % outfolder)
os.system('mv *bam* %s' % outfolder) 



