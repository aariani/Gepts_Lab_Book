#! /usr/bin/env python
### Script for aligning reads
from __future__ import print_function
import argparse
import glob
from subprocess import call

parser=argparse.ArgumentParser(prog='BWAAln', description='Align fastq file to reference genome')
parser.add_argument('-i', '--input', dest='fastqreads', help='The folder with your fastq reads')
parser.add_argument('-r', '--referencegenome', dest='ref',
                    help='The reference genome file')
parser.add_argument('-o', '--output', dest='outfolder', help='The output folder')
parser.add_argument('-q', '--qual', dest='qual', type=int, help='Minimum mapping quality')
arg=parser.parse_args()

## Check if there is some parameter missing
if 'None' in str(arg):
        parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

## sample params
fastq=arg.fastqreads.split('/')[0]
ref=arg.ref
outfolder=arg.outfolder
qual=arg.qual
reads=glob.glob('%s/*' % fastq)

for i in reads:
	print('Start aligning sample %s' % i)
	sam_out=i.split('.')[0].split('/')[1]+'_aln.sam'
	cmd='bwa mem -t 4 %s %s > %s' % (ref, i, sam_out)
	print(cmd)
	call(cmd, shell=True)

print('remove multimapping reads')
aln=glob.glob('*sam')
for i in aln:
        name=i.split('.')[0] + '_sort'
        cmd='samtools view -bShq %s %s | samtools sort - %s' % (qual, i, name)
        print(cmd)
        call(cmd, shell=True)

call('rm *sam', shell=True)

bam_files=glob.glob('*_sort.bam')
for i in bam_files:
        cmd='samtools index %s' % i
        print(cmd)
        call(cmd, shell=True)

call('mkdir %s' % outfolder, shell=True)
call('mv *bam* %s' % outfolder, shell=True)



