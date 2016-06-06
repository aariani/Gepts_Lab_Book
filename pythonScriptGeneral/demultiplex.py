#! /usr/bin/env python
##### 
##### Script for cleaning and demultiplexing GBS reads
##### implemented for restriction enzyme with a fixed recognition sequence: default == CATG (CviaII)
##### try to implement the script for handling multiple fastq.gz. files at once and output concatenated files with the same barcode
##### 
##### 
from __future__ import print_function
import gzip
import os
from optparse import OptionParser
from commands import getoutput

parser=OptionParser( description='Script for cleaning and demultiplex reads of a GBS experiment, keeping the fastq file in a compressed format')
parser.add_option('-i', '--input', dest='reads', help='your folder with the raw fastq files')
parser.add_option('-o', '--output', dest='good_reads', help='the output folder where put all the cleaned fastq files')
parser.add_option('-b', '--barcode_file', dest='bc', help='Your barcode file formatted for sabre (https://github.com/najoshi/sabre)')

(opt, args) = parser.parse_args()

f=opt.reads
bc_file=opt.bc
out_folder=opt.good_reads

if 'None' in str(opt):
	parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')	
############### Clipping chimeras
#### establishment of initial parameters

firstbases='ATG'
contaminant='AGATCGG'
#suffix=f.split('.')[0]+'_'

barcode_d={}
for line in open(opt.bc):
        line=line.strip().split()
        barcode_d[line[0]]=line[1]+'.gz'
        gzip.open(line[1]+'.gz', 'wb')

l=max([len(i) for i in barcode_d.keys()])
sample={k:0 for k in barcode_d.values()}
########################################################
##### starting preprocess
read_f=gzip.open(f, 'rb')
while True:
	name=read_f.readline()
	if name=='':
		break
	seq=read_f.readline()
	plus=read_f.readline()
	qual=read_f.readline()
	frag=seq[:l]
	index=[i for i in barcode_d.keys() if i == frag[:len(i)]]
	if len(index)==0:
		a=gzip.open('unmatched.fq.gz', 'ab')
		a.write(name+seq+plus+qual)
		a.close()
       	else:   
		seq=seq[len(index[0]):]
		qual=qual[len(index[0]):]
		a=gzip.open(barcode_d[index[0]], 'ab')
		a.write(name+seq+plus+qual)
		a.close()
	
read_f.close()
#out_f.close()
os.system('mkdir %s' % out_folder)
os.system('mv *gz %s' % out_folder)
for i in sample.keys():
	print(i, sample[i], sep='\t', file=open('stats_file2align.txt', 'a'))

