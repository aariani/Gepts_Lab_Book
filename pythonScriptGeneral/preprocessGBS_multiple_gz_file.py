#! /usr/bin/env python
##### 
##### Script for cleaning and demultiplexing GBS reads
##### implemented for restriction enzyme with a fixed recognition sequence: default == CATG (CviaII)
##### try to implement the script for handling multiple fastq.gz files at once and output concatenated files with the same barcode
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
parser.add_option('-s', '--restriction_enzyme_site', dest='site', default='CATG', help='your Restriction Enzyme recognition site (default CATG)')
parser.add_option('-l', '--min_length', dest='minlen', type='int', default='30', help='Minimum length of reads after quality trimming (default 30)')
parser.add_option('-q', '--min_quality', dest='minQ', type='int', default='20', help='Mean minimum quality (in a sliding window of 5bp) for trimming reads, assumed Sanger quality (Illumina 1.9+, default 20)')

(opt, args) = parser.parse_args()

f=opt.reads
bc_file=opt.bc
out_folder=opt.good_reads
site=opt.site
minlen=opt.minlen
minQ=opt.minQ
if 'None' in str(opt):
	parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')	
############### Clipping chimeras
def clip_chimera_adapter(s, qual, site, cont):
	if site in s:
		p=s.index(site)
		s=s[:p]+'\n'
		qual=qual[:p] +'\n'
	elif cont in s:
		p=s.index(cont)
		s=s[:p]+'\n'
		qual=qual[:p]+'\n'
	return s, qual

############### Low quality trimming
def trim_qual(seq, qual, mq):
	n_qual=[ord(x)-33 for x in qual[:-1]]
	for i in range(len(n_qual)-4): ###### trim with a sliding window of 5, the 4 at the end is for avoid error
		clip=n_qual[i:i+5]
		m=float(sum(clip)/5)
		if m <= mq:
			qual=qual[:i]+'\n'
			seq=seq[:i]+'\n'
			break
	return seq, qual	

########################################################
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
all_reads=getoutput('ls %s/*.gz' % f).split()
print(all_reads)
for i in all_reads:
	print('start preprocessing %s' % i)
	read_f=gzip.open(i, 'rb')
	while True:
		name=read_f.readline()
		if name=='':
			break
		seq=read_f.readline()
		plus=read_f.readline()
		qual=read_f.readline()
		if '1:Y:0' not in name:## process only good reads
			frag=seq[:l]
			index=[i for i in barcode_d.keys() if i == frag[:len(i)]]
			if len(index)==0:
				a=gzip.open('unmatched.fq.gz', 'ab')
				a.write(name+seq+plus+qual)
				a.close()
        		else:   
				seq=seq[len(index[0]):]
				qual=qual[len(index[0]):]
				seq,qual=clip_chimera_adapter(seq, qual, site, contaminant)
				seq, qual=trim_qual(seq, qual, minQ)
				if len(seq.strip())>=minlen:			
					if seq[:3]==firstbases:			
						sample[barcode_d[index[0]]]+=1
						a=gzip.open(barcode_d[index[0]], 'ab')
						a.write(name+seq+plus+qual)
						a.close()
	
	read_f.close()
#out_f.close()
os.system('mkdir %s' % out_folder)
os.system('mv *gz %s' % out_folder)
for i in sample.keys():
	print(i, sample[i], sep='\t', file=open('stats_file2align.txt', 'a'))

