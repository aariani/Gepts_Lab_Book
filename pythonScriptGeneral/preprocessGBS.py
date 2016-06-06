#! /usr/bin/env python
##### 
##### Script for cleaning and demultiplexing GBS reads
##### implemented for restriction enzyme with a fixed recognition sequence: default == CATG (CviaII)
##### 
##### required the installation of sickle in the system for demultiplex
##### 
from __future__ import print_function
import os
from optparse import OptionParser
from commands import getoutput

parser=OptionParser( description='Script for cleaning and demultiplexing reads of a GBS experiment. It is a wrapper for sabre, so you need to install it in your system')
parser.add_option('-i', '--input', dest='reads', help='your raw read file in fastq format')
parser.add_option('-b', '--barcode_file', dest='bc', help='Your barcode file formatted for sabre (https://github.com/najoshi/sabre)')
parser.add_option('-s', '--restriction_enzyme_site', dest='site', default='CATG', help='your Restriction Enzyme recognition site (default CATG)')
parser.add_option('-l', '--min_length', dest='minlen', type='int', default='30', help='Minimum length of reads after quality trimming (default 30)')
parser.add_option('-q', '--min_quality', dest='minQ', type='int', default='20', help='Mean minimum quality (in a sliding window of 5bp) for trimming reads, assumed Sanger quality (Illumina 1.9+, default 20)')

(opt, args) = parser.parse_args()

f=opt.reads
bc_file=opt.bc
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
################################ evaluate results
def evaluate_res(infile, o):
	x=0
	for line in open(infile): 
		x+=1
	final=x/4
	print(infile, final, sep='\t', file=open(o, 'a'))
	
########################################################
#### establishment of initial parameters

firstbases='ATG'
contaminant='AGATCGG'
########################################################
print('Starting clipping adapter,filtering and trimmin chimera')
read_f=open(f)
out_f_name='clean_'+f.split('.')[0] +'.fastq'
out_f=open(out_f_name, 'w')
while True:
	name=read_f.readline()
	if name=='':
		break
	seq=read_f.readline()
	plus=read_f.readline()
	qual=read_f.readline()
	seq,qual=clip_chimera_adapter(seq, qual, site, contaminant)
	seq, qual=trim_qual(seq, qual, minQ)
	if len(seq.strip())>=minlen:
		if '1:Y:0' not in name:## remove filtered reads
			out_f.write(name+seq+plus+qual)
read_f.close()
out_f.close()

evaluate_res(out_f_name, 'stats_first_filtering.txt')
os.system('mkdir raw_read')
os.system('mv %s raw_read/' % f)

#########################################################
print('starting demultiplexing samples')
print(out_f_name)
dem='sabre se -f %s -b %s -u unknown.fq' % (out_f_name, bc_file)
print(dem)
os.system(dem)

sample=[]
for line in open(bc_file):
	sample.append(line.strip().split()[1])
print(sample)

print('Evaluating reads number after demultiplexing')
for i in sample:
	evaluate_res(i, 'stats_demultiplex.txt')

to_align=[]
print('Keep only reads starting with RE overhang sequence and >= 30bp after demultiplexing')
for i in sample:
	dem=open(i)
	dem_outid='final_'+ i
	to_align.append(dem_outid)
	dem_out=open(dem_outid, 'w')
	while True:
		name=dem.readline()
		if name=='':break
		seq=dem.readline()
		plus=dem.readline()
		qual=dem.readline()
		if seq[:3]==firstbases:
			if len(seq.strip())>=minlen:
				dem_out.write(name+seq+plus+qual)
	
	dem.close()
	dem_out.close()

for i in to_align:
	evaluate_res(i, 'stats_file2align.txt')

for i in sample:
	os.system('rm %s' %i)

#######################
