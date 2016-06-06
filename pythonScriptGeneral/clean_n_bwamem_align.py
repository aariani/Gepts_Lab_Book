#! /usr/bin/env python
##### Script for alignment of GBS reads
##### Script for cleaning, demultiplexing and align reads against reference sequence
##### implemented for restriction enzyme with a fixed recognition sequence: default == CATG (CviaII)
##### 
##### required the installation of BWA, sabre and sickle in the system
##### 

from __future__ import print_function
import os, sys, math
from optparse import OptionParser
from commands import getoutput

parser=OptionParser( description='Script for cleaning, demultiplexing and aligning GBS reads (assumed SE) against a given reference genome. It is a wrapper for sickle, sabre and BWA, so you need to install them in your system')
parser.add_option('-i', '--input', dest='reads', help='your raw read file in fastq format')
parser.add_option('-r', '--reference', dest='ref', help='Reference genome file')
parser.add_option('-b', '--barcode_file', dest='bc', help='Your barcode file formatted for sabre (https://github.com/najoshi/sabre)')
parser.add_option('-s', '--restriction_enzyme_site', dest='site', default='CATG', help='your Restriction Enzyme recognition site')
parser.add_option('-l', '--min_length', dest='minlen', type='int', default='20', help='Minimum length of reads after quality trimming')
parser.add_option('-q', '--min_quality', dest='minQ', type='int', default='20', help='Mean minimum quality (in a sliding window of 5bp) for trimming reads, assumed Sanger quality (Illumina 1.9+)')
parser.add_option('-I', '--index_genome', dest='bwa_index', action='store_true', default=False, help='Do you want to index the genome file?')
(opt, args) = parser.parse_args()

f=opt.reads
ref=opt.ref
bc_file=opt.bc
site=opt.site
minlen=opt.minlen
minQ=opt.minQ
BWAindex=opt.bwa_index
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
########################################################
print('Starting clipping adapter and trimmin chimera')
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
	if len(seq)>=minlen:
		out_f.write(name+seq+plus+qual)
read_f.close()
out_f.close()

os.system('mkdir raw_read')
os.system('mv %s raw_read/' % f)

#########################################################
print('starting demultiplexing samples')
print(out_f_name)
dem='sabre se -f %s -b %s -u unknown -m 1' % (out_f_name, bc_file)###### BARCODE FILE NAME!!!!!!!!
print(dem)
os.system(dem)

sample=[]
for line in open(bc_file):
	sample.append(line.strip().split()[1])
print(sample)

to_align=[]
print('Keep only reads starting with RE overhang sequence')
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
			dem_out.write(name+seq+plus+qual)
	
	dem.close()
	dem_out.close()


for i in sample:
	os.system('rm %s' %i)

if BWAindex:
	print('Starting genome indexing')
	os.system('bwa index -a bwtsw %s' % ref)
	os.sytem('samtools faidx %s' % ref)

for i in to_align:
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
