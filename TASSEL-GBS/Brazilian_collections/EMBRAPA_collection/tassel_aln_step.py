#! /usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import os


parser=OptionParser(description='This is the second part of the customized TASSEL-GBS pipeline. This script aligns the GBS tags fonverted in fastq from the first parte of TASSEL-GBS pypeline')
parser.add_option('-i', '--input-fastq', dest='reads', help='Your fastq file')
parser.add_option('-r', '--ref', dest='ref', help='Your reference file')
parser.add_option('-q', '--map-qual', dest='mapqual', type='int', default='10', help='minimum read mapping quality, for avoid multiple mapping reads (default 10)')
parser.add_option('-t', '--tassel-ref-transform', dest='tasselref', action='store_true', default=False, help='Do you want to transform your Fasta reference header for Tassel? Tassel use numeric index in fasta header. If selected the new genome will be indexed for bwa automatically and a new_folder named TASSEL_INDEX will be created in the current directory')
parser.add_option('-I', '--index-ref', dest='bwa_index', action='store_true', default=False, help='Do you want to index your reference? You just need to do it once after your reference has been transformed for taxel')
(opt, args)=parser.parse_args()

#### inizialization of variables

fq=opt.reads
ref=opt.ref
mapQ=opt.mapqual
tasselref=opt.tasselref
bwa_index=opt.bwa_index

if 'None' in str(opt):
        parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

def aln_n_compress(f, r, q):
	out_sai=f.split('.')[0] +'.sai'
	out_sam=f.split('.')[0] +'.sam'
###### with bwa mem didnt works!!! modify using bwa aln
	os.system('bwa aln -t 4 %s %s > %s' % (r, f, out_sai))
	os.system('bwa samse %s %s %s > %s' % (r, out_sai, f, out_sam)) 
	out_uniq=out_sam.split('.')[0] +'_uniq.sam'
	os.system('samtools view -Shq %s %s > %s' % (q, out_sam, out_uniq)) 	
	

##### Starting pipeline with indexing and all the other stuffs
if tasselref:
	print('Starting modifying header of fasta reference sequence, the new genome will be indexed for bwa')
	bwa_index=True
	os.system('mkdir TASSEL_INDEX')
	os.system('cp %s TASSEL_INDEX/' % ref)
	new_ref='tassel_'+ref.split('/')[-1]
	#### make new ref
	x=0
	for line in open('TASSEL_INDEX/%s' % ref.split('/')[-1]):
		if '>' in line:
			x+=1
			print('>%s' % x, file=open('TASSEL_INDEX/%s' % new_ref, 'a'))
		else:
			print(line.strip(), file=open('TASSEL_INDEX/%s' % new_ref, 'a'))
	os.system('bwa index -a bwtsw TASSEL_INDEX/%s' % new_ref)
	aln_n_compress(fq, 'TASSEL_INDEX/%s' % new_ref, mapQ)		

else:
	aln_n_compress(fq, ref, mapQ)

