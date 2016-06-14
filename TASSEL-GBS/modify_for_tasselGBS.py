### Modify your Reference File for use it with TASSEL-GBS
### TASSEL-GBS accept only numerical chromosomes in the pipeline
### This script will convert the names of chromosomes/scaffolds of your reference file to numerical values
### In numerical one. It will also output the file 'aliases.txt'  with the name conversion of the chromosomes/scaffolds name in the original fasta file

from __future__ import print_function
import sys


try:
	infile=sys.argv[1]
	outfile=sys.argv[2]
	x=0
	for line in open(infile):
		if '>' in line:
			x+=1
			print('>%s' % x, file=open(outfile, 'a'))
			print(line.strip()[1:], x, file=open('aliases.txt', 'a'))
		else:
			print(line.strip(), file=open(outfile, 'a'))

except:
	print('\t\t Script for modify your reference genome for TASSEL-GBS pipeline\n\t\tUSAGE: python modify_for_tasselGBS.py REFERENCE_GENOME OUTPUT_FILE
