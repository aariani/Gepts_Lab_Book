#! /usr/bin/env python
## Script for converting vcf file to Rqtl
from __future__ import print_function
import argparse 
import glob

parser = argparse.ArgumentParser(prog='VCF2Rqtl',
                                 description='Convert VCF file to Rqtl')
parser.add_argument('-i', '--input', dest='vcf',
                    help='Your VCF input file')
parser.add_argument('--p1', dest='p1',
                    help='The ID of the first parent (as in the input file')
parser.add_argument('--p2', dest='p2',
                    help='The ID of the second parent (as in the input file')
parser.add_argument('-o', '--output', dest='outfile',
                    help='The name of your output file')
arg=parser.parse_args()

if 'None' in str(arg):
    parser.error('Input parameter missing!! Please check your command'+
                 'parameters with --h or --help')

vcf = arg.vcf
outfile=arg.outfile
p1=arg.p1
p2=arg.p2


good_calls = ['0/0', '1/1']
SNPcall = []
markerID = []

for line in open(vcf):
        if 'CHROM' in line:
                line=line.strip().split('\t')
                SNPcall.append(line[9:])
		p1p = line.index(p1)
		p2p = line.index(p2)
        elif '#' not in line:
                line =line.strip().split('\t')
                genocall=[i.split(':')[0] for i in line[9:]]
		g1 = line[p1p].split(':')[0]
		g2 = line[p2p].split(':')[0]
                if g1 in good_calls and g2 in good_calls:
                    if g1!=g2:
                        markerID.append('_'.join(line[:2]))
                        # dinamically create dictionary based on the calls of
                        # the 2 parents
                        genoD={g1:'A', g2:'B', '0/1':'-', './.':'-'}
                        SNPcall.append([genoD[i] for i in genocall])

print('ID,', ','.join(markerID), sep='', file=open(outfile, 'a'))
print(',', ','.join('1'*len(markerID)), sep='', file=open(outfile, 'a'))

for item in range(len(SNPcall[0])):
        s=[col[item] for col in SNPcall]
        print(','.join(s), file=open(outfile, 'a'))
