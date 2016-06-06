from __future__ import print_function
import sys

mapfile=sys.argv[1]	## map file change col[0] in numerical value
pedfile=sys.argv[2]	### ped file change col[5] in 1
c=sys.argv[3] ##chrom name

x=0
d_c={}
### dictionary with chr in numeric formula
for line in open(c):
	if 'Chr' in line:
		x+=1
		d_c[line.strip()]=x

for line in open(mapfile):
	line=line.strip().split('\t')
	cname=line[1].split(':')[0]
	print(d_c[cname], line[1], line[2], line[3], sep='\t', file=open('final.map', 'a'))

for line in open(pedfile):
	line=line.strip().split('\t')
	line[5]='1'
	print(line[:7])
	print('\t'.join(line), file=open('final.ped', 'a')) 


