## Create Rqtl plot from data
from __future__ import print_function
import sys

vcf = 'SNPs_for_map.recode.vcf'
outfile = 'LimaSNPs.csv'

SNPcall = []
markerID = []
#chrom = [] NO CHROMOSOMES

#genoD={'0/0':'A', '1/1':'B', '0/1': '-', './.': '-'}
x=1
for line in open(vcf):
        if 'CHROM' in line:
                line=line.strip().split('\t')
                SNPcall.append(line[9:])
		p1p = line.index('UC92')
		p2p = line.index('UCHaskell')
        elif '#' not in line:
                line =line.strip().split('\t')
             	genocall=[i.split(':')[0] for i in line[9:]]
		p1 = line[p1p].split(':')[0]
		p2 = line[p2p].split(':')[0]
		if p1!=p2:
			print(p1,p2)
			markerID.append('_'.join(line[:2]))
			genoD={p1:'A', p2:'B', '0/1':'-', './.':'-'} # dinamically create dictionary based on the calls of the 2 parents
                	SNPcall.append([genoD[i] for i in genocall])

print('ID,', ','.join(markerID), sep='', file=open(outfile, 'a'))
print(',', ','.join('1'*len(markerID)), sep='', file=open(outfile, 'a'))

for item in range(len(SNPcall[0])):
        s=[col[item] for col in SNPcall]
        print(','.join(s), file=open(outfile, 'a'))
