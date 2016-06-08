##### script for post process reads and call snp
##### reads should be in sam format
##### The reads will be filtered (-q 10 ) and sorted. 
import sys, os
from commands import getoutput

aln=getoutput('ls *sam').split('\n')
for i in aln:
	name=i.split('.')[0] + '_sort'
	cmd='samtools view -bShq 10 %s | samtools sort - %s' % (i, name)
	print(cmd)
	os.system(cmd)

######################## creates a bunch of files that finish with sort.bam


	
