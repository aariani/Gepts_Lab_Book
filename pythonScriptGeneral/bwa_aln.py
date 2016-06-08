####### script for automate alignment ousin bwa aligner and the index already prepared for 
####### tassel
####### give the suffix of the bam files and thats it!!!
import sys, os
from commands import getoutput

suffix= sys.argv[1]
demultiplex_files=getoutput('ls %s*' % suffix).split('\n')
for i in demultiplex_files:
	print(i)
	sai_out=i.split('.')[0] + '_aln.sai'
	cmd1='bwa aln -n 0.1 -o 2 -e 6 -l 20 -k 3 -t 4 ../../index/Pvulgaris_218.fa %s > %s' % (i, sai_out) 
	print(cmd1)
	os.system(cmd1)
	sam_out=sai_out.split('.')[0] + '.sam'
	cmd2='bwa samse ../../index/Pvulgaris_218.fa %s %s > %s' % (sai_out, i, sam_out)
	print(cmd2)
	os.system(cmd2)

