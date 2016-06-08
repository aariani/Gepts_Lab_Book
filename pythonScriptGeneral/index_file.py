##### subscript for indexing sam files
import os
from commands import getoutput
bam_files=getoutput('ls *recalibrated.bam').split('\n')
for i in bam_files:
	cmd='samtools index %s' % i
	print(cmd)
	os.system(cmd)
