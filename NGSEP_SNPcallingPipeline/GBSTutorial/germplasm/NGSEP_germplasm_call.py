#! /usr/bin/env python
# Script for bi-parental SNPs calling

from __future__ import print_function
import argparse
import glob
from subprocess import call

parser=argparse.ArgumentParser(prog='SNPs_calling',
                               description='Call SNPs on biparental populations with NGSEP')
parser.add_argument('-i', '--input', dest='samfiles',
                    help='The folder with your alignment files')
parser.add_argument('--NGSEP', dest='ngsep',
                    help='The path of NGSEP executable file')
parser.add_argument('-r', '--referencegenome', dest='ref',
                    help='The reference genome file')

arg=parser.parse_args()

 ## Check if there is some parameter missing
if 'None' in str(arg):
    parser.error('Input parameter missing!! Please check your command'+
                 'parameters with --h or --help')
# sample params
samfiles=arg.samfiles.split('/')[0]
ngsep=arg.ngsep
ref=arg.ref
p1=arg.p1
p2=arg.p2

# first call
all_aln = glob.glob('%s/*bam' % samfiles)

for i in all_aln:
    newID = i.split('/')[1].split('_')[0]
    cmd='java -jar ' + ngsep + ' FindVariants -h 0.0001 -noRep -noRD -noRP '\
            '-maxBaseQS 30 -minQuality 40 -maxAlnsPerStartPos 100 '\
            '-sampleId ' + newID + ' -ignoreXS %s %s NGSEP_%s' % (ref, i, newID)
    call(cmd, shell=True)

# transfer SNPs calls in mapping folder
call('mkdir mapping', shell=True)
call('mv NGSEP*vcf mapping', shell=True)

# generate HQ sites
call(r"grep '>' "+ ref + r"|cut -f 2 -d '>' > seqNames.txt", shell=True)

# merge variants
call('java -jar ' + ngsep+' MergeVariants seqNames.txt variants_list.vcf '\
     'mapping/NGSEP*vcf', shell=True)

for i in all_aln:
    newID = i.split('/')[1].split('_')[0]
    cmd = 'java -jar '+ngsep+' FindVariants -knownVariants variants_list.vcf '\
            '-h 0.0001 -noRep -noRD -noRP -maxBaseQS 30 '\
            '-maxAlnsPerStartPos 100 -sampleId '+newID+\
            ' -ignoreXS %s %s NGSEP_HQ_%s' % (ref, i, newID)
    call(cmd, shell=True)

call('mkdir genotyping', shell=True)
call('mv NGSEP_HQ*vcf genotyping', shell=True)

cmd='java -jar '+ngsep+' MergeVCF seqNames.txt '\
        'genotyping/NGSEP_HQ*vcf > finalSNPs_all.vcf'
call(cmd, shell=True)
