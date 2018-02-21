#! /usr/bin/env python3

from __future__ import print_function
import gzip
import glob
import argparse
from subprocess import call

#################################################################
def getBCInfo(bcfile):
### Get barcodes infos, open files accordingly to the Genotypes available in the barcode
### Return a barcode dictionary and the longer barcode (useful for further demultiplexing), and a dictionary used for demultiplexing stats
    barcode_d={}
    for line in open(bcfile):
        line=line.strip().split()
        barcode_d[line[0]]=line[1]+'.fq'
        open(line[1]+'.fq', 'w')
    l=max([len(i) for i in barcode_d.keys()])
    demInfo={k:0 for k in barcode_d.values()}
    return barcode_d, l, demInfo

###############################################################
def getBCindex(name,sequence,quality,barcode_d,l):
### Script for getting the barcode file index for the sequence.
### If find the barcode return the sequence and quality without the barcode sequence
    frag=sequence[:l] ## Extract the initial part of the sequence
    index=[i for i in barcode_d.keys() if i == frag[:len(i)]]
    if len(index)==0:
### Write non demultiplexed reads
        return False,False,False  ## index is not defined previously, so it will be just none, without any name associated
    else:
        index=index[0]
        sequence=sequence[len(index):]
        quality=quality[len(index):]
        return index,sequence,quality

################################################################

def clip_chimera_and_adapters(sequence, quality, REsite, adapter):
### Main function for simultaneously remove chimeras or sequencing errors (i.e. sequences that showed the RE site)
### and the adapter sequence.
### It will also do a sliding window analysis for filtering the sequences
    RE=[i in sequence for i in REsite]
    if True in RE:
        pos=[sequence.find(i) for i in REsite] ### Index RE site ## use find cause index raise an error in case 1 site is ok and the other not.
        pos=[i for i in pos if i >= 0] ### In case it did not find any modification in one of the sites it will raise a -1
        firstpos=min(pos) ## first occurrence
        sequence=sequence[:firstpos]+'\n'
        quality=quality[:firstpos] +'\n'
    elif adapter in sequence:
        pos=sequence.index(adapter)
        sequence=sequence[:pos]+'\n'
        quality=quality[:pos]+'\n'
    return sequence, quality

############### Low quality trimming

def trim_qual(sequence, quality, minQ, minlen):
    n_qual=[ord(x)-33 for x in quality[:-1]]
    for i in range(len(n_qual)-4): ###### trim with a sliding window of 5, the 4 at the end is for avoid error
        clip=n_qual[i:i+5]
        m=float(sum(clip)/5)
        if m <= minQ:
            quality=quality[:i]+'\n'
            sequence=sequence[:i]+'\n'
            break
    if len(sequence.strip()) >=minlen:
        return sequence,quality
    else:
        return False,False

###################
def check_overhang(sequence, quality, rem_sites, rmRErem):
#### script for checking the overhang sequences
    overhangLength=len(rem_sites[0])
    if sequence[:overhangLength] in rem_sites:
        if rmRErem:
            sequence=sequence[overhangLength:]
            quality=quality[overhangLength:]
        return sequence, quality
    else:
        return False,False

###################################################################

def process(sequence, quality, REsite, contaminant, minQ, minlen, rem_site, rmRErem):
    sequence, quality=clip_chimera_and_adapters(sequence, quality, REsite, contaminant)
    sequence, quality=trim_qual(sequence,quality,minQ,minlen)
    if sequence: ## check quality filtering
        sequence, quality=check_overhang(sequence, quality, rem_site, rmRErem)
		## Check overhang site
        if sequence:
            return sequence, quality
        else:
            return False,False
    else:
        return False,False


##################################################################
### Start script

parser=argparse.ArgumentParser(prog='GBSprep', description='Program for preprocessing GBS reads')
parser.add_argument('-i', '--input-folder', dest='reads', help='The folder with your raw reads in gz compressed format (REQUIRED)')
parser.add_argument('-o', '--output-folder', dest='clean_reads', help='The output folder for the final cleaned and demultiplexed reads (REQUIRED)')
parser.add_argument('-bc', '--barcode-file', dest='bc_file', help='The barcode file for your raw reads (REQUIRED)')
parser.add_argument('-s', '--restriction_enzyme_site', dest='REsites', help='The Restriction Enzyme (RE) recognition site (REQUIRED). If your RE has ambiguous nucleotide you should write all the possible sites separated by a comma (ex. For ApeKI you should use -s GCAGC,GCTGC)')
parser.add_argument('-SR', '--site-remnant', dest='RErem', help='The RE remnant site after the digestion (REQUIRED). Only reads with the remnant site after the barcode will be kept. For RE with ambiguous nucleotide you should write all the possible remnant sites separated by a comma (as in the -s parameter)')
parser.add_argument('-l', '--min-length', dest='minlen', type=int, default=30, help='Minimum length of reads after quality trimming and adapter/chimeras clipping (default: 30bp)')
parser.add_argument('-q', '--min-qual', dest='minQ',type=int, default=20, help='Mean minimum quality (in a sliding window of 5bp) for trimming reads, assumed Sanger quality (Illumina 1.8+, default: 20)')
parser.add_argument('-ad', '--adapter-contaminants', dest='contaminant', default='AGATCGG', help='The initial sequence of the adapter contaminant (default: AGATCGG)')
parser.add_argument('--remove-remnant-site', dest='rmRErem', action='store_true', default=False, help='Do you want to remove the RE remnant site from the final cleaned reads? (default: False)')
args=parser.parse_args()

## Check if there is some parameter missing
if 'None' in str(args):
    parser.error('Input parameter missing!! Please check your command line parameters with -h or --help')

## Assigning paramters to variables
reads=args.reads.split('/')[0]
clean_reads=args.clean_reads
bc=args.bc_file
REsite=args.REsites.split(',')  ## This will be a list
rem_site=args.RErem.split(',') ## And this too. Remember for the final filtering
minlen=args.minlen
minQ=args.minQ
contaminant=args.contaminant
rmRErem=args.rmRErem
### Get barcode info and open all the file
barcode_d,l,demInfo=getBCInfo(bc)

### Start preprocess the files
all_reads=glob.glob('%s/*.gz' % reads)
for i in all_reads:
    print('Start preprocessing %s file' % i)
    read_f=gzip.open(i, 'rt')
    while True:
        name=read_f.readline()
        if name=='': break  ## stop parsing the file at the end
        sequence=read_f.readline()
        plus=read_f.readline()
        quality=read_f.readline()
        if '1:Y:0' not in name:  ## keep only reads passing initial Illumina filtering
            index,sequence,quality=getBCindex(name,sequence,quality,barcode_d,l)
            if index:
                sequence, quality= process(sequence, quality, REsite, contaminant, minQ, minlen, rem_site, rmRErem)
                if sequence:
                    a=open(barcode_d[index], 'a')
                    a.write(name+sequence+plus+quality)
                    a.close()
                    demInfo[barcode_d[index]]+=1
    read_f.close()

call('mkdir %s' % clean_reads, shell=True)
call('mv *fq* %s' % clean_reads, shell=True)

print('Sample\tTotalReads', file=open('Demultiplexing_stats_%s.txt' % clean_reads, 'a'))
for i in demInfo.keys():
    print(i.split('.')[0], demInfo[i], sep='\t',
          file=open('Demultiplexing_stats_%s.txt' % clean_reads, 'a'))
