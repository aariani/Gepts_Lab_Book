#!/bin/sh

./GBSprep.py -i fastq -o demult -bc barcodesMxG.txt -s CATG -SR ATG -gz --remove-remnant-site

./bwa_aln.py -i demult -r ~/pvulgaris/assembly/bwa_index/Pvulgaris_218.fa -o aln
