#!/bin/sh
##### script for tassel4-gbs pipeline utilization
##### prepare a fastq folder with the files and the barcode files
##### This is only the first step! untill the production of the fake reads to alignment
#####
##### first step is to export the tassel pipeline path, so you don't have to write all the path everytime

TASSEL=<PATH_TO_MODIFIED_TASSEL4>
echo $PATH

mkdir tagCounts ####out for the first demultiplexing step

#### starting demultiplexing and isolating Tags
$TASSEL -Xmx12g -fork1 -FastqToTagCountPlugin -i fastq -k barcodes.csv -e cviaii -s 350000000 -o tagCounts -endPlugin -runfork1|tee logs_part1.txt

#### merge tag counts from different lanes, the -t tag extract the cnt in fastq
mkdir mergedTagCounts
$TASSEL -fork1 -MergeMultipleTagCountPlugin -i tagCounts -o mergedTagCounts/GBSTags.cnt -endPlugin -runfork1|tee logs_part2.txt

##### convert tag count to fastq for alignment
mkdir aln
$TASSEL -fork1 -TagCountToFastqPlugin -i mergedTagCounts/GBSTags.cnt -o aln/GBSTags.fq  -endPlugin -runfork1|tee logs_part3.txt



