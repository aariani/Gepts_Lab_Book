#### Master pipeline for including all the steps of TASSEL GBS

### Part1: initial filtering, tag count and conversion to fastq file for alignment steps.
### This part will produce a tagCounts directory weith the counts of all tags, a mergedTagCounts directory with the merged counts of the 2 libraries
### and an aln folder containing the collapsed tags in fastq format (GBSTags.fq file)
./tassel_gbs_part1.sh

### Alignment Step: this script will align the collapsed tags using bwa aln and bwa samse alignment. The reference genome has been already modified for tassel (with numerical chromosomes) and indexed with bwa index -a bwtsw option. Btw there is an option for indexing directly the genome from the script, and it will create a TASSEL_INDEX directory in the current directory
### The script will align all the tags and then will filter the uniquely aligned tags (mapping quality > 10). The output is stored in the al folder
### There will be a GBSTags.sam file and the GBSTags_uniq.sam with the uniquely mapped reads
 
./tassel_aln_step.py -i aln/GBSTags.fq -r ~/pvulgaris/assembly/tassel/tasselPvul.fa 

### Part2: Create topm and TBT files for subsequent SNPs calling.
### It should produce a hapmap file with the filtered SNPs by chromosome. The hapmap files could be concatenated and then converter to VCFtools as well
./tassel_gbs_part2.sh

./impute.sh 
