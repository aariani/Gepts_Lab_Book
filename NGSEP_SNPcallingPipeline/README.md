## Scripts used for NGSEP SNP calling from command line

The NGSEP pipeline was developed at CIAT. More informations on the pipeline and tutorials could be found in the [published paper](http://www.ncbi.nlm.nih.gov/pubmed/24413664) or in the [download page](https://sourceforge.net/projects/ngsep/)

There is a `scripts` folder with the scripts used in the analysis in the lab

I included also a step-by-step tutorial in the `GBSTutorial` folder where there is a detailed protocol. 

I use NGSEP just for SNPs calling, while I preprocess the reads with a [custom pipeline](https://github.com/aariani/GBSprep). The alignment is performed with BWA, and post-processed with samtools for keeping only high quality aligned reads

