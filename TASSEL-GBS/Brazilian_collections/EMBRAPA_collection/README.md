## Pipeline for processing a brazilian collection of common bean

The `masterpipe.sh` script calls the different other script

For simplicity I usually divide the TASSEL-GBS pipeline in 4 Step:

1.	**Part1**: Part before tag alignment using `tassel_gbs_part1.sh`

2.	**Alignment**: Alignment step using `tassel_aln_step.py`

3.	**Part2**: SNP calling final part using `tassel_gbs_part2.sh`

4.	**Imputation**: Imputation part with `impute.sh`

=======================================================================

The collection analyzed in this study is from the EMBRAPA germplasm of Brazilian Landraces

The genotypes analyzed are the same as the manuscript:

* [Microsatellite diversity and genetic structure among common bean (Phaseolus vulgaris L.) landraces in Brazil, a secondary center of diversity](http://link.springer.com/article/10.1007/s00122-010-1350-5/fulltext.html)
