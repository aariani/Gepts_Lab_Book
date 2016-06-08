## Pipeline for processing a brazilian collection of common bean

The `masterpipe.sh` script calls the different other script

For simplicity I usually divide the TASSEL-GBS pipeline in 4 Step:

1.	**Part1**: Part before tag alignment using `tassel_gbs_part1.sh`

2.	**Alignment**: Alignment step using `tassel_aln_step.py`

3.	**Part2**: SNP calling final part using `tassel_gbs_part2.sh`

4.	**Imputation**: Imputation part with `snpCall_Impute.sh`

=======================================================================

