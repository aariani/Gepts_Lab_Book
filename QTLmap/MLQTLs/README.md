Use of ML approacheas for QTL mapping analysis (in contrast to the other
approach since the L1 stability selection is not stable).

Rationale: Use of L1 regularization approaches for QTLs mapping by improving
current methodologies.

Took some idea from [this
manuscript](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5260053/)

Approach: Compare SVM/Lasso (for regression) or SVC/Logit for classifications

Steps:

1) Baggin by bootstrap with replacement X% (80%) of the data

2) Optimize C using CV of the bootstrap sample (try both SV and L1 general)

3) Fit the regressor and evaluate the weigths for each feature (SNPs)

4) Calculate for each SNPs how many time w > 0 over the sampling scheme (this will
refer as STABILITY SCORE (S))

5) Calculate w of each SNPs by np.sum(w) over the full sampling scheme

6) Give SNPs score by `np.sum(w)*S` (you should see how the distribution is)

7) Plot the score over the genome for identifying QTLs
