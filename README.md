# LMM DGE Pipeline


## 1.Normalization and prepare file for GEMMA

Required file <br>
Raw read counts: sample_raw_reads.txt <br>
Covariate matrix: cov_matrix.txt <br>

See normalization.R

## 2.LMM test using GEMMA 

For installment detail please see https://github.com/genetics-statistics/GEMMA

Required file <br>
gzip read counts: normalized_reads.txt.gz <br>
Phenotype: phenotype.txt <br>
Covariate matrix: cov_bim.txt <br>

Get Kinship matrix <br>
gemma -g normalized_reads.txt.gz -p phenotype.txt -c cov_bim.txt -gk 2 -notsnp -o cov_mat

LMM <br>
gemma -g normalized_reads.txt.gz -p phenotype.txt -k output/cov_mat.sXX.txt -c cov_bim.txt -lmm 4 -notsnp -o output

