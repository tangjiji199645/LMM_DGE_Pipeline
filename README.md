# LMM DGE Pipeline


## 1.Normalization and prepare file for GEMMA

Required file
Raw read counts: sample_raw_reads.txt
Covariate matrix: cov_matrix.txt

See normalization.R

## 2.LMM test using GEMMA 

For details, please see https://github.com/genetics-statistics/GEMMA

Required file
gzip read counts: normalized_reads.txt.gz
Phenotype: phenotype.txt
Covariate matrix: cov_bim.txt

Get Kinship matrix 
gemma -g normalized_reads.txt.gz -p phenotype.txt -c cov_bim.txt -gk 2 -notsnp -o cov_mat

LMM 
gemma -g normalized_reads.txt.gz -p phenotype.txt -k output/cov_mat.sXX.txt -c cov_bim.txt -lmm 4 -notsnp -o output

