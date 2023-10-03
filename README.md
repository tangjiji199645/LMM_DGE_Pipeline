# LMM DGE Pipeline

## Setup GEMMA

1. For Mac or Windows, Using docker to run GEMMA, https://github.com/genetics-statistics/GEMMA. After adding the GEMMA image into docker, open terminal to run this following command to run GEMMA:
```
 docker run -w /run -v ${local_path for your files}:/run ed5bf7499691 gemma
```

2. For Linux or HPC, download the binary format from https://github.com/genetics-statistics/GEMMA.
```
chmod u+x gemma 
./gemma 
```

## 1.Normalization and prepare file for GEMMA

Required file <br>
Raw read counts: sample_raw_reads.txt <br>
Covariate matrix: cov_matrix.txt <br>

For read counts file, the first three column is gene id, the second and third column is allele types (Ignore this here, type all A/T/C/G in one column). 

For the covariate matrix, intercept is manually required.

Use DESeq2 to normalize the raw read counts, for details, see normalization.R

## 2.LMM test using GEMMA 

For installment detail, please see https://github.com/genetics-statistics/GEMMA

Required file <br>
gzip read counts: normalized_reads.txt.gz <br>
Phenotype: phenotype.txt <br>
Covariate matrix: cov_bim.txt <br>

Get Kinship matrix <br>
gemma -g normalized_reads.txt.gz -p phenotype.txt -c cov_bim.txt -gk 2 -notsnp -o cov_mat

LMM <br>
gemma -g normalized_reads.txt.gz -p phenotype.txt -k output/cov_mat.sXX.txt -c cov_bim.txt -lmm 4 -notsnp -o output

## 3.Visualization 
Required file <br>
sample data: sample_data.txt <br>
qq plot, manhattan plot and volcano_plot

