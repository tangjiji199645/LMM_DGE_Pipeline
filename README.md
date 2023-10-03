# LMM DGE Pipeline

If you use this pipeline for published work, please cite our paper:

Tang, S., Buchman, A.S., Wang, Y. et al. [Differential gene expression analysis based on linear mixed model corrects false positive inflation for studying quantitative traits]([https://github.com/user/repo/blob/branch/other_file.md](https://www.nature.com/articles/s41598-023-43686-7)). Sci Rep 13, 16570 (2023). https://doi.org/10.1038/s41598-023-43686-7



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

Required file : <br>
Raw read counts: sample_raw_reads.txt <br>
Covariate matrix: cov_matrix.txt <br>

For read counts file, the first three column is gene id, the second and third column is allele types (Ignore this here, type all A/T/C/G in one column). 

For the covariate matrix, intercept is manually required.

Use DESeq2 to normalize the raw read counts, for details, see normalization.R

The file generated from normalization.R is :  <br>
normalized_reads.txt <br>
cov_bim.txt <br>

Use gzip command to get the compressed read counts .gz file which is required for GEMMA.

## 2.LMM test using GEMMA 

Required file: <br>
gzip read counts: normalized_reads.txt.gz <br>
Phenotype: phenotype.txt <br>
Covariate matrix: cov_bim.txt <br>

Get Kinship matrix <br>
```
gemma -g normalized_reads.txt.gz -p phenotype.txt -c cov_bim.txt -gk 2 -notsnp -o cov_mat
```
The default output file should in the output folder under the data directory, cov_mat.sXX.txt. Notice that this is generated by using gene expression data, not the really kinship matrix.


LMM <br>
```
gemma -g normalized_reads.txt.gz -p phenotype.txt -k output/cov_mat.sXX.txt -c cov_bim.txt -lmm 4 -notsnp -o output. 
```
This will conduct the DGE anlysis by LMM GEMMA approach. The default output file should in the output folder under the data directory, output.assoc.txt. <br>

p_wald, p_lrt and p_score is the p-value for wald test, likelihood ratio test and score test for the differentially expressed.

## 3.Visualization 
Required file <br>
sample data: sample_data.txt <br>
Create qq plot, manhattan plot and volcano_plot.

