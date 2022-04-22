#################################################
# load data
sample_counts<-fread("sample_raw_reads.txt")
sample_raw_read<-sample_counts[,4:100]
cov_matrix=read.table("cov_matrix.txt",header=TRUE)

#################################################
# Delete samples with low map read 
mapped_read=apply(sample_raw_read,2,sum)
summary(mapped_read)

# from https://support.bioconductor.org/p/91218/
# TPM, select genes with TPM > 0.1 in at least 20% samples
sample_raw_read_length = sample_raw_read/sample_counts$length
tpm = t( t(sample_raw_read_length) * 1e6 / colSums(sample_raw_read_length) )
tpm_0.1 = apply(tpm>0.1,1,sum)/ncol(tpm)

#delete genes 
delete_tpm = which(tpm_0.1<0.2)
sample_counts_tpm = sample_counts[-delete_tpm,]

sample_raw_read_tpm=sample_counts_tpm[,4:100]

#################################################
library(DESeq2)
###DESeq2 normalization
###col_name
col_name=data.frame(condition=cov_matrix$condition)
col_name$condition=factor(col_name$condition)
###DDS data
dds <- DESeqDataSetFromMatrix(countData = sample_raw_read_tpm,
                              colData = col_name,
                              design = ~ condition)
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
normalized_Deseq2<-counts(dds, normalized=TRUE)
normalized_Deseq2<-round(normalized_Deseq2,2)
normalized_Deseq2_log<-round(log2(normalized_Deseq2+1),2)
#################################################
# prepare GEMMA sepcific format file, for detail please see GEMMA manual
normalized_bimbam_log2 <-data.frame(geneID=sample_counts_tpm$geneID,A="A",T="T",normalized_Deseq2_log)

write.table(normalized_bimbam_log2,"normalized_reads.txt",row.names=F, col.names = F, quote = F, sep = "\t")

# add intercept 
cov_bim = cbind(intercept=1,cov_matrix)

write.table(cov_bim,"cov_bim.txt",row.names=F, col.names = F, quote = F, sep = "\t")


