library("ggplot2")
library("tidyverse")
library("ggrepel")
###qqplot##################
qq_ggplot<-function(p_val){
  qq_dat <- data.frame(obs=-log10(sort(p_val)),
                       exp=-log10( ppoints(length(p_val))))
  pd_qq <- ggplot(data=qq_dat,aes(exp,obs))+
    geom_point(alpha=0.7)+
    geom_abline()+
    xlab("Expected -log10(P-value)")+
    ylab("Observed -log10(P-value)")+ 
    theme(axis.text.x = element_text(size = 14, face = "bold"),axis.text.y = element_text(size = 14, face = "bold")) +
    theme(axis.title.x = element_text(size = 14,face = "bold", vjust = 0.5))+ 
    theme(axis.title.y = element_text(size = 14,face = "bold", vjust = 0.5))
  return(pd_qq)
}


###Manhattan plot BP position##################
bp_manhattan<-function(data){
  don <- data %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(data, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, start) %>%
    mutate( BPcum=start+tot)
  
  axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  return(list(don,axisdf))
  
}

xlab_chr=c(1:22,"X","M")
###manhattan_ggplot2##################
manhattan_ggplot2<-function(data,axisdf,label){
  ymax=round(max(-log(data$p_wald,base=10)))+1
  p<-ggplot(data, aes(x=BPcum, y=-log(p_wald,base=10))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 24 )) +
    # custom X axis:
    scale_x_continuous( label = xlab_chr, breaks= axisdf$center )  +  
    scale_y_continuous(expand = c(0,0), limits=c(0, ymax),breaks=seq(0,ymax,2)) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())+
    theme(axis.text.x = element_text(size = 14, face = "bold"),axis.text.y = element_text(size = 12, face = "bold")) +
    theme(plot.title= element_text(size = 14, face = "bold"))+
    theme(axis.title.x = element_text(size = 14,face = "bold", vjust = 0.5))+ 
    theme(axis.title.y = element_text(size = 14,face = "bold", vjust = 0.5)) + xlab("Chromosome")+ylab("-log(p-value)") +
    #theme(axis.line = element_line(colour = "grey")) +
    geom_hline(aes(yintercept=4),color="blue") +
    #geom_hline(aes(yintercept=-log(2.5*10^-6,base=10)),color="red") + 
    geom_text_repel(data=label, aes(x=BPcum, y=-log(p_wald,base=10),label=gene_name), size=4,color="blue")
  
  return(p)
}


###Volcano plot##################

volcano_plot<-function(data,title){
  p=ggplot(data=data, aes(x=beta, y=-log10(p_wald), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(size = 5,fontface="bold") +
    scale_color_manual(values=c("black","red","blue")) + ylab("-log10(P-value)") +
    geom_vline(xintercept=c(quantile(data$beta, 0.05),quantile(data$beta, 0.95)), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") + 
    theme(axis.text.x = element_text(size = 12,face = "bold", vjust = 0.5)) + 
    theme(axis.text.y = element_text(size = 12,face = "bold", vjust = 0.5)) + 
    theme(axis.title.x = element_text(size = 12,face = "bold", vjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5)) + theme(legend.position = "none") +
    theme(plot.title= element_text(size = 12, face = "bold")) + ggtitle(title)
  return(p)
}

###Load sample data##################
#23=X 24=M
sample_data<-read.table("sample_data.txt",header=TRUE)

sample_data$chr<-as.numeric(sample_data$chr)

###qqplot example##################
qq_plot<-qq_ggplot(sample_data$p_wald)

###Manhattan plot example##################
bp=bp_manhattan(sample_data)
bp_data=bp[[1]]
bp_axisdf=bp[[2]]

label=bp_data[order(bp_data$p_wald,decreasing=FALSE)[1:5],]

p_man=manhattan_ggplot2(bp_data,bp_axisdf,label)

###Volcano plot example##################
# add a column of NAs
sample_data$diffexpressed <- "0"
# if log2Foldchange > 10 and pvalue < 0.05, set as "UP" 
sample_data$diffexpressed[sample_data$beta >quantile(sample_data$beta, 0.95)  & sample_data$p_wald < 0.05] <- "1"
# if log2Foldchange < -10 and pvalue < 0.05, set as "DOWN"
sample_data$diffexpressed[sample_data$beta < quantile(sample_data$beta, 0.05) & sample_data$p_wald < 0.05] <- "2"



#UP and DOWN top 5 significant gene list
UP_gene_list=sample_data[sample_data$diffexpressed=="1",]    
DOWN_gene_list=sample_data[sample_data$diffexpressed=="2",]    

up_top_gene=UP_gene_list[order(UP_gene_list$p_wald,decreasing=FALSE)[1:5],"gene_name"]
down_top_gene=DOWN_gene_list[order(DOWN_gene_list$p_wald,decreasing=FALSE)[1:5],"gene_name"]


up_down_list=c(up_top_gene,down_top_gene)

sample_data$delabel <- NA
sample_data$delabel[sample_data$gene_name%in%up_down_list] <- sample_data[sample_data$gene_name%in%up_down_list,"gene_name"]

p_cog=volcano_plot(sample_data,"sample data")


