#March 30, 2019
#Notes for the generation of a heat map regarding gene2pheno replications

library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)

gene2pheno<-read.csv("C://Users/Peter Fiorica/Documents/Wheeler Lab/Paper Figures/G2Pheatmap.csv", header =T)
gene2pheno$tissue<-gsub("TW_","",gene2pheno$tissue)
gene2pheno$tissue<-gsub("_Elastic_Net_0.5","",gene2pheno$tissue)


gene2pheno<- gene2pheno[!gene2pheno$tissue=="DGN_WB",]
  
G2Phits<-ggplot(gene2pheno)+
  aes(x = tissue, y = gene_name, fill = -log10(pval))+
  geom_raster() +
  scale_fill_gradient2(high="royalblue3", mid="aliceblue",name = "-log(P)" ) + 
  theme_bw() +
  ggtitle("B")+
  theme(panel.grid.major.y =element_blank(),
        axis.title= element_text(size=14,face="bold"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=0.05),
        plot.title=element_text(size = 18, face="bold")
        )+   facet_grid(. ~ "Replication: PGC S-PrediXcan") + 
  theme(strip.text.x = element_text(size = 15))+
  ylab("Gene") +
  xlab("Tissue")

GTex6SCZ<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/greenflagsnames.txt" , header =T)
GTex6BD<-fread("z://bipolar_disorder/predixcan/BPgtex6greenflags.txt", header =T)

SCZgenes<-subset(GTex6SCZ, GTex6SCZ$genename== "ZBTB1" | GTex6SCZ$genename== "SMPD3"| GTex6SCZ$genename== "SH3PXD2B"|  GTex6SCZ$genename== "PRMT7"| GTex6SCZ$genename=="FLVCR1")
BDgene<- subset(GTex6BD, GTex6BD$genename== "ZNF562")

siggenes1<-rbind(SCZgenes,BDgene, fill = T)

gene<-c("ENSG00000174705.7", "ENSG00000171466.5", "ENSG00000103056.7", "ENSG00000132600.12")
beta<-c(0,0,0,0)
t<-c(0,0,0,0)
p<-c(NA,NA,NA,NA)
tissue<-c("Ovary", "Testis", "Uterus","Vagina")
genename<-c("SH3PXD2B", "ZNF562", "SMPD3", "PRMT7")
flag<-c("green", "green", "green", "green")
CHR<-c(5,NA,16,16)
BP<- c(171752185, NA, 68392457, 68344877)
gene_name<- c("ZBTB1", NA, "SMPD3", "PRMT7")
for_append<-data.frame(gene,beta,t,p,tissue, genename, flag, CHR, BP, gene_name)
siggenes<-rbind(siggenes1,for_append, fill = T)

fwrite(siggenes, "C://Users/Peter Fiorica/Documents/Wheeler Lab/siggenes4heatmap.txt", col.names = T, row.names = F, sep = "\t", quote = F)
siggenes<-fread("C://Users/Peter Fiorica/Documents/Wheeler Lab/siggenes4heatmap.txt", header =T)

siggeneshits<-ggplot(siggenes)+
  aes(x = tissue, y = genename, fill = -log10(p))+
  geom_raster() +
  scale_fill_gradient2(na.value = NA, high="royalblue3", mid="aliceblue",name = "-log(P)", limits=c(0,8.5), breaks=c(2,4,6,8)) + 
  theme_bw() +
  ggtitle("A")+
  theme(panel.grid.major.y =element_blank(),
        axis.title= element_text(size=14,face="bold"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=.05),
        plot.title=element_text(size = 14, face="bold")
  )+   facet_grid(. ~ "Discovery: GAIN PrediXcan") + 
  theme(strip.text.x = element_text(size = 15))+
  ylab("Gene") +
  xlab("Tissue")

#Side by side figure
grid.arrange(grobs= list(G2Phits + theme(legend.position="none"),siggeneshits + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())),
             nrow=1, gp=gpar(fontsize=20, fontface="bold"))
#Top-Bottom Figure
grid.arrange(grobs= list(siggeneshits, G2Phits),
             nrow=2)
