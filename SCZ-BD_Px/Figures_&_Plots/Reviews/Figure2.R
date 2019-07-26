#PRMT7 Violin Plot for Reviews
#Peter Fiorica
#19 July 2019 

library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
"%&%" = function(a,b) paste(a,b,sep="")


pheno<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/phenobinary1.txt", header =F)
pheno$V3[pheno$V3=="1"]<-"SCZ"
pheno$V3[pheno$V3=="0"]<-"Control"
colnames(pheno)<-c("V1","V2","Phenotype")


pitXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Heart_Atrial_Appendage_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
pit_gene<-data.table(FID=pitXpress$FID,IID=pitXpress$IID,Xpres=pitXpress$ENSG00000132600.12)
With_phenopit<-left_join(pit_gene , pheno , b = c("IID" = "V2"))
onlywithphenopit<-dplyr::select(With_phenopit,FID,IID,Xpres,Phenotype)


PRMT7title<-expression(paste(italic("PRMT7"), " (Heart Atrial Appendage)"))
PRMT7 <-ggplot(data = onlywithphenopit, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(PRMT7title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()


tiff("Z://AA_nonGAIN_SCZ/PeerJReviews/plots/violin_plot_prmt7_figure2.tiff", width=837, height=678)
PRMT7
dev.off()
