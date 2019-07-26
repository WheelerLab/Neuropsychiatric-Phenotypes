#Plotting PrediXcan SNPs MAFs for 1000G and CAAPA
#Peter Fiorica
#5 July 2019
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
"%&%" = function(a,b) paste(a,b,sep="")


#GTEx6SNPs<-fread("/home/peter/AA_nonGAIN_SCZ/PrediXcan/pullingSNPsFromPredictors/RSIDSGTEXV6_NODUPLICATES.TXT", header =F)
GTEx7SNPs<-fread("/home/peter/wl3backup/RSIDSGTEXV7_NODUPLICATES.TXT", header =F)
GTEx7SNPs<-data.frame(GTEx7SNPs)
colnames(GTEx7SNPs)<-c("SNP")
GWAS_1000G<- fread("/home/peter/wl3backup/1000G_gwas_maf.frq", header =T)
GWAS_CAAPA<- fread("/home/peter/CAAPAnonGAINSCZ.assoc.dosage", header =T)


for (i in 1:22){
  a<-fread("/home/peter/wl3backup/imputation_summary/1000G/chr"%&%i%&%".info", header = T)
  if(exists("SNPs1000G")){
    SNPs1000G<-bind_rows(SNPs1000G,a)
  }else{
    SNPs1000G <- a
  }
}

for (j in 1:22){
  b <- fread ("/home/peter/wl3backup/imputation_summary/CAAPA/chr"%&%j%&%".info", header =T)
  if(exists("CAAPA_SNPs")){
   CAAPA_SNPs<-bind_rows(CAAPA_SNPs,b)
  }else{
    CAAPA_SNPs<-b
  }
}

bothrefSNPs<-left_join(SNPs1000G,CAAPA_SNPs, by = "SNP")

unfiltered_A<-ggplot(data=bothrefSNPs, aes(x=(MAF.x),y=MAF.y))+
  geom_point(alpha=0.05)+
  geom_abline(col="royalblue")+
  ylab("CAAPA")+xlab("1000G")+
  facet_wrap(~as.character("Unfiltered"), nrow =1) + 
  ggtitle("A")+
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=15)
        )+
  theme_bw(20)

bothSNPs<-left_join(GWAS_CAAPA,GWAS_1000G, by ="SNP")

both_px<-left_join(GTEx7SNPs,bothSNPs, by = "SNP")

px_plot<- ggplot(data=both_px, aes(x=(1-FRQ),y=MAF))+
  geom_point(alpha=0.05)+
  geom_abline(col="royalblue")+
  xlab("CAAPA")+ylab("1000G")+
  facet_wrap(~as.character("GTEx PrediXcan SNPs"), nrow = 1) + 
  ggtitle("B")+
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=15)
  )+
  theme_bw(20)

gwas_plot<- ggplot(data=bothSNPs, aes(x=(1-FRQ),y=MAF))+
  geom_point(alpha=0.05)+
  geom_abline(col = "royalblue")+
  xlab("CAAPA")+ylab("1000G")+ ggtitle("C")+
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=15)
  ) +
  facet_wrap(~as.character("GWAS SNPs"),nrow=1) + 
  theme_bw(20)


tiff("/home/peter/mafs4allsnps1.tiff", height = 700, width = 2100, compression = "lzw")
grid.arrange(grobs= list(unfiltered_A, px_plot,gwas_plot), ncol=3, nrow =1)
dev.off()

png("/home/peter/mafs4allsnps1.png", height = 700, width = 2100)
grid.arrange(grobs= list(unfiltered_A, px_plot,gwas_plot), ncol=3, nrow =1)
dev.off()

