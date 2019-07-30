# ashr application to PrediXcan results

```
library(data.table)
library(dplyr)
library(ashr)

"%&%" = function(a,b) paste(a,b,sep="")
brain_tissues<-c('Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex','Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia')
```

`ashr` is an R package built for adaptive shrinkage to calculate local false sign rate (lfsr).  lfsr is an Empirical Bayes approach for large-scale hypothesis testing.  It follows the assumption that the distribution of the effects is unimodal at zero.  This method uses effect size and the standard error to calculate a probability of getting a sign of an effect wrong.
See [Matthews Biostatistics 2017](https://academic.oup.com/biostatistics/article/18/2/275/2557030)

```
GTEv7.1<-fread("Z://PrediXcan/tissuelistV7",header =F)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub(".{35}$","",GTEv7.1$V1)
GTEv7<-GTEv7.1$V1

for (i in GTEv7){
  one<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_"%&% i %&%"_imputed_europeans_tw_0.5_signif.db_association.txt", header = T)
  one$tissue<-i
  if(exists("predixcan_results1")){
    predixcan_results1<-bind_rows(predixcan_results1,one)
  }else{
    predixcan_results1<-one
  }
}

ash_pX_half_uni1 <- ash(predixcan_results1$beta, predixcan_results1$`se(beta)`, mixcompdist = 'halfuniform', method='fdr')

scz <- mutate(predixcan_results1, ash_halfuni_pX=ash_pX_half_uni1$result$lfdr)

scz<-scz %>%
  arrange(ash_halfuni_pX)


fwrite(scz,"z://AA_nonGAIN_SCZ/PeerJReviews/gtex7_all_ashr_results.txt", col.names=T, row.names=F, quote =F, sep = "\t")
```
