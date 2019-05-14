# Filtering PrediXcan Association Results
After we predict expression and associate it with out phenotype, we need to filter by p-value.  To do this, we need to read in all of the association files.  From here, we need to look at our total number of gene-tissue tests and divide by the number of tissues.  This will give us an average number of tests per tissue, so we do not favor one tissue over another.

Here is an example for the tissues from GTEx Version 7:

```{r}
GTEX7<-fread("z://PrediXcan/tissuelistV7",header=FALSE)
GTEX7<-GTEX7$V1
genenames<-fread("Z://PrediXcan/BP_Chrome.txt",header =T)# A list of genes that include their ENSG ID, their gene name, the chromosome, and position.

for(i in GTEX7){
		a<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_" %&% tis %&% "_association.txt",header = T,stringsAsFactors = FALSE)
  a$tissue<-tis
  tisswnames7<-left_join(a,genenames,by=c("gene"))
  if(exists("alltiss7")){
    alltiss7<-bind_rows(alltiss7,tisswnames7)
  }else{
    alltiss7<-tisswnames7
  }
}

row.has.na7 <- apply(alltiss7, 1, function(x){any(is.na(x))})
sum(row.has.na7)
newtiss7<-alltiss7[!row.has.na7,] #quickly remove NAs

significant7<-subset(newtiss7, p<= (0.05/5179))#5179 is the average number of tests per tissue.  We got this number by dividing 248592 (The total number of gene-tissue tests) by 48 (the number of tissues in GTEx V7))
#No significant genes predicted with GTEx V7
fwrite(alltiss7, "/home/peter/AA_nonGAIN_SCZ/PrediXcan/GTExV7AllResults.txt", col.names = T, row.names = F, quote=F, sep="\t")
```
The method above is an easy way to filter PrediXcan results for any set of models.  This is the same method we used for filtering our MESA and DLPFC results, but there were not significant genes,

GTEx V6 provides and interesting situation because not all genes were reliably predicted.  (See [Im Lab blog](http://hakyimlab.org/post/2017/v7-v6p-analysis/)).  As a result, we can only use genes that were flagged as green. 
```
database_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
for (i in database_tissues){
  predixcanresult<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/1000GTW_"%&%i %&% "_0.5.db_association.txt", header =T)
  predixcanresult$tissue <- i
  if(exists("total")){
    total<-bind_rows(total,predixcanresult)
  }else(total<-predixcanresult)
}

flags<-fread("Z://PrediXcan/flags.txt", header=T)
gtex6Wflags<-left_join(total,flags, by= c("gene","tissue"="model"))

row.has.naflags <- apply(gtex6Wflags, 1, function(x){any(is.na(x))})
sum(row.has.naflags)
newgtex6vflags<-gtex6Wflags[!row.has.naflags,]

noredflags<-newgtex6vflags[!(newgtex6vflags$flag=="red"),]
greenflags<-noredflags[!(noredflags$flag=="yellow"),]#THE NUMBER OF LINES IN THIS FILE ARE THE TOTAL NUMBER OF GENES WE WILL BE USING
#greenflags are the viable association results.
significant6<-subset(total, p<(0.05/(119507/44))
```
