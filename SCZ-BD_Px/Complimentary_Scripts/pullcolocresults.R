#Script to Pull COLOC P1-P4 values for PRMT7 and SMPD3

library(dplyr)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

database_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

for (i in database_tissues){
  a<-fread(sprintf("zcat %s","/home/peter/AA_nonGAIN_SCZ/PrediXcan/Colocalization/COLOC_results/SCZ_" %&% i %&% ".txt.gz"), header = T)
  a$tissue<-i
  c<-subset(a,gene_id == "ENSG00000132600.12" | gene_id == "ENSG00000103056.7")
  if(exists("b")){
    b<-rbind(b,c)
  }else{
    b<-c
  }
  print("Completed " %&% i)
}

fwrite(b,"/home/peter/AA_nonGAIN_SCZ/PrediXcan/Colocalization/SCZ_GWAS_PRMT7_SMPD3.txt", col.names = T, row.names = F, sep="\t", quote =F)
