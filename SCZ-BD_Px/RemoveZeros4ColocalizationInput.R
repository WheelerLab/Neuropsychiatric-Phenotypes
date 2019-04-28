library(dplyr)
library(data.table)
library(R.utils)
"%&%" = function(a,b) paste(a,b,sep="")

database_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

for (i in database_tissues){
  a<-fread(sprintf("zcat %s","/home/angela/Ad_PX_pipe_TEST_POPS/peter_SCZ/COLOC_input/SCZ_GWAS_" %&% i %&% ".txt.gz"), header = T)
  a[a==0]<-NA
  row.has.na <- apply(a, 1, function(x){any(is.na(x))})
  b<-a[!row.has.na,]
  fwrite(b,"/home/peter/AA_nonGAIN_SCZ/PrediXcan/Colocalization/COLOC_input/SCZ_GWAS_" %&% i %&% ".txt", col.names = T, row.names = F, sep="\t", quote =F)
  gzip("/home/peter/AA_nonGAIN_SCZ/PrediXcan/Colocalization/COLOC_input/SCZ_GWAS_" %&% i %&% ".txt", overwrite=T)
  print("Completed " %&% i)
  }
