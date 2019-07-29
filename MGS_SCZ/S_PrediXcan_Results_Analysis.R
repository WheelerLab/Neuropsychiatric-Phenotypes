#Reading in S-PrediXcan data for PGC_SCZ Results
#Peter Fiorica
#25 June 2019

library(dplyr)
library(data.table)

database_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
"%&%" = function(a,b) paste(a,b,sep="")

for (tis in database_tissues){
  a<-fread("Z://AA_nonGAIN_SCZ/PGC_SCZ_Summary_Stats/SPrediXcan_Results/TW_"%&% tis %&%"_S-PrediXcan_Results_SCZ.txt", header =T)
  a$tissue<-tis
  if(exists("MetaXcan_Results")){
    MetaXcan_Results<-bind_rows(MetaXcan_Results,a)
  }else(MetaXcan_Results<-a)
}

#This works for every tissue except Esophagus muscularis.
#The script yielded a segmentation when I ran it for chr2ls 
#In this tissue, the MetaXcan script will yield a segmentation fault on the eleventh chromosome read.
#Note that this is not necessarily Chr 11, but whatever chromosome is read in as the 11th.
#esophagus_spredixcan.sh: line 1: 26070 Segmentation fault

#What happens if I only read in 10 chromosomes?

#Read in chr 10-19
#INFO - Processing input gwas
#INFO - Successfully parsed input gwas in 14.9019060135 seconds
#INFO - Started metaxcan process
#INFO - Loading model from: /home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/TW_Esophagus_Muscularis_0.5.db
#INFO - Loading covariance data from: /home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/TW_Esophagus_Muscularis.txt.gz
#INFO - Processing input gwas
#INFO - Started metaxcan association
#INFO - 10 % of model's snps found so far in the gwas study
#INFO - 20 % of model's snps found so far in the gwas study
#INFO - 30 % of model's snps found so far in the gwas study
#INFO - 40 % of model's snps found so far in the gwas study
#INFO - 50 % of model's snps found so far in the gwas study
#INFO - 60 % of model's snps found so far in the gwas study
#INFO - 70 % of model's snps found so far in the gwas study
#INFO - 78 % of model's snps used
#INFO - Sucessfully processed metaxcan association in 31.2243731022 seconds

#So I tried doing this for other chromosomes one by one, but chromosome 2 is the one giving me the issue. . .
#I'm gonna try splitting it in half. . . 

#gunzip -cpgc_scz_gwas_chr2.txt.gz  | split -l 200000 - pgc_scz_gwas_chr2.txt.gz.part

for (tis in database_tissues){
  a<-fread("Z://bipolar_disorder/PGC_BIP_SummaryStats/output/TW_"%&% tis %&%"_S-PrediXcan_Results_BIP.txt", header =T)
  a$tissue<-tis
  if(exists("MetaXcan_Results")){
    MetaXcan_Results<-bind_rows(MetaXcan_Results,a)
  }else(MetaXcan_Results<-a)
}
