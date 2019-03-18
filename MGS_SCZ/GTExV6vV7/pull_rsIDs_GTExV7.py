#pull rsids from GTEx V7 db files

import sqlite3
import pandas
pops = ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary","Artery_Tibial","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Transformed_fibroblasts","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Muscle_Skeletal","Minor_Salivary_Gland","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"]

RSID = []
for pop in pops:
    conn = sqlite3.connect("/home/wheelerlab3/Data/PrediXcan_db/GTEx-V7_HapMap-2017-11-29/gtex_v7_" + pop + "_imputed_europeans_tw_0.5_signif.db") #connect to .db file
    c = conn.cursor()
    c.execute('SELECT (`rsid`) FROM WEIGHTS;') #run search through python
    for row in c:
      row = row[0].split(" ") #tuple is just one element, so split into two
      RSID.append([row[0], pop]) #extracts rsIDs
    print(pop)
    conn.close()
genename_R2_df = pandas.DataFrame(RSID, columns = ['rsid', 'tissue']) #make list of lists into df
genename_R2_df.to_csv("/home/peter/AA_nonGAIN_SCZ/PrediXcan/pullingSNPsFromPredictors/rsidsGTExV7.txt", index = False) #print to file
