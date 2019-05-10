# Running PrediXcan with a Residual Phenotype
Original PrediXcan does not have a flag for running a logistic regression with covariates.  We have to use the covariates to account for population structure.  We also need to use a logistic regression because our phenotype is case-control.   Below is a short set of commands to perform a logistic regression of our phenotype.

```
library(dplyr)
library(data.table)
eigenvec<-fread("Z://AA_nonGAIN_SCZ/GWAS/PCAforCovariates.evec", header=F)
pheno<-fread("Z://AA_nonGAIN_SCZ/QCSteps/QCStep2/QCStep2.fam",header=F)
fit<-glm(pheno$V6~eigenvec$V3 + eigenvec$V4 + eigenvec$V5 + eigenvec$V6 + eigenvec$V7 + eigenvec$V8 + eigenvec$V9 + eigenvec$V10 + eigenvec$V11 + eigenvec$V12)

nophenofam<-fread("Z://AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/UMichFiltered/1000GFilteredPlink.fam", header =F)
nophenofam$V5<-pheno$V5
nophenofam$V6<-pheno$V6
fwrite(nophenofam, "Z://AA_nonGAIN_SCZ/GWAS/Plinkbfiles/1000Gfilteredrs.fam", col.names = F, row.names = F, sep = " ", quote =F)

newphenotype<-fit$residuals
residualphenos<-data.table(x=pheno$V1,y=pheno$V2,z=newphenotype)
colnames(residualphenos)<-c("FID","IID","Pheno")
write.table(residualphenos,file="Z://AA_nonGAIN_SCZ/GWAS/residualphenoMGS.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
```
Now that we have our residual phenotype, we can perform PrediXcan.
Below is the short script I made to run it.
```
for tiss in `cat tissuelist3`; do
    echo $tiss
    PrediXcan.py --predict --assoc --linear --dosages /home/peter/AA_nonGAIN_SCZ/PrediXcan/dosages/CAAPA/ --dosages_prefix CAAPApredixcandosage --samples /home/peter/AA_nonGAIN_SCZ/PrediXcan/dosages/samples.txt --weights /home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/${tiss} --pheno /home/peter/AA_nonGAIN_SCZ/GWAS/residualphenoMGS.txt --output_prefix /home/peter/AA_nonGAIN_SCZ/PrediXcan/output/CAAPA/CAAPA_${tiss}
done
#tissuelist 3 is a list of the 44 GTEx version 6 tissues.
#I repeated this command with a list of GTEx version 7, MESA predictors, and DLFPC
```
After PrediXcan runs, we have predicted expression and association files.  We can filter our association files to look for significantly associated genes.
