# GWAS Quality Control for MGS GAIN data

### Step 0: Merging the nonGAIN dataset with the GAIN dataset, setting heterozygous haploid genotypes as missing, performing sex check
```
plink --bfile /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep0/premerge_wGAIN --bmerge /home/peter/AA_GAIN_SCZ/AA_bfiles/Merged_AA_GAIN.bed /home/peter/AA_GAIN_SCZ/AA_bfiles/Merged_AA_GAIN.bim /home/peter/AA_GAIN_SCZ/AA_bfiles/Merged_AA_GAIN.fam --make-bed --out QCStep0

plink --bfile QCStep0 --set-hh-missing --make-bed --out NoHH/QCStep0NoHH

plink --bfile NoHH/QCStep0NoHH --check-sex --missing --out QCStep0SexCheck
```

### Step 1: Identifying Unfiltered Genotyping Rate 
```
plink --bfile ../QCStep0/NoHH/QCStep0NoHH  --missing --out QCStep1
```
#### Step 1A: Plotting unfiltered genotyping rate
These steps are completed in R
```
"%&%"=function(a,b) paste(a,b,sep="")
my.dir<-"Z://AA_nonGAIN_SCZ/QCSteps/"
lmiss<-fread(my.dir%&%"QCStep1/QCStep1.lmiss",header = T)
hist(lmiss$F_MISS) #This creates a histogram of the missingness of the data before we filter by genotyping rate.
dim(lmiss)[1] #This tells us the number of SNPs we are working with before filtering by genotyping rate
table(lmiss$F_MISS<0.01)
table(lmiss$F_MISS<0.02)
sum(lmiss$F_MISS<0.01)/(dim(lmiss)[1])
sum(lmiss$F_MISS<0.02)/(dim(lmiss)[1])#The percent of SNPs have a genotyping call rate of 98%
```

### Step 2: Filtering SNPs by Genotyping Rate
```
plink --bfile /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep0/NoHH/QCStep0NoHH --geno 0.01 --make-bed --out /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep2/QCStep2
```
 
### Step 3: Identifying Filtered Genotype Rate
```
plink --bfile ../QCStep2/QCStep2 --missing --out QCStep3
```

#### Step 3A: Plotting Filtered Genotyping Rate
These steps are completed in R
```
newimiss<-fread(my.dir%&%"QCStep3/QCStep3.imiss")
hist(newimiss$F_MISS)
newlmiss<-fread(my.dir%&%"QCStep3/QCStep3.lmiss")
hist(newlmiss$F_MISS)
dim(newlmiss)[1]
```

### Step 4: Filtering by HWE
```
plink --bfile ../QCStep2/QCStep2 --hardy --out QCStep4
```

#### Step 4A: Plotting HWE Frequencies and Removing SNPs outside of HWE
These steps are completed in R
```
hwe<-fread(my.dir%&%"QCStep4/QCStep4.hwe",header =T)
summary(hwe$P)
hist(hwe$P)
table(hwe$P<1e-06)
table(hwe$P<1e-06)/sum(table(hwe$P<1e-06))
```

### Step 5: IBD Pruning
```
plink --bfile /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep2/QCStep2 --indep-pairwise 50 5 0.3 --out QCStep5a
#This step removed over 400,000 SNPs
```

#### Step 5A: Extracting individuals with excess IBD
```
plink --bfile /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep2/QCStep2 --extract ../QCStep5A/QCStep5a.prune.in --genome --out QCStep5B
```
The next steps are completed in R
```
ibd<-fread(my.dir %&% "QCStep5/QCStep5B/QCStep5B.genome",header=T)
ggplot(data = ibd, aes(x=Z0,y=Z1))+geom_point(alpha=1/4)+theme_bw()
##Now we can check for duplicates in the data
dups<-data.frame()
for( i in 1:dim(ibd)[1]){
  if(as.character(ibd$IID1[i])==as.character(ibd$IID2[i])){
    dups<-rbind(dups,ibd[i,])
  }
}
dim(dups)
#It looks like there are some individuals clustered at or around (0,0), suggesting that these individuals may be identical twins.  These individuals will need to be removed.
hapmap <- filter(ibd,grepl('NA',IID1))
#No hapmap individuals. No surpise.
toExclude <- c(as.character(dups$IID1),as.character(hapmap$IID1))
a <- as.character(ibd$IID1) %in% toExclude
others <- ibd[a==FALSE,]
#Isolating individuals that need to be removed.
toremove<-filter(others,PI_HAT>=0.2)
write.table(toremove,my.dir%&%"QCStep5/QCStep5B/Relate.to.remove.txt",quote=FALSE, row.names = FALSE)
```
#### Step 5B: Extracting individuals with excess heterozygosity
```
plink --bfile ../../QCStep2/QCStep2 --extract ../QCStep5A/QCStep5a.prune.in --het --out QCStep5c
```
The next steps are completed in R
```
HET<-fread(my.dir%&%"QCStep5/QCStep5C/QCStep5c.het",header =T)
h=HET$"N(NM)"-HET$"O(HOM)"/HET$"N(NM)"
oldpar = par(mfrow=c(1,2))
hist(h,50)
hist(HET$F,50)
summary(HET$F)
abline(v=mean(HET$F)+6*sd(HET$F),col="red")
abline(v=mean(HET$F)-6*sd(HET$F),col="red")

sortHET <- HET[order(HET$F),]
outliers <- data.table()

for(i in 1:length(sortHET$F)){
  if(sortHET[i,6] > (mean(sortHET$F)+3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
  if(sortHET[i,6] < (mean(sortHET$F)-3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
}

hetoutliers <- select(outliers, FID, IID)
dim(hetoutliers) #This tells us how many outliers there are.

write.table(allexclude2, file = my.dir%&%"QCStep5/QCStep5C/HetOutliers.txt", quote = F, col.names = T, row.names = F)
```

#### Step 5C: Removing individuals with excess heterozygosity or relatedness
```
plink --bfile ../../QCStep2/QCStep2 --extract ../QCStep5A/QCStep5a.prune.in --remove ../QCStep5B/Relate.to.remove.txt --genome --out QCStep5D
plink --bfile ../../QCStep2/QCStep2 --extract ../QCStep5A/QCStep5a.prune.in --remove ../QCStep5B/Relate.to.remove.txt --make-bed --out QCStep5D
```
The next steps are completed in R
```
IBD<-fread(my.dir %&% "QCStep5/QCStep5D/QCStep5D.genome",header=T)
ggplot(data = IBD, aes(x=Z0,y=Z1))+geom_point(alpha=1/4)+theme_bw()
```

#### Step 5D: Second Heterozygosity Check
```
plink --bfile ../QCStep5D/QCStep5D  --het --extract ../QCStep5A/QCStep5a.prune.in --remove ../QCStep5B/Relate.to.remove.txt --out QCStep5E
```
The next steps are completed in R...We really just repeat the steps from 5B
```
HET<-fread(my.dir%&%"QCStep5/QCStep5E/QCStep5E.het",header =T)
h=HET$"N(NM)"-HET$"O(HOM)"/HET$"N(NM)"
oldpar = par(mfrow=c(1,2))
hist(h,50)
hist(HET$F,50)
summary(HET$F)
abline(v=mean(HET$F)+6*sd(HET$F),col="red")
abline(v=mean(HET$F)-6*sd(HET$F),col="red")

sortHET <- HET[order(HET$F),]
outliers <- data.table()

for(i in 1:length(sortHET$F)){
  if(sortHET[i,6] > (mean(sortHET$F)+3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
  if(sortHET[i,6] < (mean(sortHET$F)-3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
}

hetoutliers <- select(outliers, FID, IID)
dim(hetoutliers)
```
#### Step 5E: Remove Heterozygosity Outliers
```
plink --bfile ../QCStep5D/QCStep5D --remove ../QCStep5C/HetOutliers.txt --make-bed --out QCStep5F
```

## Principal Component Analysis
### Step 6: Merge with HapMap Data
Dependent on the genotyping platoform used (i.e. Affymatrix, Illumina, etc.), the SNP identifiers are recorded differently (rsID, SNP_A-#, AFFY-SNP-#). The hapmap individuals are recorded by rsID, and they share different positions from SNP_A-#.  This means that we have to change the MGS data identifiers to rsID.

The next steps are completed in R
```
TotalSNPs<-fread("Z://AA_GAIN_SCZ/summarystatistics.txt",header=T)
#The summary statistics from the dbGaP data containing information that contains an rsID and position for each SNP_A-#
SelectSNPs<-dplyr::select(TotalSNPs,"MarkerAccession","ChrID", "ChrPosition","SubmittedSNPID")
#Isolates the positions we need to merge with a .bim file
bim<-fread("Z://AA_nonGAIN_SCZ/QCSteps/QCStep5/QCStep5F/QCStep5F.bim", header = F)
#Reading in the last .bim from QC
mergedbim <- left_join(bim, SelectSNPs, by = c("V1" = "ChrID", "V4" = "ChrPosition"))
mergedbim2 <- mutate(mergedbim,snp=ifelse(is.na(`Marker accession`),V2, `Marker accession`))
newbim <- dplyr::select(mergedbim2,V1,snp,V3,V4,V5,V6)
filetest2<-newbim[!duplicated.data.frame(newbim),]

write.table(filetest2,"Z:/AA_nonGAIN_SCZ/QCSteps/QCStep5/QCStep5F/newbim.bim",quote=F, sep="\t",row.names=F,col.names=F)
```
From here, we return to commandline to synchronize and merge with plink
#### Step 6A: Synchronize PLINK bfiles
```
plink --fam QCStep5F.fam --bed QCStep5F.bed --bim newbim.bim --make-bed --out newbfiles
```

#### Step 6B: Merge GWAS data with HapMap bfiles
Note that the bfiles are stil in Hg18 build.
Additionally, this attempt will fail, but it will output a .missnp file.
```
plink --bfile ../newbfiles --bmerge ../HAPMAP3_hg18/HM3_ASN_CEU_YRI_Unrelated_hg18_noAmbig.bed ../HAPMAP3_hg18/HM3_ASN_CEU_YRI_Unrelated_hg18_noAmbig.bim ../HAPMAP3_hg18/HM3_ASN_CEU_YRI_Unrelated_hg18_noAmbig.fam --make-bed --out QCStep6A
```

#### Step 6C: Exclude missSNPs
```
plink  --bfile ../HAPMAP3_hg18/HM3_ASN_CEU_YRI_Unrelated_hg18_noAmbig --exclude ../QCStep6A/QCStep6A-merge.missnp --make-bed --out QCStep6B
```

#### Step 6D: Merge attempt 2
```
plink --bfile ../newbfiles --bmerge ../QCStep6B/QCStep6B.bed ../QCStep6B/QCStep6B.bim ../QCStep6B/QCStep6B.fam --make-bed --out QCStep6C
```

#### Step 6E: Filter by Minor Allele Frequency and Genotyping Rate
```
plink --bfile ../QCStep6C/QCStep6C --geno 0.01 --maf 0.05 --chr 1-22 --make-bed --out QCStep6D
```

#### Step 6F: IBD Pruning with HapMap data
```
plink --bfile ../QCStep6D/QCStep6D --indep-pairwise 50 5 0.3 --recode --out QCStep6E
awk '{print $1,$2,$3,$4,$5,1}' ../QCStep6D/QCStep6D.fam > QCStep6E.fam
```

#### Step 6G: Running principal component analysis with smartpca or PLINK
```
perl ../make_par_file.pl ../QCStep6E/QCStep6E 0 > QCStep6F.par
smartpca -p QCStep6F.par
#I also want to see how plinks PCA compares to EIGENSOFT pca,
#So I also ran the next command:
plink --bfile QCStep6D --pca --out QCStep6D_PCA
#PLINK will output eigenvectors and eigenvalues that can be used for the next steps, but I chose to continue using the smartpca ourput
```

### Step 7: Plotting PCA data
The next steps are in R
```
pcsgenonox0.01<-read.table(my.dir%&%"QCStep6/QCStep6E/QCStep6E.evec",skip=1)
pcdfgenonox0.01 <- data.frame(popinfo, pcsgenonox0.01[,2:11]) %>% rename (PC1=V2,PC2=V3,PC3=V4,PC4=V5,PC5=V6,PC6=V7,PC7=V8,PC8=V9,PC9=V10,PC10=V11)
gwasgenonox0.01 <- filter(pcdfgenonox0.01,pop=='GWAS')
hm3genonox0.01 <- filter(pcdfgenonox0.01, grepl('NA',IID))
evalgenonox0.01<-scan(my.dir%&%"QCStep6/QCStep6E/QCStep6E.eval")[1:10]
round(eval/sum(evalgenonox0.01))#Calculate the percent explained by each PC

ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() + scale_colour_brewer(palette="Set1")
ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC1,y=PC3,col=pop,shape=pop))+geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() + scale_colour_brewer(palette="Set1")
ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC2,y=PC3,col=pop,shape=pop))+geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() + scale_colour_brewer(palette="Set1")
```

The groups should cluster into four distinct populations for the plot of PC1 vs. PC2  Three small populations should be easily identifiable: YRI, CEU, and ASN.  The GWAS population should represent a band between YRI and CEU since African American individuals are genetically admixed.
