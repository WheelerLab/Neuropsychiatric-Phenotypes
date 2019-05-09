# GWAS Quality Control for MGS GAIN data (In progress)

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



