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
```
#These steps are completed in R
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
#These steps are completed in R
newimiss<-fread(my.dir%&%"QCStep3/QCStep3.imiss")
hist(newimiss$F_MISS)
newlmiss<-fread(my.dir%&%"QCStep3/QCStep3.lmiss")
hist(newlmiss$F_MISS)
dim(newlmiss)[1]
```

### Step 4: Filtering by 
```
plink --bfile ../QCStep2/QCStep2 --hardy --out QCStep4
```


