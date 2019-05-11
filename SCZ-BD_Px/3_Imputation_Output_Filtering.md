# Filtering Imputation Output
Once imputation has been ran, we need to filter the SNPs for r<sup>2</sup> \> 0.8 and MAF \> .01.  Below are the generic PLINK commands than are ran to convert one chromosome file from the imputation server to PLINK bfiles.

```
plink --vcf /home/peter/AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/chr22.dose.vcf.gz --biallelic-only --make-bed --out 1000Gc22biallelic
plink --bfile 1000Gc22biallelic --write-snplist --out all22snps1000G
#This writes us a list of all the SNPs on the chromosome.
cat all22snps1000G.snplist | sort | uniq -d > duplicated22snps1000G.snplist
#This creates a list of duplicates based on the snplist
plink --bfile 1000Gc22biallelic --exclude  duplicated22snps1000G.snplist --make-bed --out 1000Gc22nodups
#this removes the duplicates on the snplist
plink --bfile 1000Gc22nodups --qual-scores /home/peter/AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/chr22.info 7 1 1 --qual-threshold 0.8 --maf 0.01 --make-bed  --out /home/peter/AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/UMichFiltered/chr22
#This  makes the final bfile set for the chromsome filted by r2 and maf
```
Now all of the chromosome files should be in bfile format. Merge them into one bfile with a command similar to the following:
```
plink --bfile chr1 --merge-list ListOfBfiles.txt --make-bed --out 1000GFilteredPlink
```

From here, there are two steps that we need to do before going further. 1). We need to update the SNP names from chr:position to rsID.  PrediXcan uses rsIDs in its prediction models and these will be necessary going forward.  2) Get the files into PrediXcan dosage format.
```
#This updates the chr:position to rsID
plink --bfile /home/peter/AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/UMichFiltered/1000GFilteredPlink --update-name /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423_no_dups.txt --make-bed --out 1000GrsFiltered
#This converts PLINK files to PrediXcan dosages.
python plinktodosages4predixcan.py --bfile /home/peter/AA_nonGAIN_SCZ/Imputation/UMichResults/1000G/UMich1000G/UMichFiltered/1000GrsFiltered --out /home/peter/AA_nonGAIN_SCZ/PrediXcan/dosages/1000G/1000GdosagesChr
```
The output of these commands are ready for GWAS and PrediXcan, respectively.
Note: [plinktodosages4predixcan.py](https://github.com/hakyim/PrediXcan/blob/master/Software/convert_plink_to_dosage.py) is available in the [Im Lab's PrediXcan github repository](https://github.com/hakyim/PrediXcan). 
