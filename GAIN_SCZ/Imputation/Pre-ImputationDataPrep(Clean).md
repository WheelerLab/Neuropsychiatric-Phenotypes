---
title: "Pre-Imputation Data Prep (Clean)"
author: "Peter Fiorica"
date: "December 5, 2017"
output: html_document
---

```
~$plink --bfile /home/peter/Documents/Imputation/preImputation/QCStep2Liftover/hg19 --freq --out newfreq

```

The next .pl document will was copied from the wheelab1 directory (path: /home/wheelerlab1/Data/preImputation-check/HRC-1000G-check-bim.pl)
```
~$perl HRC-1000G-check-bim.pl -b /home/peter/Documents/Imputation/preImputation/QCStep2Liftover/hg19.bim -f newfreq.frq -r /home/wheelerlab1/Data/preImputation-check/all.caapa.sorted.txt -h
```
In the above command, we made the file "Run-plink.sh". This file contains every plink command that we need for the next command.
```
~$bash Run-plink.sh
```
The command above executes a series of plink commands contained in "Run-plink.sh".  That command will output bfiles for chromosomes 1-22.  Once these bfiles have been created, we want to recode them to .vcf files.  The following command will do that for us.
```
$plink --bfile hg19-updated-chr1 --recode vcf --out hg19-updated-chr1
##Execute this command on each chromsome up and including 22
$plink --bfile hg19-updated-chr22 --recode vcf --out hg19-updated-chr22
```
Since the .vcf files are large, we want to zip and sort them before we upload them.  The following command will zip and sort them:
```
~$vcf-sort hg19-updated-chr22.vcf | bgzip -c > hg19-updated-chr22.vcf.gz
##Execute this command on each chromsome up to and including 22
~$vcf-sort hg19-updated-chr22.vcf | bgzip -c > hg19-updated-chr22.vcf.gz
```
