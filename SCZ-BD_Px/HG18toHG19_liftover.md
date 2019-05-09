#Liftover from hg18 to hg19

Following QC, we have to change the genome build from hg18 to hg19 for impuation with the University of Michican Imputation Server. We know that the build is in hg18 because when looking up SNPs on the UCSC genome browser, the SNPs are only found with their identifier and respective position under hg18.

```
nano newfile
plink --bfile /home/peter/AA_nonGAIN_SCZ/QCSteps/QCStep2/QCStep2 --recode --out /home/peter/AA_nonGAIN_SCZ/liftover/newfile

python LiftMap.py -m newfile.map -p newfile.ped -o new

plink --file new --make-bed --out hg19_forImputationPrep
```
