# Preparing Data to be Uploaded to the University of Michigan Imputation Server

```
plink --bfile /home/peter/AA_nonGAIN_SCZ/liftover/hg19_forImputationPrep --freq --out newfreq
perl /home/peter/AA_GAIN_SCZ/Imputation/preImputation/HRC-1000G-check-bim.pl -b /home/peter/AA_nonGAIN_SCZ/liftover/hg19_forImputationPrep.bim -f newfreq.frq -r all.caapa.sorted.txt -h
```

The command above outputs the script `Run-plink.sh`.

Execute this script so that it performs a series of indpendent PLINK commands on each chromosome/

```
bash plink2vcf.sh #converts from PLINK format to .vcf file
bash vcfzip.sh #zips the .vcf file
```
The .vcf files are now ready to be uploaded for imputation.
