plink --bfile hg19_forImputationPrep --exclude Exclude-hg19_forImputationPrep-HRC.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-hg19_forImputationPrep-HRC.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-hg19_forImputationPrep-HRC.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-hg19_forImputationPrep-HRC.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --out hg19_forImputationPrep-updated
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 1 --out hg19_forImputationPrep-updated-chr1
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 2 --out hg19_forImputationPrep-updated-chr2
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 3 --out hg19_forImputationPrep-updated-chr3
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 4 --out hg19_forImputationPrep-updated-chr4
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 5 --out hg19_forImputationPrep-updated-chr5
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 6 --out hg19_forImputationPrep-updated-chr6
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 7 --out hg19_forImputationPrep-updated-chr7
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 8 --out hg19_forImputationPrep-updated-chr8
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 9 --out hg19_forImputationPrep-updated-chr9
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 10 --out hg19_forImputationPrep-updated-chr10
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 11 --out hg19_forImputationPrep-updated-chr11
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 12 --out hg19_forImputationPrep-updated-chr12
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 13 --out hg19_forImputationPrep-updated-chr13
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 14 --out hg19_forImputationPrep-updated-chr14
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 15 --out hg19_forImputationPrep-updated-chr15
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 16 --out hg19_forImputationPrep-updated-chr16
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 17 --out hg19_forImputationPrep-updated-chr17
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 18 --out hg19_forImputationPrep-updated-chr18
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 19 --out hg19_forImputationPrep-updated-chr19
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 20 --out hg19_forImputationPrep-updated-chr20
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 21 --out hg19_forImputationPrep-updated-chr21
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 22 --out hg19_forImputationPrep-updated-chr22
plink --bfile hg19_forImputationPrep-updated --reference-allele Force-Allele1-hg19_forImputationPrep-HRC.txt --make-bed --chr 23 --out hg19_forImputationPrep-updated-chr23
rm TEMP*
