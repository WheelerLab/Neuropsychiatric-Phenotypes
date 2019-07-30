# Running MulTiXcan

Clone this [directory](https://github.com/hakyimlab/MetaXcan) from the Im lab, and the scripts needed are inside.

```
python MulTiXcan.py 
--expression_folder /home/peter/AA_nonGAIN_SCZ/PrediXcan/output/1000G/ #Folder with all of you expression files (MESA, GTEx v6, GTEx v7, etc.) of predictors
--expression_pattern "1000GTW_(.*)_0.5.db_predicted_expression.txt" #Pattern for predicted expression files, so that Multixcan can go through all of the tissues efficiently
--input_phenos_file /home/peter/AA_nonGAIN_SCZ/PrediXcan/titleresidualphenoMGS.txt  #Phenotype file.  It had a header as FID 	IID 	Pheno
--input_phenos_column Pheno  #Specify the phenotype column header 
--mode linear #This is the regression type
--output /home/peter/AA_nonGAIN_SCZ/MulTiXcanOutputnonGAIN #This is your output file.  You only get one output file.
```

MultiXcan will not match your predicted expression IIDs and phenotypes automatically.
This means that if you have your individuals ordered as 1, 2, 3, 4, 5 in your expression file,
you need to have the same in your phenotype file.
You unfortunately can't have extra individuals in either your phenotype or gene expression files.
