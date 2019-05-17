library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
"%&%" = function(a,b) paste(a,b,sep="")

pheno<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/phenobinary1.txt", header =F)
pheno$V3[pheno$V3=="1"]<-"SCZ"
pheno$V3[pheno$V3=="0"]<-"Control"
colnames(pheno)<-c("V1","V2","Phenotype")

#SMPD3
lymphocyteXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/1000GTW_Cells_EBV-transformed_lymphocytes_0.5.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000103056.7

lymphocyte_gene<-data.table(FID=lymphocyteXpress$FID,IID=lymphocyteXpress$IID,Xpres=lymphocyteXpress$ENSG00000103056.7)
With_phenolymph<-left_join(lymphocyte_gene , pheno , b = c("IID" = "V2"))
onlywithphenolymph<-dplyr::select(With_phenolymph,FID,IID,Xpres,Phenotype)

#SH3PXD2B
tibialnerveXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/1000GTW_Nerve_Tibial_0.5.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000174705.7

tibialnerve_gene<-data.table(FID=tibialnerveXpress$FID,IID=tibialnerveXpress$IID,Xpres=tibialnerveXpress$ENSG00000174705.7)
With_phenotibialnerve<-left_join(tibialnerve_gene , pheno , b = c("IID" = "V2"))
onlywithphenotibialnerve<-dplyr::select(With_phenotibialnerve,FID,IID,Xpres,Phenotype)

Esophagus_MuscularisXPres<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/1000GTW_Esophagus_Muscularis_0.5.db_predicted_expression.txt", header =T)
#The gene we are looking for is ENSG00000162769.8

#FLVCR1
Esophagus_MuscularisXPres<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/1000GTW_Esophagus_Muscularis_0.5.db_predicted_expression.txt", header =T)
#The gene we are looking for is ENSG00000162769.8

esoph_gene<-data.table(FID=Esophagus_MuscularisXPres$FID,IID=Esophagus_MuscularisXPres$IID,Xpres=Esophagus_MuscularisXPres$ENSG00000162769.8)
With_phenoesoph<-left_join(esoph_gene , pheno , b = c("IID" = "V2"))
onlywithphenoesoph<-dplyr::select(With_phenoesoph,FID,IID,Xpres,Phenotype)

#RP11-645c24.5
colonXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Colon_Sigmoid_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000260306.1

colon_gene<-data.table(FID=colonXpress$FID,IID=colonXpress$IID,Xpres=colonXpress$ENSG00000260306.1)
With_phenocolon<-left_join(colon_gene , pheno , b = c("IID" = "V2"))
onlywithphenocolon<-dplyr::select(With_phenocolon,FID,IID,Xpres,Phenotype)

pancXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Pancreas_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000260306.1

panc_gene<-data.table(FID=pancXpress$FID,IID=pancXpress$IID,Xpres=pancXpress$ENSG00000260306.1)
With_phenopanc<-left_join(panc_gene , pheno , b = c("IID" = "V2"))
onlywithphenopanc<-dplyr::select(With_phenopanc,FID,IID,Xpres,Phenotype)

#PRMT7
pitXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Pituitary_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000132600.12

pit_gene<-data.table(FID=pitXpress$FID,IID=pitXpress$IID,Xpres=pitXpress$ENSG00000132600.12)
With_phenopit<-left_join(pit_gene , pheno , b = c("IID" = "V2"))
onlywithphenopit<-dplyr::select(With_phenopit,FID,IID,Xpres,Phenotype)

#ZBTB1
stom7Xpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Stomach_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000126804.9

stom7_gene<-data.table(FID=stom7Xpress$FID,IID=stom7Xpress$IID,Xpres=stom7Xpress$ENSG00000126804.9)
With_phenostom7<-left_join(stom7_gene , pheno , b = c("IID" = "V2"))
onlywithphenostom7<-dplyr::select(With_phenostom7,FID,IID,Xpres,Phenotype)

#AC093843.1
ThyroidXpress<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/v7_GTEx/1000G_gtex_v7_Thyroid_imputed_europeans_tw_0.5_signif.db_predicted_expression.txt",header = T)
#The gene we are interested in is ENSG00000224819.1

thyroid_gene<-data.table(FID=thyroidXpress$FID,IID=thyroidXpress$IID,Xpres=thyroidXpress$ENSG00000224819.1)
With_phenothyroid<-left_join(thyroid_gene , pheno , b = c("IID" = "V2"))

#ZNF562
phenoBP<-fread("Z://bipolar_disorder/predixcan/pheno.txt", header =F)
phenoBP$V3[phenoBP$V3=="0"]<-"Control"
phenoBP$V3[phenoBP$V3=="1"]<-"BD"
colnames(phenoBP)<-c("V1","V2","Phenotype")

cerebellum_gene<-data.table(FID=cerebellumXpress$FID,IID=cerebellumXpress$IID,Xpres=cerebellumXpress$ENSG00000171466.5)
With_phenocerebellum<-left_join(cerebellum_gene , phenoBP , b = c("IID" = "V2"))
onlywithphenocerebellum<-dplyr::select(With_phenocerebellum,FID,IID,Xpres,Phenotype)
ggplot(data = onlywithphenocerebellum, aes(x=as.character(Phenotype),y=Xpres)) + geom_boxplot(fill='navy', color="purple") + ggtitle("ZNF562") +xlab("Case/Control Status") + ylab("Predicted Gene Expression")
onlywithphenocerebellum$Phenotype<- factor(onlywithphenocerebellum$Phenotype)

ZNF562 <-ggplot(data = onlywithphenocerebellum, aes(x=reorder(as.character(Phenotype),-Xpres),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(ZNF562title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("mediumseagreen","royalblue3"))) + theme_bw()

smpd3title<- expression(paste(italic("SMPD3"), " (Lymphocytes v6)"))
SMPD3<-ggplot(data = onlywithphenolymph, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(smpd3title) +xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

sh3pxd2btitle<- expression(paste(italic("SH3PXD2B"), " (Tibial Nerve v6)"))
SH3PXD2B<-ggplot(data = onlywithphenotibialnerve, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(sh3pxd2btitle) +xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

RP11.645c24.5title<- expression(paste(italic("RP11.645c24.5"), " (Sigmoid Colon v7)"))
RP11.645c24.5 <-ggplot(data = onlywithphenocolon, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(RP11.645c24.5title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

RP11.645c24.5_panc_title<-expression(paste(italic("RP11.645c24.5"), " (Pancreas v7)"))
RP11.645c24.5_panc <-ggplot(data = onlywithphenopanc, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(RP11.645c24.5_panc_title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

PRMT7title<-expression(paste(italic("PRMT7"), " (Heart Atrial Appendage v7)"))
PRMT7 <-ggplot(data = onlywithphenopit, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(PRMT7title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

ZBTB1title<-expression(paste(italic("ZBTB1"), " (Stomach v7)"))
ZBTB1_v7 <-ggplot(data = onlywithphenostom7, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(ZBTB1title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

AC093843.1title<-expression(paste(italic("AC093843.1"), " (Thyroid v7)"))
AC093843.1 <-ggplot(data = onlywithphenothyroid, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(AC093843.1title)+ xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

FLVCR1title<-expression(paste(italic("FLVCR1"), " (Esophagus Muscularis v6)"))
FLVCR1<-ggplot(data = onlywithphenoesoph, aes(x=as.character(Phenotype),y=Xpres, fill=Phenotype)) +
  geom_violin() +geom_boxplot(fill="white", width = .125) + ggtitle(FLVCR1title) + xlab("Phenotype") + ylab("Predicted Gene Expression")+ scale_fill_manual(values=alpha(c("royalblue3","sandybrown"))) + theme_bw()

my.legend<-g_legend(FLVCR1)

onlywithphenothyroid<-dplyr::select(With_phenothyroid,FID,IID,Xpres,Phenotype)


grid.arrange(grobs= list(SMPD3 + theme(legend.position="none", axis.title.y = element_blank()),
                         SH3PXD2B + theme(legend.position="none", axis.title.y = element_blank()),
                         FLVCR1 + theme(legend.position="none", axis.title.y = element_blank()),
                         PRMT7+ theme(legend.position="none", axis.title.y = element_blank()),
                         AC093843.1 + theme(legend.position="none", axis.title.y = element_blank()),
                         ZBTB1_v7 + theme(legend.position="none", axis.title.y = element_blank()), 
                         RP11.645c24.5 + theme(legend.position="none", axis.title.y = element_blank()),
                         RP11.645c24.5_panc + theme(legend.position="none", axis.title.y = element_blank()), 
                         ZNF562 + theme(legend.position="none", axis.title.y = element_blank())),
             nrows=3,
             left = textGrob("Predicted Gene Expression", rot = 90, vjust = 1, gp=gpar(fontsize=16)), my.legend)