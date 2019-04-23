library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)
"%&%" = function(a, b) paste(a, b, sep = "")


hapmappopinfo<-read.table("z://AA_nonGAIN_SCZ/QCSteps/QCStep6/HAPMAP3_hg18/MyDirectory_pop_HM3_hg18_forPCA.txt")%>%select(V1,V3)
colnames(hapmappopinfo) <- c("pop","IID")
fam<-read.table("z://bipolar_disorder/QC/qcstep6/QCStep6E.fam",header =F)%>%select(V1,V2)
colnames(fam) <- c("FID","IID")
pcs<-read.table("z://bipolar_disorder/QC/qcstep6/QCStep6E.evec",header =F)
eval<-read.table("z://bipolar_disorder/QC/qcstep6/QCStep6E.eval")
popinfo <- left_join(fam,hapmappopinfo,by="IID")
popinfo <-mutate(popinfo, pop=ifelse(is.na(pop),'GWAS',as.character(pop)))
pcdf <- data.frame(popinfo, pcs[,2:11]) %>% rename (PC1=V2,PC2=V3,PC3=V4,PC4=V5,PC5=V6,PC6=V7,PC7=V8,PC8=V9,PC9=V10,PC10=V11)
gwas <- filter(pcdf,pop=='GWAS')
hm3 <- filter(pcdf, grepl('NA',IID))
table(popinfo$pop)


a<- ggplot() + geom_point(data=gwas,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ ggtitle("A") + scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

b<- ggplot() + geom_point(data=gwas,aes(x=PC1,y=PC3,col=pop,shape=pop))+geom_point(data=hm3,aes(x=PC1,y=PC3,col=pop,shape=pop))+ ggtitle("B") + scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

c<- ggplot() + geom_point(data=gwas,aes(x=PC2,y=PC3,col=pop,shape=pop))+geom_point(data=hm3,aes(x=PC2,y=PC3,col=pop,shape=pop))+  ggtitle("C") + scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

grid.arrange(grobs= list(a,b,c), ncol=3, top=textGrob("Bipolar Disorder Study Individuals Plotted Against HapMap Phase 3", gp=gpar(fontsize=20, fontface="bold")))




fam<-read.table("Z://AA_nonGAIN_SCZ/QCSteps/QCStep6/QCStep6E/QCStep6E.fam")%>%select(V1,V2)#This is the last .fam file that was made in the QC process
colnames(fam) <- c("FID","IID")
popinfo <- left_join(fam,hapmappopinfo,by="IID")
popinfo <-mutate(popinfo, pop=ifelse(is.na(pop),'GWAS',as.character(pop)))
table(popinfo$pop)

pcsgenonox0.01<-read.table("Z://AA_nonGAIN_SCZ/QCSteps/QCStep6/QCStep6E/QCStep6E_noX_geno0.01.evec",skip=1)
pcdfgenonox0.01 <- data.frame(popinfo, pcsgenonox0.01[,2:11]) %>% rename (PC1=V2,PC2=V3,PC3=V4,PC4=V5,PC5=V6,PC6=V7,PC7=V8,PC8=V9,PC9=V10,PC10=V11)
gwasgenonox0.01 <- filter(pcdfgenonox0.01,pop=='GWAS')
hm3genonox0.01 <- filter(pcdfgenonox0.01, grepl('NA',IID))
evalgenonox0.01<-scan("Z://AA_nonGAIN_SCZ/QCSteps/QCStep6/QCStep6E/QCStep6E_noX_geno0.01.eval")[1:10]
round(eval/sum(evalgenonox0.01))#Calculate the percent explained by each PC

d <- ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ ggtitle("A") +  scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

e <- ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC1,y=PC3,col=pop,shape=pop))+ geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ ggtitle("B") +  scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

f <-ggplot() + geom_point(data=gwasgenonox0.01,aes(x=PC2,y=PC3,col=pop,shape=pop))+geom_point(data=hm3genonox0.01,aes(x=PC1,y=PC2,col=pop,shape=pop))+ ggtitle("C") +  scale_colour_manual(values = c("gray56", "mediumseagreen", "royalblue3","sandybrown")) + theme_bw(20)

grid.arrange(grobs= list(d,e,f), ncol=3, top=textGrob("Schizophrenia Study Individuals Plotted Against HapMap Phase 3", gp=gpar(fontsize=20, fontface="bold")))