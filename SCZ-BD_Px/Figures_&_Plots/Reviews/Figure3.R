#Making heatmap for PRMT7 in GTEx 7
#Peter Fiorica
#18 July 2019

library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)

GTEv7.1<-fread("Z://PrediXcan/tissuelistV7",header =F)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub(".{35}$","",GTEv7.1$V1)
GTEv7<-GTEv7.1$V1

PRMT7_pgc<-fread("Z://AA_nonGAIN_SCZ/PeerJReviews/PRMT7_spredixcanGTExV7.txt", header = T)
siggenes<-fread("C://Users/Peter Fiorica/Documents/Wheeler Lab/siggenes4heatmap.txt", header =T)
PRMT7_gain<-subset(siggenes,genename=="PRMT7")


PRMT7_GAIN<-dplyr::select(PRMT7_gain, t,p,tissue)
PRMT7_GAIN$data<-"GAIN"
PRMT7_PGC<-dplyr::select(PRMT7_pgc, zscore, pvalue ,tissue)
PRMT7_PGC$data<- "PGC"
colnames(PRMT7_PGC)<-colnames(PRMT7_GAIN)
total<-bind_rows(PRMT7_GAIN,PRMT7_PGC)

t<-c(NA,NA,NA,NA,NA)
p<-c(NA,NA,NA,NA,NA)
tissue<-c("Brain_Amygdala","Brain_Spinal_cord_cervical_c-1" , "Brain_Substantia_nigra", "Ovary", "Uterus")
data<-c("GAIN", "GAIN","PGC","PGC", "PGC")
for_append<-data.table(t,p,tissue,data)

total<-bind_rows(total,for_append)



#PRMT7_pvalue<-ggplot(total)+
             aes(x = tissue, y = data, fill = -log10(p))+
             geom_raster() +
             scale_fill_gradient2(na.value = NA, high="royalblue3", mid="aliceblue",name = "-log(P)", limits=c(0,8.5), breaks=c(2,4,6,8)) + 
             theme_bw() +
             ggtitle("A")+
             theme(panel.grid.major.y =element_blank(),
                   axis.title= element_text(size=14,face="bold"),
                   axis.text.x = element_text(angle=90, hjust = 1, vjust=.05),
                   plot.title=element_text(size = 14, face="bold")
             )+   #facet_grid(. ~ "P-Value") + 
             #theme(strip.text.x = element_text(size = 15))+
             ylab("Data Source") +
             xlab("Tissue")

#PRMT7_test<- ggplot(total)+
  aes(x = tissue, y = data, fill = t) +
  geom_raster() +
  scale_fill_gradient2(na.value = NA, high="royalblue3", mid="white", low = "red",name = "T-Score" ) + 
  theme_bw() +
  ggtitle("B")+
  theme(panel.grid.major.y =element_blank(),
        axis.title= element_text(size=14,face="bold"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=.05),
        plot.title=element_text(size = 14, face="bold")
  )+   #facet_grid(. ~ "T-Score") + 
  #theme(strip.text.x = element_text(size = 15))+
  ylab("Data Source") +
  xlab("Tissue")

#grid.arrange(grobs= list(PRMT7_pvalue, PRMT7_test),nrow=2)


tiff("Z://AA_nonGAIN_SCZ/PeerJReviews/plots/Bubble_plot_Figure3.tiff", width=1200, height = 500,compression = c("lzw"))

ggplot(data=total, aes(y=data, x = tissue, size = -log(total$p), color = total$t))+
  geom_point(alpha=0.6)+
  scale_y_discrete(limits=c("PGC","GAIN"))+
  scale_size_continuous(range =c(0,15), breaks=c(2,4,6,8), name = "-log(P)")+
  scale_color_gradient2(na.value = NA, high="royalblue3", mid="white", low = "red",name = "T-Score" ) + 
  theme_bw()+
  theme(panel.grid.major.y =element_blank(),
        axis.title= element_text(size=14,face="bold"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=.05),
        plot.title=element_text(size = 14, face="bold", hjust = 0.5)
        )+ 
  ylab("Data Source") +
  xlab("Tissue")

dev.off()
