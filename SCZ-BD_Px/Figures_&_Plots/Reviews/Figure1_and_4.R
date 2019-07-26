#March 30, 2019

#Making a manhattan Plot using ggplot2

#Notes taken from https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html and 
#https://gist.github.com/slowkow/9041570
library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrepel)

"%&%" = function(a,b) paste(a,b,sep="")


tissue<-fread("Z://wl3backup/AllTissueAssociationGTExVersion7.txt",header = TRUE, stringsAsFactors = TRUE)
#Read in GTEx V7

don7 <- tissue %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tissue, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf7 = don7 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

SCZ_GTExv7<- ggplot(don7, aes(x=BPcum, y=-log10(p))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("royalblue3", "sandybrown", "mediumseagreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18," ", 20, " ", 22), breaks= axisdf7$center,expand = c(0, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw(16) +
  xlab("Chromosome") + ylab("-log(P)") + ggtitle("A") +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) + 
  geom_label_repel(aes(label=ifelse(don7$p < 6.0e-06 & don7$p>5.480e-06, as.character(don7$gene_name),"")),box.padding = unit(0.35, "lines"))

gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3, color = "royalblue3") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color ="Sandybrown") +
    geom_line(aes(expected, clower), linetype = 2, color ="sandybrown") +
    #geom_hline(yintercept= -log10(9.65e-6), color ="mediumseagreen") +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(16) +
    ggtitle("B")+
    theme(legend.position="none",
          panel.grid=element_blank()
    ) + 
    geom_label_repel(aes(x=expected, y=observed, label=ifelse(observed==observed[3], "PRMT7","")),nudge_x=-0.25, box.padding = unit(0.35, "lines"))
}

qqgtex7<-gg_qqplot(tissue$p, ci=0.95)


#tiff("Z://AA_nonGAIN_SCZ/PeerJReviews/plots/SCZ_man7_Figure1.tiff", width=1250, height = 625,compression = c("lzw"))
#grid.arrange(grobs= list(SCZ_GTExv7, qqgtex7), ncols=2, nrow=1)
#dev.off()

png("Z://wl3backup/SCZ_man7_Figure1.png", width=1250, height = 625)
grid.arrange(grobs= list(SCZ_GTExv7, qqgtex7), ncols=2, nrow=1)
dev.off()

GTEv7.1<-fread("Z://wl3backup/tissuelistV7",header =F)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub("^.*?_","",GTEv7.1$V1)
GTEv7.1$V1<-gsub(".{35}$","",GTEv7.1$V1)
GTEv7<-GTEv7.1$V1

genenames<-fread("Z://BP_Chrome.txt", header = T)

for (i in GTEv7){
  a<-fread("Z://wl3backup/bip_px_gtex7/1000G_BP_V7_correctphenogtex_v7_"%&%i%&%"_imputed_europeans_tw_0.5_signif.db_association.txt", header = T)
  a$tissue<-i
  if(exists("alltiss")){
    alltiss<-bind_rows(alltiss, a)
  }else{
    alltiss<-a
  }
}

allw_names<-left_join(alltiss,genenames, by = "gene")
allw_names<-allw_names[complete.cases(allw_names),]

fwrite(allw_names, "Z://all_gtex7_results.txt", col.names=T, row.names = F, sep="\t", quote = F)
#Making Plots for Bipolar Disorder Data
tissueBP<-fread("Z://all_gtex7_results.txt",header = TRUE, stringsAsFactors = TRUE)

donBP <- tissueBP %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tissueBP, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdfBP = donBP %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

BP_GTEx7<- ggplot(donBP, aes(x=BPcum, y=-log10(p))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("royalblue3", "sandybrown", "mediumseagreen"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16," ",18," ", 20, " ", 22), breaks= axisdfBP$center,expand = c(0, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,5.5)) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw(16) +
  xlab("Chromosome") + ylab("-log(P)") + ggtitle("A") +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) #+ facet_wrap(~"Bipolar Disorder Gene Associations", nrow=1)
  # PrediXcan Line
  #geom_hline(yintercept = -log10(1.84e-5), color= "red")

#tiff(filename = "Z://AA_nonGAIN_SCZ/PeerJReviews/plots/BD_manhattan_Figure4.tiff",750,750)
#BP_GTEx7
#dev.off()


gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3, color = "royalblue3") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color ="Sandybrown") +
    geom_line(aes(expected, clower), linetype = 2, color ="sandybrown") +
    #geom_hline(yintercept= -log10(1.84e-5), color ="mediumseagreen") +
    scale_y_continuous(limits = c(0,5.25) ) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(16) +
    ggtitle("B")+
    theme(legend.position="none",
          panel.grid=element_blank()
    ) #+ facet_wrap(~"Cerebellum-GTEx Version 6", nrow=1)
}



bpqq<-gg_qqplot(tissueBP$p, ci =0.95)


png("Z://wl3backup/BD_man7_Figure4.png",width=1250, height = 625 )
grid.arrange(grobs= list(BP_GTEx7,bpqq), ncol=2)
dev.off()
