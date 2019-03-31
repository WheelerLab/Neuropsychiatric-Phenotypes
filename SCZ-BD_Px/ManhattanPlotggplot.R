#March 30, 2019

#Making a manhattan Plot using ggplot2

#Notes taken from https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html and 
#https://gist.github.com/slowkow/9041570
library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)

manhattantotal4plot<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/greenflagsnames.txt", header = TRUE, stringsAsFactors = TRUE)
#Read in GTEx V6
tissue<-fread("Z://AA_nonGAIN_SCZ/PrediXcan/output/1000G/AllTissueAssociationGTExVersion7.txt",header = TRUE, stringsAsFactors = TRUE)
#Read in GTEx V7

don <- manhattantotal4plot %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(manhattantotal4plot, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

SCZ_GTExV6<-ggplot(don, aes(x=BPcum, y=-log10(p))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("royalblue3", "sandybrown", "mediumseagreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18," ", 20, " ", 22), breaks= axisdf$center, expand = c(0,0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw(16) +
 xlab("Chromosome") + ylab("-log(P)") + ggtitle("A") +
  theme(legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    ) + facet_wrap(~"GTEx Version 6", nrow=1)+
  # PrediXcan Line
  geom_hline(yintercept = -log10(1.84e-5), color= "red")


don7 <- tissue %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
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
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18," ", 20, " ", 22), breaks= axisdf$center,expand = c(0, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw(16) +
  xlab("Chromosome") + ylab("-log(P)") + ggtitle("B") +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) + facet_wrap(~"GTEx Version 7", nrow=1)+
  # PrediXcan Line
  geom_hline(yintercept = -log10(9.65e-6), color= "red")


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
    geom_hline(yintercept= -log10(1.84e-5), color ="mediumseagreen") +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(16) +
    ggtitle("C")+
    theme(legend.position="none",
          panel.grid=element_blank()
    ) + facet_wrap(~"GTEx Version 6", nrow=1)
}

qqgtex6<-gg_qqplot(manhattantotal4plot$p, ci= 0.95)

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
    geom_hline(yintercept= -log10(9.65e-6), color ="mediumseagreen") +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(16) +
    ggtitle("D")+
    theme(legend.position="none",
          panel.grid=element_blank()
    ) + facet_wrap(~"GTEx Version 7", nrow=1)
}

qqgtex7<-gg_qqplot(tissue$p, ci=0.95)

grid.arrange(grobs= list(SCZ_GTExV6, SCZ_GTExv7, qqgtex6, qqgtex7), nrows=2, top=textGrob("Scizophrenia Gene-Tissue Associations", gp=gpar(fontsize=20, fontface="bold")))

#Making Plots for Bipolar Disorder Data
tissueBP<-fread("Z://bipolar_disorder/predixcan/output/AllTissueAssociationGTEx Version6.txt",header = TRUE, stringsAsFactors = TRUE)

donBP <- tissueBP %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tissueBP, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdfBP = donBP %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

BP_GTEx6<- ggplot(donBP, aes(x=BPcum, y=-log10(p))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("royalblue3", "sandybrown", "mediumseagreen"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16," ",18," ", 20, " ", 22), breaks= axisdfBP$center,expand = c(0, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw(16) +
  xlab("Chromosome") + ylab("-log(P)") + ggtitle("A") +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) + facet_wrap(~"Bipolar Disorder Gene Associations", nrow=1)+
  # PrediXcan Line
  geom_hline(yintercept = -log10(1.84e-5), color= "red")


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
    geom_hline(yintercept= -log10(1.84e-5), color ="mediumseagreen") +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(16) +
    ggtitle("B")+
    theme(legend.position="none",
          panel.grid=element_blank()
    ) + facet_wrap(~"Cerebellum-GTEx Version 6", nrow=1)
}



cerebellum<-fread("Z://bipolar_disorder/predixcan/output/GTEx_V6/1000G_BP_correctphenoTW_Brain_Cerebellum_0.5.db_association.txt", header =T)
row.has.na <- apply(cerebellum, 1, function(x){any(is.na(x))})
sum(row.has.na)

cerebellumnona<-cerebellum[!row.has.na,]
cerebellumqq<-gg_qqplot(cerebellumnona$p, ci =0.95)

grid.arrange(grobs= list(BP_GTEx6,cerebellumqq), ncol=2, top=textGrob("Bipolar Disorder Gene-Tissue Associations", gp=gpar(fontsize=20, fontface="bold")))
