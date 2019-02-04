#!/usr/bin/env Rscript
library(ggplot2)
library(grid)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

setwd(args[1])

D<-read.csv(args[2],head=FALSE)
#D$AdPi<-D$V6*D$V5/(D$V3-D$V2+1)
D$AdDxy<-D$V8*D$V5/(D$V3-D$V2+1)
D$AdFst<-D$V9*D$V5/(D$V3-D$V2+1)

grob <- grobTree(textGrob(as.character(args[5]), x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))

png(args[3], width = 6, height = 2, units = "in", res = 300)
ggplot() + geom_line(data=D, aes(x=V4, y=log10(AdDxy)),colour="orange") + 
  geom_line(data=D, aes(x=V4, y=log10(AdFst)),colour="red") +
  scale_y_continuous(limits = c(-3, -1),breaks = c(-3,-2,-1),labels=c(0.001,0.01,0.1)) +
  scale_x_continuous(breaks=seq(0,ceiling(max(D$V3)/1000000)*1000000,by = 1000000),labels=seq(0,ceiling(max(D$V3)/1000000),by = 1)) +
  labs(x = "Position on chromosome (Mb)", y = NULL, title = NULL) +
  annotation_custom(grob)
dev.off()