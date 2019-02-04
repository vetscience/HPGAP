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
Pi<-read.table(args[2])
Theta<-read.table(args[3])

harmonicNumber = 0
numChromosomes = as.numeric(args[7])
for (i in 1:(numChromosomes - 1)) {
  harmonicNumber = harmonicNumber + 1.0/i
}
Theta$Theta <- Theta$V3/(as.numeric(args[8])*harmonicNumber)

grob <- grobTree(textGrob(as.character(args[9]), x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))
png(args[5], width = 6, height = 2, units = "in", res = 300)
ggplot() + geom_line(data=Theta, aes(x=V2, y=log10(Theta)),colour="orange") + 
  geom_line(data=Pi, aes(x=V2, y=log10(V5)),colour="red") +
  scale_y_continuous(limits = c(-5, 0),breaks = c(-5,-4,-3,-2,-1,0),labels=c("10-5","10-4",0.001,0.01,0.1,1)) +
  scale_x_continuous(breaks=seq(0,ceiling(max(Theta$V2)/1000000)*1000000,by = 1000000),labels=seq(0,ceiling(max(Theta$V2)/1000000),by = 1)) +
  labs(x = "Position on chromosome (Mb)", y = NULL, title = NULL)+
  annotation_custom(grob)
dev.off()

Tj<-read.table(args[4])
png(args[6], width = 6, height = 2, units = "in", res = 300)
ggplot() + geom_line(data=Tj, aes(x=V2, y=V4)) + 
  geom_hline(aes(yintercept = 0),colour="black", linetype="dashed") + 
  scale_y_continuous(limits = c(-3, 4),breaks = c(-3,-2,-1,0,1,2,3,4),labels=c(-3,-2,-1,0,1,2,3,4)) +
  scale_x_continuous(breaks=seq(0,ceiling(max(Theta$V2)/1000000)*1000000,by = 1000000),labels=seq(0,ceiling(max(Theta$V2)/1000000),by = 1)) +
  labs(x = "Position on chromosome (Mb)", y = NULL, title = NULL)+
  annotation_custom(grob)
dev.off()