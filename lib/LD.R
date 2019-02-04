#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

list<-read.table(args[2])
dflist <- vector(mode='list', length = length(list$V1))
for(i in list$V1){
  dfr <- read.delim(as.character(i),sep="",header=F,check.names=F,stringsAsFactors=F);
  colnames(dfr) <- c("dist","rsq")
  dfr$lgdist=log10(dfr$dist);dfr$distc <- cut(dfr$lgdist,breaks=seq(from=min(dfr$lgdist)-1,to=max(dfr$lgdist)+1,by=0.05))
  dfr1 <- dfr %>% group_by(distc) %>% dplyr::summarise(mean=mean(rsq),median=median(rsq))
  dflist[[i]]<- dfr1 %>% mutate(start=as.double(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.double(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
}
#grob <- grobTree(textGrob("chrIII", x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))
png(args[3],width = 6, height = 6, units = "in", res = 300)
p<-ggplot();
for(i in list$V1){
  p <- p + geom_line(data=dflist[[i]],aes(x=end,y=mean),size=0.3,alpha=0.9,colour="red")
}
p <- p + labs(x="Distance between sites (log10 bp)",y=expression("Linkage Disequilibrium"~(r^{2}))) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6),labels=c(expression(10^{0}),expression(10^{1}),expression(10^{2}),expression(10^{3}),expression(10^{4}),expression(10^{5}),expression(10^{6})))+
  theme_bw()
print(p);
#annotation_custom(grob)
dev.off()