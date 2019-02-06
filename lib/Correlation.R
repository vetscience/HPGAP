args = commandArgs(trailingOnly=TRUE)

args[2]="/home/darcy/PopGen_WorkFlow/Example/05.IntraPopulation/Slidingwindow/final.South.allstat.txt"
setwd=args[1]
input=args[2]

df<-read.table(input, header=TRUE)
#df<-na.omit(df)


piVSgc_content<-c(cor(df$pi,df$gc_content,method="spearman",use = "complete.obs"),
                  cor.test(df$pi,df$gc_content,method="spearman",use = "complete.obs")$p.value)
piVScoding_percentage<-c(cor(df$pi,df$coding_percentage,method="spearman",use = "complete.obs"),
                         cor.test(df$pi,df$coding_percentage,method="spearman",use = "complete.obs")$p.value)

piVSr2<-c(cor(df$pi,df$r2,method="spearman",use = "complete.obs"),
          cor.test(df$pi,df$r2,method="spearman",use = "complete.obs")$p.value)
piVSpnps<-c(cor(df$pi,df$pnps,method="spearman",use = "complete.obs"),
            cor.test(df$pi,df$pnps,method="spearman",use = "complete.obs")$p.value)
piVStheta<-c(cor(df$pi,df$theta,method="spearman",use = "complete.obs"),
            cor.test(df$pi,df$theta,method="spearman",use = "complete.obs")$p.value)
piVStajimsD<-c(cor(df$pi,df$tajimsD,method="spearman",use = "complete.obs"),
             cor.test(df$pi,df$tajimsD,method="spearman",use = "complete.obs")$p.value)
r2VStajimsD<-c(cor(df$r2,df$tajimsD,method="spearman",use = "complete.obs"),
               cor.test(df$r2,df$tajimsD,method="spearman",use = "complete.obs")$p.value)
psVSpnps<-c(cor(df$ps,df$pnps,method="spearman",use = "complete.obs"),
               cor.test(df$ps,df$pnps,method="spearman",use = "complete.obs")$p.value)
stat.result<-rbind(
  summary(na.omit(df$gc_content)),
  summary(na.omit(df$coding_percentage)),
  summary(na.omit(df$pi)),
  summary(na.omit(df$theta)),
  summary(na.omit(df$tajimsD)),
  summary(na.omit(df$r2)),
  summary(na.omit(df$ps)),
  summary(na.omit(df$pn)),
  summary(na.omit(df$pnps))
)

rownames(stat.result)<-c("gc_content","coding_percentage","pi","theta","tajimsD","r2","ps","pn","pnps")
colnames(stat.result) <- c("Min", "1st Qu",  "Median",  "Mean", "3rd Qu", "Max")
cor.result<-rbind(piVSgc_content,
                  piVScoding_percentage,
                  piVSr2,
                  piVSpnps,
                  piVStheta,
                  piVStajimsD,
                  r2VStajimsD,
                  psVSpnps)

colnames(cor.result)<-c("roh","p-value")
print(stat.result)
print(cor.result)

x <- c(input, ".stat.table") 
stat.result.file<-paste(x, collapse = "")
x <- c(input, ".cor.table") 
cor.result.file<-paste (x, collapse = "")
write.table(stat.result,row.names=TRUE,col.names=NA,quote=FALSE,file = stat.result.file,sep="\t")
write.table(cor.result,row.names=TRUE,col.names=NA,quote=FALSE,file = cor.result.file,sep="\t")

