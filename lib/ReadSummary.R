library(ggplot2)
library(gridExtra)
setwd("/home/darcy/PopGen_WorkFlow/Cb_Validation/01.QualityControl/Reports")
sample <- "JU1341"
dir <- "/home/darcy/PopGen_WorkFlow/Cb_Validation/01.QualityControl/Reports/Samples/"

samplelist<-list()
maplist<-list()
for (sample in list.files(dir)){
  print(sample)
  ############ Length################
  r <- "/(\\w+)_lhist"
  lib.id<-gsub(r, "\\1", regmatches(Sys.glob(paste(dir, sample, "/*lhist.txt", sep="")),gregexpr(r,Sys.glob(paste(dir, sample, "/*lhist.txt", sep="")))))
  raw.lengthList <- lapply(Sys.glob(paste(dir, sample, "/*lhist.txt", sep="")), read.table)
  raw.No.reads<-sapply(1:length(raw.lengthList[]),function(x) sum(raw.lengthList[[x]]$V2))
  raw.No.bases<-sapply(1:length(raw.lengthList[]),function(x) sum(as.numeric(raw.lengthList[[x]]$V2)*as.numeric(raw.lengthList[[x]]$V1)))
  #raw.No.reads <- c(raw.No.reads, sum(raw.No.reads)); raw.No.bases <- c(raw.No.bases, sum(raw.No.bases))

  lengthList <- lapply(Sys.glob(paste(dir, sample, "/*lhist.filt.txt", sep="")), read.table)
  No.reads<-sapply(1:length(lengthList[]),function(x) sum(lengthList[[x]]$V2))
  No.bases<-sapply(1:length(lengthList[]),function(x) sum(as.numeric(lengthList[[x]]$V2)*as.numeric(lengthList[[x]]$V1)))
  #No.reads <- c(No.reads, sum(No.reads)); No.bases <- c(No.bases, sum(No.bases))
  df<-data.frame(raw.No.reads,raw.No.bases,No.reads,No.bases)
  mylist<-list(lib.id=lib.id,raw.No.reads=raw.No.reads,raw.No.bases=raw.No.bases,No.reads=No.reads,No.bases=No.bases)

  ############Distribution of Phred Quality#########
  raw.aqhistList <- lapply(Sys.glob(paste(dir, sample, "/*_aqhist.txt", sep="")), read.table)
  raw.aq.count.sum.r1<-sapply(1:length(raw.aqhistList[]),FUN = function(x) sum(as.numeric(raw.aqhistList[[x]]$V2)))
  raw.aq.count.sum.r2<-sapply(1:length(raw.aqhistList[]),FUN = function(x) sum(as.numeric(raw.aqhistList[[x]]$V4)))
  raw.aq.sum.r1<-sapply(1:length(raw.aqhistList[]),FUN = function(x) sum(as.numeric(raw.aqhistList[[x]]$V2)*as.numeric(raw.aqhistList[[x]]$V1)))
  raw.aq.sum.r2<-sapply(1:length(raw.aqhistList[]),FUN = function(x) sum(as.numeric(raw.aqhistList[[x]]$V4)*as.numeric(raw.aqhistList[[x]]$V1)))
  raw.aq.avg.r1<-raw.aq.sum.r1/raw.aq.count.sum.r1
  raw.aq.avg.r2<-raw.aq.sum.r2/raw.aq.count.sum.r2
  #raw.aq.avg <- paste(as.character(round(raw.aq.avg.r1,2)),as.character(round(raw.aq.avg.r2,2)),sep = "/")
  mylist[["raw.aq.count.sum.r1"]]<-raw.aq.count.sum.r1; mylist[["raw.aq.count.sum.r2"]]<-raw.aq.count.sum.r1
  mylist[["raw.aq.avg.r1"]]<-raw.aq.avg.r1; mylist[["raw.aq.avg.r2"]]<-raw.aq.avg.r2
  
  aqhistList <- lapply(Sys.glob(paste(dir, sample, "/*_aqhist.filt.txt", sep="")), read.table)
  aq.count.sum.r1<-sapply(1:length(aqhistList[]),FUN = function(x) sum(as.numeric(aqhistList[[x]]$V2)))
  aq.count.sum.r2<-sapply(1:length(aqhistList[]),FUN = function(x) sum(as.numeric(aqhistList[[x]]$V4)))
  aq.sum.r1<-sapply(1:length(aqhistList[]),FUN = function(x) sum(as.numeric(aqhistList[[x]]$V2)*as.numeric(aqhistList[[x]]$V1)))
  aq.sum.r2<-sapply(1:length(aqhistList[]),FUN = function(x) sum(as.numeric(aqhistList[[x]]$V4)*as.numeric(aqhistList[[x]]$V1)))
  aq.avg.r1<-aq.sum.r1/aq.count.sum.r1
  aq.avg.r2<-aq.sum.r2/aq.count.sum.r2

  #aq.avg <- paste(as.character(round(aq.avg.r1,2)),as.character(round(aq.avg.r2,2)),sep = "/")
  mylist[["aq.count.sum.r1"]]<-aq.count.sum.r1; mylist[["aq.count.sum.r2"]]<-aq.count.sum.r1
  mylist[["aq.avg.r1"]]<-aq.avg.r1; mylist[["aq.avg.r2"]]<-aq.avg.r2

  ############Distribution of GC content#########
  raw.gchistList <- lapply(Sys.glob(paste(dir, sample, "/*_gchist.txt", sep="")), read.table)
  raw.gc.sum<-sapply(1:length(raw.gchistList[]),FUN = function(x) sum(as.numeric(raw.gchistList[[x]]$V2)*as.numeric(raw.gchistList[[x]]$V1)))
  raw.gc.count.sum<-sapply(1:length(raw.gchistList[]),FUN = function(x) sum(as.numeric(raw.gchistList[[x]]$V2)))
  raw.gc.avg<-raw.gc.sum/raw.gc.count.sum
  mylist[["raw.gc.avg"]]<-raw.gc.avg

  gchistList <- lapply(Sys.glob(paste(dir, sample, "/*_gchist.filt.txt", sep="")), read.table)
  gc.sum<-sapply(1:length(gchistList[]),FUN = function(x) sum(as.numeric(gchistList[[x]]$V2)*as.numeric(gchistList[[x]]$V1)))
  gc.count.sum<-sapply(1:length(gchistList[]),FUN = function(x) sum(as.numeric(gchistList[[x]]$V2)))
  gc.avg<-gc.sum/gc.count.sum
  mylist[["gc.avg"]]<-gc.avg
  
  # turn mylist into dataframe
  mydf<-do.call(cbind.data.frame, mylist)
  # drop some colomns for simplicity and attach the dataframe to sample list
  drops <- c("raw.aq.count.sum.r1","raw.aq.count.sum.r2","aq.count.sum.r1","aq.count.sum.r2")
  samplelist[[sample]]<-mydf[ , !(names(mydf) %in% drops)]
  
  COV<-read.table(paste(dir, sample,"/COV.stat.txt", sep=""))
  #Avg.cov<-round(sum(as.numeric(COV$V3)*COV$V2)/sum(as.numeric(COV$V3)))
  #Reference.length<-562768885
  #Coverage<-round(sum(as.numeric(COV$V3))/Reference.length*100,2)
  #Coverage10<-round(sum(as.numeric(COV$V3[10:length(COV$V3)]))/Reference.length*100,2)

  ## mapping list
  SN <- read.table(paste(dir, sample, "/SN.stat.txt", sep=""))
  rownames(SN)<-SN$V1
  mapdf<-data.frame("sample.id"=sample,
                    "raw.data"=sum(mydf$raw.No.bases),
                    "filtered.data"=sum(mydf$No.bases),
                    "reads.mapped"= SN["readsmapped:",2],
                    "Percentage.mapped"=round(as.numeric((SN["readsmapped:",2])/sum(mydf$No.reads)*100),2),
                    "error.rate"=SN["errorrate:",2],
                    "insertsize.avg"=SN["insertsizeaverage:",2],
                    "insertsize.std"=SN["insertsizestandarddeviation:",2],
                    "covered.area"=round(sum(as.numeric(COV$V3))),
                    "covered.area.avg"=round(sum(as.numeric(COV$V3)*COV$V2)/sum(as.numeric(COV$V3)))
                    )
  maplist[[sample]]<-mapdf
  
  alldf<-data.frame("sample.id"=sample,
                    "raw.data"=sum(mydf$raw.No.bases),
                    "filtered.data"=sum(mydf$No.bases),
                    "reads.mapped"= SN["readsmapped:",2],
                    "Percentage.mapped"=round(as.numeric((SN["readsmapped:",2])/sum(mydf$No.reads)*100),2),
                    "error.rate"=SN["errorrate:",2],
                    "insertsize.avg"=SN["insertsizeaverage:",2],
                    "insertsize.std"=SN["insertsizestandarddeviation:",2],
                    "covered.area"=round(sum(as.numeric(COV$V3))),
                    "covered.area.avg"=round(sum(as.numeric(COV$V3)*COV$V2)/sum(as.numeric(COV$V3)))
  )
  alllist[[sample]]<-mapdf
}

alldf<-do.call("rbind", maplist)
write.table(alldf, "alldf.stat.xls", sep="\t")
