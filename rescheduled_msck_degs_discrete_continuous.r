# 
# setwd("/data/szj/src/r/XEJ/redegseq/input/rescheduled")
# list.files()
# # compare_tcga_lihc_normal_gtex=readRDS(file = "compare_tcga_lihc_normal_gtex.RDS")
# # rownames(compare_tcga_lihc_normal_gtex)=compare_tcga_lihc_normal_gtex$GeneGeneName
# # which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_normal_gtex) ,ignore.case = T,perl = T))
# # compare_tcga_lihc_normal_gtex=compare_tcga_lihc_normal_gtex[,which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_normal_gtex) ,ignore.case = T,perl = T))]
# 
# compare_tcga_lihc_tumor_gtex=readRDS(file = "compare_tcga_lihc_tumor_gtex.RDS")
# rownames(compare_tcga_lihc_tumor_gtex)=compare_tcga_lihc_tumor_gtex$GeneGeneName
# which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_tumor_gtex) ,ignore.case = T,perl = T))
# compare_tcga_lihc_tumor_gtex=compare_tcga_lihc_tumor_gtex[,which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_tumor_gtex) ,ignore.case = T,perl = T))]
# 
# # compare_tcga_lihc_tumor_normal=readRDS(file = "compare_tcga_lihc_tumor_normal.RDS")
# # rownames(compare_tcga_lihc_tumor_normal)=compare_tcga_lihc_tumor_normal$GeneGeneName
# # which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_tumor_normal) ,ignore.case = T,perl = T))
# # compare_tcga_lihc_tumor_normal=compare_tcga_lihc_tumor_normal[,which(grepl(pattern ="^TCGA|^GTEX" ,x =colnames(compare_tcga_lihc_tumor_normal) ,ignore.case = T,perl = T))]
# 
# # str(compare_tcga_lihc_normal_gtex)
# str(compare_tcga_lihc_tumor_gtex)
# # str(compare_tcga_lihc_tumor_normal)
# 
# #查看样本数
# length(which(grepl(pattern = "^TCGA",x =colnames(compare_tcga_lihc_tumor_gtex) ,ignore.case = T,perl = T)))
# #那就分为10份吧，窗口大小为1:30
# colnames(compare_tcga_lihc_tumor_gtex)
# #查看ZNT1的分布
# library(ggpubr)
# data_znt1=as.data.frame(t(compare_tcga_lihc_tumor_gtex[which(grepl(pattern = "^SLC30A1_",x = rownames(compare_tcga_lihc_tumor_gtex),ignore.case = T,perl = T)),]))
# data_znt1[which(grepl(pattern = "^TCGA",x = rownames(data_znt1),ignore.case = T,perl = T)),"GROUP"]="TCGA_LIHC"
# data_znt1[which(grepl(pattern = "^GTEX",x = rownames(data_znt1),ignore.case = T,perl = T)),"GROUP"]="GTEX"
# quantile(unlist(data_znt1$SLC30A1_7779))
# ggdensity(data = data_znt1,x ="SLC30A1_7779" ,color = "GROUP")
# 
# ##对ZNT1列进行拆分
# data_GTEX=compare_tcga_lihc_tumor_gtex[,which(grepl(pattern = "^GTEX",x = colnames(compare_tcga_lihc_tumor_gtex),ignore.case = T,perl = T))]
# data_TCGA=compare_tcga_lihc_tumor_gtex[,which(grepl(pattern = "^TCGA",x = colnames(compare_tcga_lihc_tumor_gtex),ignore.case = T,perl = T))]
# dim(data_GTEX)
# dim(data_TCGA)
# data_znt1_TCGA=as.data.frame(t(compare_tcga_lihc_tumor_gtex[which(grepl(pattern = "^SLC30A1_",x = rownames(compare_tcga_lihc_tumor_gtex),ignore.case = T,perl = T)),which(grepl(pattern = "^TCGA",x = colnames(compare_tcga_lihc_tumor_gtex),ignore.case = T,perl = T))]))
# data_TCGA=data_TCGA[,order(data_znt1_TCGA$SLC30A1_7779)]
# for(colindex in c(1:ncol(data_TCGA))){
#     print(data_TCGA[which(grepl(pattern = "^SLC30A1_",x =rownames(data_TCGA) ,ignore.case = T,perl = T)),colindex])
# }
# 
# ##创建文件夹
setwd("D:/temp/20180805/output/")
list.files()
for(ncomb in c(1:30)){
    if(!dir.exists(paste("_",ncomb,sep = ""))){
        dir.create(paste("_",ncomb,sep = ""))
    }
}
list.dirs()
ls()
# 
# ###数据保存
# save.image(file = "20180805.Rdata")
setwd("D:/temp/20180805/")
load("20180805.Rdata")
rm("a","begin","colindex")
###进行任务组合
tasks=data.frame()
for(width in seq(5,5,by = 1)){
  for(begin in c(1:(ncol(data_TCGA)-width+1))){
    newrowindex=nrow(tasks)+1
    range_window=c(begin:(begin+width-1))
      tasks[newrowindex,"width"]=width
      tasks[newrowindex,"colindex"]=paste(range_window,collapse="_")
      tasks[newrowindex,"outdir"]=paste("_",width,sep="")
  }
}
dim(tasks)

###启动日志系统
# Import the log4r package.
library('log4r')
# Create a new logger object with create.logger().
logger <- create.logger()
# Set the logger's file output.
logfile(logger) <- 'rescheduled_msck_degs_discrete_continuous.log'
# Set the current level of the logger.
level(logger) <- 'INFO'

###开启并行化系统
gc()
library(parallel)
library(LaF)
library(edgeR)
parallel::detectCores()
cls=parallel::makePSOCKcluster(5)
parallel::clusterExport(cl = cls,varlist = ls())

##独立的关闭并行化系统
parallel::stopCluster(cls)

a=c(1,2,3,4,5)
getwd()

###开始并行化计算
parApply(cl = cls,X =tasks,MARGIN = 1,FUN= function(x){
    ##控制流结束--->>>
      library(LaF)
      COMMAND=LaF::get_lines(filename = "D:/temp/20180805/cmd.txt",line_numbers = 1)
        print(COMMAND)
      while(grepl(x = COMMAND,pattern = "SLEEP",perl = T)){
          time=stringi::stri_extract_first(str = COMMAND,regex = "\\d+")
          time=as.numeric(time)
          if(time>1000){
              return("EXIT SWIFTLY!")
          }
          log4r::info(logger = logger,message = paste("Begin to Sleep",time))
          Sys.sleep(time)
        COMMAND=LaF::get_lines(filename = "D:/temp/20180805/cmd.txt",line_numbers = 1)
       }
    ##控制流结束<<<---
    outdir=x["outdir"]
    colindex=x["colindex"]
    log4r::info(logger = logger,message = colindex)
    ##组合matrix
    newMatrix=NULL
    ##选取TCGA_TUMOR的列
    index_selected_column=na.omit(unlist(stringi::stri_split(str = x["colindex"],regex = "_")))
    index_selected_column=as.numeric(index_selected_column)
    newMatrix_TCGA_TUMOR=data_TCGA[,c(index_selected_column)]
    newMatrix_TCGA_TUMOR[,"GeneGeneName"]=rownames(newMatrix_TCGA_TUMOR)
    newMatrix_GTEX=data_GTEX
    newMatrix_GTEX[,"GeneGeneName"]=rownames(newMatrix_GTEX)
    newMatrix=merge(x = newMatrix_GTEX,y =newMatrix_TCGA_TUMOR ,by = "GeneGeneName")
    rownames(newMatrix)=newMatrix[,"GeneGeneName"]
    newMatrix=newMatrix[,-which(grepl(pattern = "GeneGeneName",x = colnames(newMatrix),ignore.case = T,perl = T))]
    log4r::info(logger = logger,message = dim(newMatrix))
    gc()
    ###开始执行edgeR
    library(edgeR)
    group_tumorVSgtex=c()
    for(thecolindex in c(1:ncol(newMatrix))){
        if(grepl(pattern = "^GTEX",x = colnames(newMatrix)[thecolindex],ignore.case = T,perl = T)){
            group_tumorVSgtex[thecolindex]="GTEX"
        }else{
            group_tumorVSgtex[thecolindex]="TCGA_TUMOR"
        }
    }
#     str(group_tumorVSgtex)

    y <- DGEList(counts=newMatrix, group=as.factor(group_tumorVSgtex))
log4r::info(logger = logger,message = "1")####################################################
    thethredshold<-min(summary(as.factor(group_tumorVSgtex)))*0.6
log4r::info(logger = logger,message = "2")####################################################
    keep <- rowSums(cpm(y)>1) >= thethredshold
log4r::info(logger = logger,message = "3")####################################################
    y <- y[keep, , keep.lib.sizes=FALSE]
log4r::info(logger = logger,message = "4")####################################################
    y <- calcNormFactors(y)
log4r::info(logger = logger,message = "5")####################################################
    y <- estimateDisp(y)
log4r::info(logger = logger,message = "6")####################################################
    et <- exactTest(y,pair = c("GTEX","TCGA_TUMOR"))
log4r::info(logger = logger,message = "7")####################################################
    tTag <- topTags(et, n=nrow(y))
log4r::info(logger = logger,message = "8")####################################################
    tTag <- as.data.frame(tTag)
log4r::info(logger = logger,message = "9")###################################################
    all_outdir=paste("D:/temp/20180805/output/",outdir,"/",colindex,"_TCGA_LIHC_TUMOR_vs_GTEX_edgeR.csv",sep = "")
log4r::info(logger = logger,message = "10")####################################################
log4r::info(logger = logger,message = all_outdir)####################################################
    write.csv(tTag,file = all_outdir)
log4r::info(logger = logger,message = "11")####################################################
    log4r::info(logger = logger,message = paste("Write to",all_outdir))
log4r::info(logger = logger,message = "12")####################################################
    gc()
})

# ls()
tasks
# str(data_GTEX)
# str(data_TCGA)

###开始计算TCGA_LIHC VS GTEX的吧
# ls()
# library(edgeR)
# group_tumorVSgtex=c()
# for(colindex in c(1:ncol(compare_tcga_lihc_tumor_gtex))){
#     if(grepl(pattern = "^GTEX",x = colnames(compare_tcga_lihc_tumor_gtex)[colindex],ignore.case = T,perl = T)){
#         group_tumorVSgtex[colindex]="GTEX"
#     }else{
#         group_tumorVSgtex[colindex]="TCGA_TUMOR"
#     }
# }
# str(group_tumorVSgtex)
# y <- DGEList(counts=compare_tcga_lihc_tumor_gtex, group=as.factor(group_tumorVSgtex))
# y$samples
# message(thethredshold<-min(summary(as.factor(group_tumorVSgtex)))*0.6)
# keep <- rowSums(cpm(y)>1) >= thethredshold
# y <- y[keep, , keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# y <- estimateDisp(y)
# et <- exactTest(y,pair = c("GTEX","TCGA_TUMOR"))
# tTag <- topTags(et, n=nrow(y))
# tTag <- as.data.frame(tTag)
# write.csv(tTag,file = "/data/szj/src/r/XEJ/redegseq/output/msck/edgeR/TCGA_LIHC_TUMOR_vs_GTEX_edgeR.csv")
# save.image(file = "/data/szj/src/r/XEJ/redegseq/output/msck/edgeR/compare_tcga_lihc_tumor_gtex.Rdata")
