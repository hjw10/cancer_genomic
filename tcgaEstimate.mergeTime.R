setwd("E:/1-生信上课/肿瘤基因组学/实验课/TCGA-WEB/5time")
score2<-read.table("scores.txt",header = T,sep ="\t",stringsAsFactors=F)
rn<-dim(score2)[1]
library(stringi)
library(stringr)
for (i in 1:rn){
  score2$ID2[i]<-str_sub(score2$ID[i],start = 1,end=12)
}
cl_df_frame_zong$ID2<-as.character(cl_df_frame_zong$ID2)
colnames(cl_df_frame_zong)[1]<-"ID2"
library(dplyr)
#left_join在dplyr包
score_clinc<-left_join(x=score2,y=cl_df_frame_zong,by="ID2")
#接下来就要过滤掉临床信息缺失的患者

score_clinc2<-filter(score_clinc,futime!="unkown")

save.image("tcgaestimatemergetime.RData")

