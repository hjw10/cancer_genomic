library(tidyverse)
chr_pq=read.table("chr_pq.txt",header = F,sep = "\t",stringsAsFactors = F) 
#该文件由band文件转换得到，参考我的简书文章：染色体区段6p21.31在哪？https://www.jianshu.com/p/477a07192fe6
colnames(chr_pq)=c("chr","arm","cutoff")
chr_pq$chr=paste("chr",chr_pq$chr,sep = "")

cnv_regions=read.table("./try2/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",header = T,sep = "\t",stringsAsFactors = F)
cnv_regions=cnv_regions%>%filter(state!=3)
cnv_regions$cell_group_name=str_replace(cnv_regions$cell_group_name,"all.*observations\\.","")
cnv_regions$cnv_type=""
for (i in 1:nrow(cnv_regions)) {
  tmp_chr_pq=chr_pq%>%filter(chr==cnv_regions[i,"chr"])
  if(cnv_regions[i,"start"] <= tmp_chr_pq[1,"cutoff"] & cnv_regions[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
    cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
  }else if (cnv_regions[i,"start"] >= tmp_chr_pq[2,"cutoff"] & cnv_regions[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
    cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
  }else {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"chr"],"p,",cnv_regions[i,"chr"],"q",sep = "")
  }
  if (cnv_regions[i,"state"] < 3) {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_loss",sep = "")
  }
  if (cnv_regions[i,"state"] > 3) {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_gain",sep = "")
  }
  if (str_detect(cnv_regions[i,"cnv_type"],",")) {
    tmp1=strsplit(cnv_regions[i,"cnv_type"],",")[[1]][1]
    tmp2=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
    tmp3=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
    cnv_regions[i,"cnv_type"]=paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
    rm(list = c("tmp1","tmp2","tmp3"))
  }
}

#CNV 有8类的代码类似
cnv_1.1.1.1=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.1"))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.1.1.2=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.2"))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.1.2  =unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.2"  ))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.2.1.1=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.1"))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.2.1.2=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.2"))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.2.2.1=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.1"))[,"cnv_type"],",",collapse = ""),",")[[1]])
cnv_1.2.2.2=unique(strsplit(paste0(as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.2"))[,"cnv_type"],",",collapse = ""),",")[[1]])

###cnv_1.1.1
cnv_1.1.1=intersect(cnv_1.1.1.1,cnv_1.1.1.2)
cnv_1.1.1.1_uniq=setdiff(cnv_1.1.1.1,cnv_1.1.1)
cnv_1.1.1.2_uniq=setdiff(cnv_1.1.1.2,cnv_1.1.1)
###cnv_1.2.1
cnv_1.2.1=intersect(cnv_1.2.1.1,cnv_1.2.1.2)
cnv_1.2.1.1_uniq=setdiff(cnv_1.2.1.1,cnv_1.2.1)
cnv_1.2.1.2_uniq=setdiff(cnv_1.2.1.2,cnv_1.2.1)
###cnv_1.2.2
cnv_1.2.2=intersect(cnv_1.2.2.1,cnv_1.2.2.2)
cnv_1.2.2.1_uniq=setdiff(cnv_1.2.2.1,cnv_1.2.2)
cnv_1.2.2.2_uniq=setdiff(cnv_1.2.2.2,cnv_1.2.2)
###cnv_1.1
cnv_1.1=intersect(cnv_1.1.1,cnv_1.1.2)
cnv_1.1.1_uniq=setdiff(cnv_1.1.1,cnv_1.1)
cnv_1.1.2_uniq=setdiff(cnv_1.1.2,cnv_1.1)
###cnv_1.2
cnv_1.2=intersect(cnv_1.2.1,cnv_1.2.2)
cnv_1.2.1_uniq=setdiff(cnv_1.2.1,cnv_1.2)
cnv_1.2.2_uniq=setdiff(cnv_1.2.2,cnv_1.2)
###cnv_1
cnv_1=intersect(cnv_1.1,cnv_1.2)
cnv_1.1_uniq=setdiff(cnv_1.1,cnv_1)
cnv_1.2_uniq=setdiff(cnv_1.2,cnv_1)

cat("cnv_1:",cnv_1,"\n",file="out.txt",append = T)
cat("cnv_1.1_uniq:",cnv_1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2_uniq:",cnv_1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1_uniq:",cnv_1.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.2_uniq:",cnv_1.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1_uniq:",cnv_1.2.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2_uniq:",cnv_1.2.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1.1_uniq:",cnv_1.1.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1.2_uniq:",cnv_1.1.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1.1_uniq:",cnv_1.2.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1.2_uniq:",cnv_1.2.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2.1_uniq:",cnv_1.2.2.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2.2_uniq:",cnv_1.2.2.2_uniq,"\n",file="out.txt",append = T)

