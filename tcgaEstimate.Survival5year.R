setwd("E:/1-生信上课/肿瘤基因组学/实验课/TCGA-WEB/time")
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
                                
library(survival)

rt$futime=rt$futime/365                                        #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)
  #pValue=format(pValue, scientific = TRUE)
  
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)
  
  pdf(file=paste(gene,".survival-5years.pdf",sep=""),
      width=6,
      height=6)
  plot(fit, 
       lwd=2,
       col=c("red","blue"),
       xlab="Time (year)",
       #加入xlim参数,此步关键
       xlim=c(0,5),
       mark.time=T,
       ylab="Survival rate",
       main=paste(gene,"(p=", pValue ,")",sep="") )
  legend("topright", 
         c("High","Low"), 
         lwd=2, 
         col=c("red","blue"))
  dev.off()
}
write.table(outTab,file="survival-5years.xls",sep="\t",row.names=F,quote=F)





#####求真实的P值
#setwd("E:/1-生信上课/肿瘤基因组学/实验课/TCGA-WEB/time")
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
library(survival)                                 

rt$futime=rt$futime/365   #如果以月为单位，除以30；以年为单位，除以365

library(dplyr)
###filter函数在dplyr包中，取生存期小于等于3年的数据
rt1<-filter(rt,futime<=3)


rt1$group<-ifelse(rt1$StromalScore<median(rt1$StromalScore),"Low","High")
sur<-Surv(time = rt1$futime,event = rt1$fustat)
mod<-survfit(sur~group,data=rt1)
survdiff(sur~group,data=rt1)

rt1$group2<-ifelse(rt1$ImmuneScore<median(rt1$ImmuneScore),"Low","High")
sur<-Surv(time = rt1$futime,event = rt1$fustat)
mod<-survfit(sur~group2,data=rt1)
survdiff(sur~group2,data=rt1)
