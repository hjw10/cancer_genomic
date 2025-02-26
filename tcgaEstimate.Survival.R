
#install.packages("survival")

setwd("C:\\Users\\lexb4\\Desktop\\tcgaEstimate\\09.survival")    #工作目录（需修改）

library(survival)
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
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

  pdf(file=paste(gene,".survival.pdf",sep=""),
      width=6,
      height=6)
  plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste(gene,"(p=", pValue ,")",sep="") )
  legend("topright", 
       c("High","Low"), 
       lwd=2, 
       col=c("red","blue"))
  dev.off()
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)
