
library(MCPcounter)
library(pheatmap)
genes <- data.table::fread("genes.txt",data.table = F)

probesets <- data.table::fread("probesets.txt",data.table = F,header = F)
exprMat <- read.table("exprMat.txt", header=T, sep="\t", check.names=F,row.names = 1)
### 成功！！！
results<- MCPcounter.estimate(exprMat,
                              featuresType= "HUGO_symbols",
                              probesets=probesets,
                              genes=genes)

heatmap(as.matrix(results),col=colorRampPalette(c("blue","white","red"))(100)) 
