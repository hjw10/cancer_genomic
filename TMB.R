library(maftools)
laml.maf='./TCGA_maf_vaf.gz'
laml = read.maf(maf = laml.maf)
#data@data$VAF <- data@data$t_alt_count/data@data$t_depth
#k=data.frame(data$Tumor_Sample_Barcode,data$het)
summaryT=data.frame(getSampleSummary(laml))
getGeneSummary(laml)
getClinicalData(laml)
getFields(laml)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw=F)

###onco
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
#查看变异类型对应的颜色
print(vc_cols)

oncoplot(maf = laml, colors = vc_cols, top = 20)
oncoplot(maf = laml, colors = vc_cols,genes = c('DNMT3A','NPM1', 'RUNX1','TTN','CTNNB1','CTNNB2','CTNNB3'))
######
#titv 函数将SNP分类为 Transitions and Transversions ，并以各种方式返回汇总表的列表。汇总的数据还可以可视化为显示六个不同转换的总体分布的boxplot图，以及显示每个样本中的转换分数的堆叠条形图。
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
#任意组基因
oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1','TTN','CTNNB1','CTNNB2','CTNNB3'),colors=vc_cols)

###瀑布图
rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.6)
##单基因点位点突变
lollipopPlot(maf = laml, gene = 'TP53', showMutationRate = TRUE,AACol = 'HGVSp_Short')
#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))



plotVaf(maf = laml, vafCol = 'vaf')
#关键TMB计算
TMB_table=tmb(maf=laml, captureSize = 45, logScale = TRUE)#这个地方50是他们给的，一般是45m
