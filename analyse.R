# setwd("..")
#counts = read.csv(".\\LUAD_Counts_expMatrix.csv", header = T,row.names = 1)
#LUAD_Counts_expMatrix = counts
library(TCGAbiolinks)
if (0){
  quary = GDCquery(project = "TCGA-LUAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTseq - Counts")
  
  samplesDown = getResults(query,cols = c("cases"))
}

samplesDown = names(LUAD_Counts_expMatrix)
dataSmTP = TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
dataSmNT = TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")

Counts = data.frame(c(LUAD_Counts_expMatrix[,dataSmNT],LUAD_Counts_expMatrix[,dataSmTP]))
rownames(Counts) = row.names(LUAD_Counts_expMatrix)
colnames(Counts) = c(dataSmNT, dataSmTP)

library("edgeR")
group = c(rep(1,59),rep(2,533))
y = DGEList(counts=Counts, group = group)

keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=FALSE]

y = calcNormFactors(y)

y = estimateDisp(y)

et = exactTest(y)

et = topTags(et, n=10000000)

et = as.data.frame(et)
et = cbind(rownames(et),et)

colnames(et) = c("gene_id","log2FoldChange","log2CPM","PValue","FDR")
write.table(et,"./all_LUAD_DEG.xls", sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE,na="")

etSig = et[which(et$PValue<0.5 & abs(et$log2FoldChange)>1),]

etSig[which(etSig$log2FoldChange >0),"up_down"] = "UP"
etSig[which(etSig$log2FoldChange <0),"up_down"] = "DOWN"

write.table(etSig, "./LUAD_DEG.xls",sep = "\t",col.names = TRUE,row.names = F,quote = F,na="")

TCGAVisualize_volcano(etSig$log2FoldChange, etSig$PValue,
                      filename = "volcano.pdf", xlab = "logFC",
                      names = rownames(etSig), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = rownames(etSig)[which(abs(etSig$log2FoldChange) >= 100)],
                      highlight.color = "orange")

library(AnnotationHub)	#library导入需要使用的数据包
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(dplyr)
library(ggplot2)


EG2Ensembl=toTable(org.Hs.egENSEMBL)	 #将ENTREZID和ENSEMBL对应的数据存入该变量
tmp = merge(Ensembl_ID_To_Genename,etSig,by="gene_id")
f=tmp$Ensembl_ID	#list转化为字符向量
geneLists=data.frame(ensembl_id=f)
results=merge(geneLists,EG2Ensembl,by='ensembl_id',all.x=T)
id=na.omit(results$gene_id)  #提取出非NA的ENTREZID

ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "MF", pvalueCutoff = 0.01, readable= TRUE) #GO富集分析
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图
plotGOgraph(ego) 	#GO图，看不清楚可以尝试左上角另存为pdf

#KEGG分析#
ekk <- enrichKEGG(gene= id,organism  = 'hsa', qvalueCutoff = 0.05)	 #KEGG富集分析
dotplot(ekk,font.size=8)	# 画气泡图


library(pheatmap)
pheatmap(head(LUAD_Counts_expMatrix,10), #表达数据
         cluster_rows = T,#行聚类
         cluster_cols = T,#列聚类
         
         annotation_legend=TRUE, # 显示样本分类
         show_rownames = T,# 显示行名
         show_colnames = T,# 显示列名
         scale = "row", #对行标准化
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100) # 热图基准颜色
)
