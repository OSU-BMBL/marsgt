library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
inputFolder = '/fs/ess/PCON0022/DMT/GSE201402/BC_Pathway/GSEA_Cluster10/'
df = read.table(paste(inputFolder,"gene_score.txt",sep=''),col.names=c('SYMBOL','logFC')) #读入txt

df_id<-bitr(df$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Mm.eg.db")#对应的物种，小鼠的是org.Mm.eg.db

df_all<-merge(df,df_id,by="SYMBOL",all=F)#使用merge合并
# df_all <- df_all[1:3000,]
head(df_all) #再看看数据
dim(df_all) #因为有一部分没转换成功，所以数量就少了。

df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$logFC #把foldchange按照从大到小提取出来
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(gene_fc)

GO <- gseGO(
  gene_fc, #gene_fc
  ont = "ALL",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Mm.eg.db,#小鼠注释基因
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",#p值校正方法
)

sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
head(sortGO)
dim(sortGO)

gseaplot2_modified <- function(GO, id, title, color = "black", base_size = 18,
                               ES_geom = "line", subplots = c(1, 2), rel_heights = c(1.5, .1, .5)) {
  plot <- gseaplot(GO, id, title, color, base_size, by = "runningScore")
  pvalue <- sortGO[sortGO$ID == id, "pvalue"]
  qvalue <- sortGO[sortGO$ID == id, "qvalue"]
  plot +
    ggplot2::annotate("text", x = Inf, y = Inf,
                      label = paste("p-value =", round(pvalue, 4), "\nq-value =", round(qvalue, 4)),
                      hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
}

p <- gseaplot2_modified(GO, sortGO$ID[1], title = 'Structural Constituent of Eye Lens',
                   color = "black", base_size = 18, subplots = c(1, 2),
                   rel_heights = c(1.5, .1, .5))
svg (file=paste0("/users/PCON0022/duanmaoteng/Final_HGT/GSE201402/DeepMars-main/figures/Structural_Constituent_of_Eye_Lens.svg"),width = 10, height = 6)
print(p)
dev.off()
