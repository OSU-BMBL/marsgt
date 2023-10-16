library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
inputFolder = '/fs/ess/PCON0022/DMT/GSE201402/BC_Pathway/GSEA_Cluster10/'
df = read.table(paste(inputFolder,"gene_score.txt",sep=''),col.names=c('SYMBOL','logFC')) # Read in txt

df_id<-bitr(df$SYMBOL, # The column to be converted is the SYMBOL column of df dataframe
            fromType = "SYMBOL", # ID type to be converted
            toType = "ENTREZID", # ID type to convert to
            OrgDb = "org.Mm.eg.db") # Corresponding species, mouse is org.Mm.eg.db

df_all<-merge(df,df_id,by="SYMBOL",all=F) # Merge using merge
head(df_all) # Check the data again
dim(df_all) # There are some that didn't convert successfully, so the count is reduced.

df_all_sort <- df_all[order(df_all$logFC, decreasing = T),] # Sort by logFC in descending order
gene_fc = df_all_sort$logFC # Extract foldchange from largest to smallest
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID # Add corresponding ENTREZID to the extracted foldchange
head(gene_fc)

GO <- gseGO(
  gene_fc, 
  ont = "ALL",# "BP", "MF", "CC" or "ALL"
  OrgDb = org.Mm.eg.db, # Mouse annotation gene
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", # p-value adjustment method
)

sortGO<-GO[order(GO$enrichmentScore, decreasing = T),] # Sort by enrichment score from high to low
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
