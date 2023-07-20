# BiocManager::install("slingshot")
library(Seurat)
library(ggplot2)
library(slingshot)
library(RColorBrewer)
library(viridis)

rna <- Matrix::readMM('B_Gene_Cell.mtx')
gene_names <- read.csv("Gene_names.tsv",header =FALSE)
cell_names <- read.table('Cell_names.txt',sep =',',header = F)
cell_type <- read.table('Cell_types.txt',sep =',',header = F)
# cell_type <- lapply(cell_type, as.character)
cell_type$V1 <- as.character(cell_type$V1)
B_umap <- as.matrix(read.table('fa_B_umap.txt',sep =',',header = F))
# B_umap <- as.matrix(read.table('B_umap.txt',sep =',',header = F))

rna <- as(rna, "dgCMatrix")
rna <- t(rna)
colnames(rna) <- cell_names$V1
rownames(rna) <- gene_names$V1
# rownames(cell_type) <- cell_names$V1
# rownames(gene_names) <- gene_names$V1
# names(gene_names)<-c('gene_short_name')

B_cell <- CreateSeuratObject(counts = rna)
B_cell@meta.data$labels <- as.factor(cell_type$V1)
colnames(B_umap) <- c('UMAP1','UMAP2')
rownames(B_umap) <- cell_names$V1

# Define a new S4 class, this class has a cell.embeddings slot
setClass(
  "UMAP",
  representation(
    cell.embeddings = "matrix"
  )
)

# Create an instance of this class
umap_instance <- new("UMAP")
# Assign a value to the cell.embeddings slot
umap_instance@cell.embeddings <- B_umap
# Assign this instance to B_cell@reductions$wnn.umap
B_cell@reductions$wnn.umap <- umap_instance

# B_cell@reductions$wnn.umap@cell.embeddings <- B_umap
colnames(B_umap) <- paste0("UMAP_", 1:2)
B_cell[["umap"]] <- CreateDimReducObject(embeddings = B_umap, key = "UMAP_", assay = DefaultAssay(B_cell))

# slingClusterLabels(sc)
sc <- as.SingleCellExperiment(B_cell)
sc <- slingshot(sc, clusterLabels = "labels", 
                reducedDim = "UMAP", 
                start.clus = '10',
                end.clus = c('6'),
                extend = 'pc1',
)

# slingPseudotime(sc)

ip_file = 'plot_save/'
svg (file=paste(ip_file,"slingshot_c.svg",sep='/'),width=12, height=8)
plot(reducedDims(sc)$UMAP,
     col = c("#FFACAA","#7FD2FF","#FF9D1E","#894FC6")[as.factor(sc$labels)],pch=16,
     asp=1)
lines(SlingshotDataSet(sc), lwd=2,col = 'Black',type = 'c')
dev.off()

svg (file=paste(ip_file,"slingshot_psedotime.svg",sep='/'),width=12, height=8)
nc <- 3
pt <- slingPseudotime(sc)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sc), col = colors, pch = 16, cex = 0.5, main = i)
  lines(SlingshotDataSet(sc), lwd=2,col = 'Black',type = 'c')
}
