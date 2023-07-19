library(Seurat)
library(ggplot2)
library(UCell)
library(Hmisc)
library(tidyr)
library(dplyr)
library(ggpubr)

rds <- readRDS("pbmc_multimodal.downsampled20k.Tcell.seurat.RNA.rds")
rna <- Matrix::readMM('B_Gene_Cell.mtx')
gene_names <- read.csv("Gene_names.tsv",header =FALSE)
cell_names <- read.table('Cell_names.txt',sep =',',header = F)
cell_type <- read.table('Cell_types.txt',sep =',',header = F)
# cell_type <- lapply(cell_type, as.character)
cell_type$V1 <- as.character(cell_type$V1)
B_umap <- as.matrix(read.table('fa_B_umap.txt',sep =',',header = F))
pathway <- read.csv("/fs/ess/PCON0022/LiJingxian/Data/pathway/KEGG/human/pathway_human.csv")

Cellular_Processes <- unique(pathway$pathway_function[pathway$biofun_name=="09140 Cellular Processes"])
Human_Diseases <- unique(pathway$pathway_function[pathway$biofun_name=="09160 Human Diseases"])

# pathway_combine <- 
rna <- as(rna, "dgCMatrix")
rna <- t(rna)
colnames(rna) <- cell_names$V1
rownames(rna) <- gene_names$V1
rownames(cell_type) <- cell_names$V1
rownames(gene_names) <- gene_names$V1
names(gene_names)<-c('gene_short_name')

B_cell <- CreateSeuratObject(counts = rna)
B_cell$cell_type <- cell_type$V1
B_cell<- NormalizeData(B_cell, normalization.method = "LogNormalize", scale.factor = 10000)
B_cell <- FindVariableFeatures(B_cell, selection.method = "vst", nfeatures = 2000)
B_cell <- ScaleData(B_cell,features=VariableFeatures(B_cell)) 
B_cell <- RunPCA(B_cell, features = VariableFeatures(object = B_cell))
B_cell@meta.data
B_cell <- FindNeighbors(B_cell, dims = 1:10)
B_cell <- FindClusters(B_cell, resolution = 0.5)
B_cell <- RunUMAP(B_cell, dims = 1:10)

colnames(B_umap) <- c('UMAP1','UMAP2')
rownames(B_umap) <- cell_names$V1
B_cell@reductions$wnn.umap <- rds@reductions$wnn.umap
B_cell@reductions$wnn.umap@cell.embeddings <- B_umap
colnames(B_umap) <- paste0("UMAP_", 1:2)
B_cell[["umap"]] <- CreateDimReducObject(embeddings = B_umap, key = "UMAP_", assay = DefaultAssay(B_cell))


## PI3K signaling pathway
PI3K_signaling_gene <- c('RPS6', 'SYK', 'EIF4B', 'PKN3', 'TLR2', 'MAP2K1', 'RAC1', 'GRB2',
                         'SOS1', 'PTK2', 'PRKAA1', 'NFKB1', 'PTEN', 'HRAS', 'PIK3CA',
                         'MDM2', 'STK11', 'JAK1', 'MTOR', 'CD19', 'MTCP1', 'BCL2')
print(length(PI3K_signaling_gene))
markers <- list()
PI3K_signaling_gene_list <- c(PI3K_signaling_gene)
pathway_name <- trimws(toString(strsplit("PI3K pathway",'\\[')[[1]][1]))
pathway_name <- gsub("-","_",pathway_name)
pathway_name <- gsub(" ","_",pathway_name)
# pathway_name <- 'Apoptosis_signaling_pathway'
markers$V1 <- c(PI3K_signaling_gene_list)
names(markers) <- pathway_name
signature.names <- paste0(names(markers),"_UCell")
# 主函数
B_cell <- AddModuleScore_UCell(B_cell, features = markers)
B_cell$PI3K_pathway_UCell <- B_cell$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
B_cells_0 <- subset(B_cell, subset = cell_type == "0")
B_cells_6 <- subset(B_cell, subset = cell_type == "6")
B_cells_10 <- subset(B_cell, subset = cell_type == "10")
B_cells_13 <- subset(B_cell, subset = cell_type == "13")
PI3K_B_cells_0_value <- B_cells_0$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
PI3K_B_cells_6_value <- B_cells_6$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
PI3K_B_cells_10_value <- B_cells_10$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
PI3K_B_cells_13_value <- B_cells_13$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)

## MAPK signaling pathway
MAPK_signaling_gene <- c('MAP2K6', 'GNA12', 'IRAK1', 'RASGRF1', 'TGFBR1', 'MAP2K1',
                         'MAP3K4', 'RASA2', 'MAP3K7', 'BRAF', 'MAP3K1', 'MAP4K2', 'TAB2',
                         'RAC1', 'TRAF6', 'RAPGEF2', 'TAOK3', 'SOS1', 'TRAF2', 'MAP3K5',
                         'STK3', 'MAP2K2', 'ARAF', 'RAP1A', 'TGFB1', 'PAK1')
print(length(MAPK_signaling_gene))
markers <- list()
MAPK_signaling_gene_list <- c(MAPK_signaling_gene)
pathway_name <- trimws(toString(strsplit("MAPK pathway",'\\[')[[1]][1]))
pathway_name <- gsub("-","_",pathway_name)
pathway_name <- gsub(" ","_",pathway_name)
# pathway_name <- 'Apoptosis_signaling_pathway'
markers$V1 <- c(MAPK_signaling_gene_list)
names(markers) <- pathway_name
signature.names <- paste0(names(markers),"_UCell")
# 主函数
B_cell <- AddModuleScore_UCell(B_cell, features = markers)
B_cell$MAPK_pathway_UCell <- B_cell$MAPK_pathway_UCell/max(B_cell$MAPK_pathway_UCell)
B_cells_0 <- subset(B_cell, subset = cell_type == "0")
B_cells_6 <- subset(B_cell, subset = cell_type == "6")
B_cells_10 <- subset(B_cell, subset = cell_type == "10")
B_cells_13 <- subset(B_cell, subset = cell_type == "13")
MAPK_B_cells_0_value <- B_cells_0$MAPK_pathway_UCell/max(B_cell$MAPK_pathway_UCell)
MAPK_B_cells_6_value <- B_cells_6$MAPK_pathway_UCell/max(B_cell$MAPK_pathway_UCell)
MAPK_B_cells_10_value <- B_cells_10$MAPK_pathway_UCell/max(B_cell$MAPK_pathway_UCell)
MAPK_B_cells_13_value <- B_cells_13$MAPK_pathway_UCell/max(B_cell$MAPK_pathway_UCell)

## PD signaling pathway
PD_signaling_gene <- c('PIK3CA', 'MAP2K6', 'JAK1', 'MAP3K3', 'LCK', 'PDCD1', 'ZAP70',
                       'MTOR', 'MAPK14', 'MAP2K1', 'NFKB1', 'TRAF6', 'HRAS', 'STAT1',
                       'NFATC1', 'IFNGR1')
print(length(PD_signaling_gene))
markers <- list()
PD_signaling_gene_list <- c(PD_signaling_gene)
pathway_name <- trimws(toString(strsplit("PD pathway",'\\[')[[1]][1]))
pathway_name <- gsub("-","_",pathway_name)
pathway_name <- gsub(" ","_",pathway_name)
# pathway_name <- 'PD_signaling_pathway'
markers$V1 <- c(PD_signaling_gene_list)
names(markers) <- pathway_name
signature.names <- paste0(names(markers),"_UCell")
# 主函数
B_cell <- AddModuleScore_UCell(B_cell, features = markers)
B_cell$PD_pathway_UCell <- B_cell$PD_pathway_UCell/max(B_cell$PD_pathway_UCell)
B_cells_0 <- subset(B_cell, subset = cell_type == "0")
B_cells_6 <- subset(B_cell, subset = cell_type == "6")
B_cells_10 <- subset(B_cell, subset = cell_type == "10")
B_cells_13 <- subset(B_cell, subset = cell_type == "13")
PD_B_cells_0_value <- B_cells_0$PD_pathway_UCell/max(B_cell$PD_pathway_UCell)
PD_B_cells_6_value <- B_cells_6$PD_pathway_UCell/max(B_cell$PD_pathway_UCell)
PD_B_cells_10_value <- B_cells_10$PD_pathway_UCell/max(B_cell$PD_pathway_UCell)
PD_B_cells_13_value <- B_cells_13$PD_pathway_UCell/max(B_cell$PD_pathway_UCell)

B_cells_0_value <- B_cells_0$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
B_cells_6_value <- B_cells_6$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
B_cells_10_value <- B_cells_10$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)
B_cells_13_value <- B_cells_13$PI3K_pathway_UCell/max(B_cell$PI3K_pathway_UCell)

# 创建数据框
data <- list(PI3K_B_cells_10_value, PI3K_B_cells_13_value, PI3K_B_cells_0_value, PI3K_B_cells_6_value,
             MAPK_B_cells_10_value, MAPK_B_cells_13_value, MAPK_B_cells_0_value, MAPK_B_cells_6_value,
             PD_B_cells_10_value, PD_B_cells_13_value, PD_B_cells_0_value, PD_B_cells_6_value)

# 为每个向量创建数据框
data_frames <- lapply(1:length(data), function(i) {
  pathway <- c("PI3K", "MAPK", "PD")[(i - 1) %/% 4 + 1]
  B_cells <- c("B_cells_10", "B_cells_13", "B_cells_0", "B_cells_6")[(i - 1) %% 4 + 1]
  data.frame(value = data[[i]], pathway = pathway, B_cells = B_cells)
})

# 将所有数据框合并
data_long <- do.call(rbind, data_frames)

# 指定 group 变量的顺序
data_long$group <- paste(data_long$pathway, data_long$B_cells, sep = "_")
data_long$group <- factor(data_long$group, levels = unique(data_long$group))

# 修改 position_dodge() 函数的 width 参数以减小组内类别之间的距离
dodge <- position_dodge(width = 0.1)

# 使用 ggplot2 绘制箱线图
p <- ggplot(data_long, aes(x = group, y = value, fill = B_cells)) +
  geom_boxplot(position = dodge, width = 0.5) +
  theme_minimal() +
  labs(title = "Boxplots of Pathways by B_cells Categories", x = "Groups", y = "Values") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("B_cells_0" = "#FFACAA", "B_cells_6" = "#894FC6", "B_cells_10" = '#7FD2FF', "B_cells_13" = "#FF9D1E")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off") +
  # 去掉背景中的网格虚线
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- p +
  stat_compare_means(comparisons = list(
    c("PI3K_B_cells_13", "PI3K_B_cells_0"), c("PI3K_B_cells_13", "PI3K_B_cells_6"), c("PI3K_B_cells_13", "PI3K_B_cells_10"),
    c("MAPK_B_cells_13", "MAPK_B_cells_0"), c("MAPK_B_cells_13", "MAPK_B_cells_6"), c("MAPK_B_cells_13", "MAPK_B_cells_10"),
    c("PD_B_cells_13", "PD_B_cells_0"), c("PD_B_cells_13", "PD_B_cells_6"), c("PD_B_cells_13", "PD_B_cells_10")),
    method = "t.test",label = "p.format", label.y = c(1.1, 1.3, 1.2,1.1, 1.3, 1.2,1.1, 1.3, 1.2,1.1, 1.3, 1.2))
print(p)
