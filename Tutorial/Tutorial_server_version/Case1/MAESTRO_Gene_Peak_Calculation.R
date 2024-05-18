library(reticulate)
library(Matrix)
library(MAESTRO)

# Load RDS files -- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201402
retina_atac <- readRDS('GSM6062731_retina_10x.atac.RDS')
multiomics_data <- readRDS('GSM6062732_multiomics_mouse_wt_GeneExpression_Peaks_lists.RDS')
retina_rna <- readRDS('GSM6062732_retina_10x.rna.RDS')

# Extract expression data
atac_counts <- retina_atac@assays$ATAC@counts
rna_counts <- retina_rna@assays$RNA@counts

# Extract gene, cell, and peak names
gene_names_rna <- rownames(rna_counts)
cell_names_rna <- colnames(rna_counts)
cell_names_atac <- colnames(atac_counts)
peak_names <- rownames(atac_counts)

# Save ATAC-seq data as MTX format
writeMM(obj = atac_counts, file = "retina_atac_counts.mtx")

# Save RNA-seq data as MTX format
writeMM(obj = rna_counts, file = "retina_rna_counts.mtx")

# Save gene, cell, and peak names to specified files
write.table(gene_names_rna, file = 'Gene_names.txt', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cell_names_rna, file = 'Cell_names.txt', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cell_names_atac, file = 'ATAC_Cell_names.txt', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(peak_names, file = 'Peak_names.txt', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)


rna <- Matrix::readMM('Gene_Cell.mtx')
atac <- Matrix::readMM('Peak_Cell.mtx')
gene_names <- read.table('Gene_names.txt',sep =',',header = F)
cell_names <- read.table('Cell_names.txt',sep =',',header = F)
atac_cell_names <- read.table('Cell_names.txt',sep =',',header = F)
peak_names <- read.table('Peak_names.txt',sep =',',header = F)
colnames(rna) <- cell_names$V1
rownames(rna) <- gene_names$V1
colnames(atac) <- atac_cell_names$V1
rownames(atac) <- peak_names$V1
rna <- rna[,intersect(colnames(rna),colnames(atac))]
atac <-atac[,intersect(colnames(rna),colnames(atac))]

organism = 'GRCm38' 
gene_count_matrix = rna
peak_count_matrix = atac
print(dim(gene_count_matrix))
print(dim(peak_count_matrix))
gene_peak_list <- list()

for (i in 1:length(seq(1,dim(peak_count_matrix)[1],10000))){
  print(i)
  pbmc_peak <- peak_count_matrix[seq(1,dim(peak_count_matrix)[1],10000)[i]:min(seq(1,dim(peak_count_matrix)[1],10000)[i]+9999,dim(peak_count_matrix)[1]),]
  n <- nrow(pbmc_peak)
  print(n)
  dia <- diag(n)
  dia <- as(dia, "dgCMatrix")
  rownames(dia) <- rownames(pbmc_peak)
  colnames(dia) <- 1:ncol(dia)
  gene_peak <-
    ATACCalculateGenescore(dia,
                           organism = organism,
                           decaydistance = 10000,
                           model = "Enhanced")
  gene_peak_list[i] <- gene_peak
  
}

Gene_Peak <- c()
for (m in gene_peak_list){
  Gene_Peak <- cbind(Gene_Peak,m)
}
colnames(Gene_Peak) <- rownames(peak_count_matrix)
Gene_names <- intersect(rownames(Gene_Peak),rownames(gene_count_matrix))
Gene_Peak <- Gene_Peak[Gene_names,]
gene_count_matrix <- gene_count_matrix[Gene_names,]

write(x=Gene_Peak@Dimnames[[1]],file = 'Gene_names.tsv')
write(x=Gene_Peak@Dimnames[[2]],file = 'Peak_names.tsv')
writeMM(obj = Gene_Peak ,file = 'Gene_Peak.mtx')
write(x=colnames(gene_count_matrix),file = 'Cell_names.tsv')      
writeMM(obj = gene_count_matrix ,file = 'Gene_Cell.mtx')
writeMM(obj = peak_count_matrix ,file = 'Peak_Cell.mtx')
