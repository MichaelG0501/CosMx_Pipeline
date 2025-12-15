library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)           # needed for dgCMatrix / CsparseMatrix
library(ggplot2)
library(patchwork)
library(parallel)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/CosMx_Pipeline/cosmx_outs")

n_clusters = 8

meta <- as.data.frame(data.table::fread("/rds/general/project/spatialtranscriptomics/live/CosMx/Bidcell/new/cosmx_metadata_final.csv"))
rownames(meta) <- meta$V1
meta$V1 <- NULL
rownames(meta) <- sapply(strsplit(rownames(meta), "-"), `[`, 1)
meta <- meta[meta$QC == FALSE, ]
meta$coreid <- paste0(meta$ID, "_", meta$position)

data <- data.table::fread(paste0("/rds/general/project/spatialtranscriptomics/live/CosMx/Bidcell/new/expr_mtx_TMA", 1, ".csv"))
genes <- data[[1]]
data <- as.matrix(data[, -1])
rownames(data) <- genes

data2 <- data.table::fread(paste0("/rds/general/project/spatialtranscriptomics/live/CosMx/Bidcell/new/expr_mtx_TMA", 2, ".csv"))
genes2 <- data2[[1]]
data2 <- as.matrix(data2[, -1])
rownames(data2) <- genes2

data <- cbind(data, data2)
rm(data2)
gc()

meta <- meta[meta$n_genes_by_counts > 100, ]
data <- data[, colnames(data) %in% rownames(meta)]
meta <- meta[rownames(meta) %in% colnames(data), ]

cosmx <- CreateSeuratObject(data, meta.data = meta)
rm(data)
gc()

CPM <- apply(cosmx@assays$RNA$counts, 2, function(x) (x / sum(x)) * 1e6)
CPM <- as(CPM, "dgCMatrix")
cosmx@assays$RNA$CPM <- CPM
expr <- log2((CPM / 10) + 1)
expr <- as(expr, "CsparseMatrix")
cosmx@assays$RNA$data <- expr

cosmx_list <- SplitObject(cosmx, split.by = "coreid")
print(names(cosmx_list))
names(cosmx_list) <- unique(cosmx@meta.data$coreid)
print(names(cosmx_list))
rm(cosmx)
gc()

out_dir <- "by_samples"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
mclapply(names(cosmx_list), function(nm) {
  sample_dir <- file.path(out_dir, nm)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
  saveRDS(cosmx_list[[nm]], file.path(sample_dir, paste0(nm, ".rds")), compress = FALSE)
  NULL
}, mc.cores = n_clusters, mc.preschedule = FALSE)
