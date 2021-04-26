library(Seurat)


ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

function_path <- "scripts/functions.R"

source(function_path)
set.seed(0)


data_dir <- getwd()
#input_dir <- paste0(data_dir, "/output/allSamples/")
#input_dir <- paste0(data_dir, "/output/allSamples_GerminalBCells/")
#input_dir <- paste0(data_dir, "/output/allSamples_BCells/")

#seurat_object <- readRDS(paste0(input_dir, "files/all_samples_dubRm_seurat.rda"))
#seurat_object <- readRDS(paste0(input_dir, "files/allSamples_GerminalBCells_seurat.rda"))
#seurat_object <- readRDS(paste0(input_dir, "files/allSamples_BCells_seurat.rda"))

# cluster_file <- paste0(input_dir, "files/cluster_info.txt")
# umap_file <- paste0(input_dir, "files/UMAP_info.txt")
# harmony_file <- paste0(input_dir, "files/harmony_info.txt")
# pca_file <- paste0(input_dir, "files/pca_info.txt")

seurat_object <- readRDS(snakemake@input[[1]])
cluster_file <- snakemake@output[[1]]
cell_file <- snakemake@output[[2]]
pca_file <- snakemake@output[[3]]
umap_file <- snakemake@output[[4]]
harmony_file <- snakemake@output[[5]]


cell_type_info <- seurat_object[["cell_type"]]
#cluster_info <- seurat_object[["new_cluster"]]
cluster_info <- Idents(seurat_object)
UMAP_info <- Embeddings(seurat_object, reduction = "umap")
harmony_info <- Embeddings(seurat_object, reduction = "harmony")
pca_info <- Embeddings(seurat_object, reduction = "pca")

write.table(cell_type_info, file = cell_file, sep = ",",
	row.names = TRUE, col.names = TRUE)

write.table(cluster_info, file = cluster_file, sep = ",",
	row.names = TRUE, col.names = TRUE)

write.table(UMAP_info, file = umap_file, sep = ",",
	row.names = TRUE, col.names = TRUE)

write.table(harmony_info, file = harmony_file, sep = ",",
	row.names = TRUE, col.names = TRUE)

write.table(pca_info, file = pca_file, sep = ",",
	row.names = TRUE, col.names = TRUE)
