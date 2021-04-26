library(Seurat)
library(dplyr)
library(pheatmap)
library(viridis)
library(sctransform)
library(DoubletFinder)

data_input <- snakemake@input
function_path <- snakemake@params[[1]]
data_list <- snakemake@params[[2]]
save_path <- snakemake@params[[3]]

# data_input <- c("Tonsil1a_count/outs", "Tonsil1b_count/outs",
#               "Tonsil2a_count/outs", "Tonsil2b_count/outs",
#               "Tonsil3a_count/outs", "Tonsil3b_count/outs")
# function_path <- "scripts/functions.R"
# data_list <- list(Tonsil1 = c("Tonsil1a", "Tonsil1b"),
#                   Tonsil2 = c("Tonsil2a", "Tonsil2b"),
#                   Tonsil3 = c("Tonsil3a", "Tonsil3b"))

print(data_list)
tonsil_list <- names(data_list)

sample_list <- unlist(data_list)

names(data_input) <- sample_list
sample_list <- setNames(sample_list, sample_list)
tonsil_list <- setNames(tonsil_list, tonsil_list)

source(function_path)
set.seed(0)


data_dir <- getwd()
#save_dir <- paste0(data_dir, "/output/mergedTonsil/")
save_dir <- paste0(data_dir, save_path)

lapply(tonsil_list, function(x){
  print(data_list[[x]][1])
  print(data_list[[x]][2])
  })

# tonsil_list <- setNames(tonsil_list)

# Set arguments for normalization and finding variable genes with Seurat
min_cells <- 3
min_features <- 200
project_name <- "scRNA_tonsil"
selection_method <- "vst"
nfeatures <- 2000
### Here I set percent mito quite high to deal with the percent of reads that are
# mitochondrial being high
percent_mito <- 20
ngene <- c(200, 7500)
nADT <- 4000
initial_analysis <- FALSE

if(initial_analysis){
  seurat_list_single <- lapply(sample_list,
    function(x) create_seurat_object_snakemake(sample = x,
      count_path = data_input[[x]], ADT = TRUE, min_features = min_features,
      min_cells = min_cells))


  seurat_list <- lapply(tonsil_list, function(x)
    initial_processing_merge_two(sample1 = data_list[[x]][1],
      sample2 = data_list[[x]][2], seurat_list = seurat_list_single,
      project_name = project_name, merge_name = x,
      selection_method = selection_method, nfeatures = nfeatures,
      percent_mito = percent_mito, ngene = ngene,
      nADT = nADT, sctransform = TRUE, ADT = TRUE))


  seurat_list <- lapply(tonsil_list, function(x)
    PCA_dimRed(sample_object = seurat_list[[x]],
      sample_name = x, save_dir = save_dir))

  saveRDS(seurat_list, file = paste0(save_dir, "files/PCA_seurat.rda"))

  ndims <- c(Tonsil1 = 11, Tonsil2 = 13, Tonsil3 = 12)

  seurat_list <- lapply(tonsil_list,
    function(x) group_cells(sample_object = seurat_list[[x]],
      sample_name = x, nPCs = ndims[[x]], plot_ADT = TRUE,
      save_dir = save_dir))

  seurat_list <- lapply(tonsil_list, function(x) doublet_finder(seurat_list[[x]],
    x, nPCs = ndims[[x]], sct = TRUE, expted_doublets = 0.039))

  lapply(tonsil_list, function(x) {
    saveRDS(seurat_list[[x]], file = paste0(save_dir, "files/",
      x, "_seurat.rda"))
  })

} else {
  seurat_list <- lapply(tonsil_list, function(x){
    return(readRDS(file = paste0(save_dir,
      "files/", x, "_seurat.rda")))
    })
  ndims <- c(Tonsil1 = 11, Tonsil2 = 13, Tonsil3 = 12)
}



plot_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_ADT",
  "nFeature_ADT", "MKI67", "TOP2A")

# # pdf(paste0(save_dir, "images/UMAP_all_samples.pdf"))

# # plotDimRed(seurat_object, col_by = "cluster")

# # plotDimRed(seurat_object, col_by = "orig.ident")

# # ADT_markers <- rownames(seurat_object[["ADT"]])

# # lapply(ADT_markers, function(x) plotDimRed(seurat_object, col_by = x))

# # dev.off()


save_plot <- lapply(tonsil_list, function(x) {
  paste0(save_dir, "/images/umap_", x, "_qc.pdf")
  })

lapply(tonsil_list, function(x) plotDimRed(sample_object = seurat_list[[x]],
  save_plot = save_plot[[x]], col_by = plot_features, return_plot = FALSE))


seurat_cell_markers <- lapply(tonsil_list, function(x) cell_markers(seurat_list[[x]], x,
  save_dir = save_dir))

# seurat_list <- lapply(tonsil_list, function(x) significant_markers(seurat_list[[x]],
#   sig_logFC = 1))

# lapply(tonsil_list, function(x) {
#     saveRDS(seurat_list[[x]], file = paste0(save_dir, "files/",
#       x, "_seurat.rda"))
#     })


# Remove dublet cells
seurat_list_dedub <- lapply(tonsil_list, function(x)
  remove_doublets(seurat_list[[x]]))

# Repeat DE with dedub cells
ADT_genes <- rownames(GetAssayData(object = seurat_list_dedub[[1]], slot = "data",
      assay = "ADT"))

seurat_list_dedub <- lapply(tonsil_list, function(x){
  seurat_object <- seurat_list[[x]]
  seurat_object@misc$DE <- NULL
  return(seurat_object)
  })

plot_features_new <- c("cluster", "orig.ident", ADT_genes, plot_features)

save_plot_dedub <- lapply(tonsil_list, function(x) {
  paste0(save_dir, "/images/umap_", x, "_dubRm.pdf")
  })

lapply(tonsil_list, function(x) plotDimRed(sample_object = seurat_list_dedub[[x]],
  save_plot = save_plot_dedub[[x]], col_by = plot_features_new, return_plot = FALSE))


dedub_cell_markers <- lapply(tonsil_list, function(x) cell_markers(seurat_list_dedub[[x]], x,
  save_dir = save_dir))

# seurat_list_dedub <- lapply(tonsil_list, function(x) significant_markers(seurat_list_dedub[[x]],
#   sig_logFC = 1))

lapply(tonsil_list, function(x) {
    saveRDS(seurat_list_dedub[[x]], file = paste0(save_dir, "files/",
      x, "_dubRm_seurat.rda"))
    })


