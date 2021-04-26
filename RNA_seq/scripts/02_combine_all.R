library(Seurat)
library(dplyr)
library(pheatmap)
library(viridis)
library(harmony)

print("combine_all")
data_input <- snakemake@input
function_path <- snakemake@params[[1]]
tonsil_list <- snakemake@params[[2]]
output_file <- snakemake@output[[1]]
save_dir <- snakemake@params[[3]]

# data_input <- c("output/mergedTonsil_snakemake/files/Tonsil1_dubRm_seurat.rda",
#                 "output/mergedTonsil_snakemake/files/Tonsil2_dubRm_seurat.rda",
#                 "output/mergedTonsil_snakemake/files/Tonsil3_dubRm_seurat.rda")
# function_path <- "scripts/functions.R"
# tonsil_list <- c(Tonsil1 = c("Tonsil1a", "Tonsil1b"),
#                  Tonsil2 = c("Tonsil2a", "Tonsil2b"),
#                  Tonsil3 = c("Tonsil3a", "Tonsil3b"))
# output_file <- "output/allSamples_snakemake/files/allSamples_dubRm_seurat.rda"
# save_dir <- "output/allSamples_snakemake/"

print(tonsil_list)
tonsil_list <- names(tonsil_list)
names(data_input) <- tonsil_list
print(data_input)

tonsil_list <- setNames(tonsil_list, tonsil_list)

source(function_path)
set.seed(0)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))


# Set arguments for normalization and finding variable genes with Seurat
selection_method <- "vst"


first_analysis <- TRUE
sctransform <- TRUE
harmony_combine <- TRUE

# snakemake mt calls
#vars_to_regress <- "percent.mt"
# resolution <- 0.9
# ndims <- 27
# seed <- 0
# nfeatures <- 2000

# Snakemake nomt calls
vars_to_regress <- NULL
resolution <- 0.8
ndims <- 27
seed <- 0
nfeatures <- 3000
normalized <- TRUE

if(first_analysis){

  seurat_list_dubrm <- lapply(data_input, function(x){
    readRDS(file = x)
    })

  print(seurat_list_dubrm)

  if (harmony_combine){

    if (!normalized){

      seurat_object_dubrm <- merge_objects(seurat_list = seurat_list_dubrm,
        sample_name = "all_samples_dubRm",
        nfeatures = nfeatures, vars_to_regress = vars_to_regress,
        selection_method = selection_method)
      saveRDS(seurat_object_dubrm, file = paste0(save_dir, "files/all_samples_sct.rda"))
    } else {
      seurat_object_dubrm <- readRDS(file = paste0(save_dir, "files/all_samples_sct.rda"))
    }

    # Remove the features that are skewing results but have little
    # biological significance
    variable_features_dubRm <- VariableFeatures(seurat_object_dubrm)
    to_remove <- c(variable_features_dubRm[grepl("IGKC", variable_features_dubRm)],
                   variable_features_dubRm[grepl("IGLC", variable_features_dubRm)],
                   variable_features_dubRm[grepl("IGLV", variable_features_dubRm)],
                   variable_features_dubRm[grepl("HLA", variable_features_dubRm)],
                   variable_features_dubRm[grepl("IGH", variable_features_dubRm)])
    variable_features_dubRm <- variable_features_dubRm[!(variable_features_dubRm %in% 
                                                       to_remove)]

    set.seed(0)
    seurat_object_dubrm <- run_harmony(seurat_object_merge = seurat_object_dubrm,
      save_dir = save_dir, sample_name = "all_samples_dubRm", resolution = resolution,
      plot_ADT = TRUE, ndims = ndims, variable_features = variable_features_dubRm,
      seed = seed)
    # pdf(paste0(save_dir, "images/umap_test3.pdf"))
    # lapply(0:10, function(x){
    #   print(x)
    #   set.seed(0)
    #   seurat_object_dubrm2 <- run_harmony(seurat_object_merge = seurat_object_dubrm,
    #     save_dir = save_dir, sample_name = "all_samples_dubRm", resolution = 0.9,
    #     plot_ADT = TRUE, ndims = 27, variable_features = variable_features_dubRm,
    #     seed = x)
    #   p <- plotDimRed(seurat_object_dubrm2, col_by = "cluster")
    #   print(p[[1]] + ggplot2::ggtitle(paste0("seed_", x)))
    #   })

    # dev.off()

  } else{


    seurat_object_dubrm <- integrate_samples(object_list = seurat_list_dubrm,
      sctransform = sctransform)

    seurat_object_dubrm <- PCA_dimRed(sample_object = seurat_object_dubrm,
      sample_name = "all_samples_dubRm", save_dir = save_dir)

    ndims <- 13

    seurat_object_dubrm <- group_cells(sample_object = seurat_object_dubrm,
      sample_name = "all_samples_dubRm", nPCs = ndims, plot_ADT = TRUE,
      save_dir = save_dir, resolution = 0.9)
  }

  saveRDS(seurat_object_dubrm, file = paste0(save_dir, "files/all_samples_tmp.rda"))
} else {
  seurat_object_dubrm <- readRDS(file = paste0(save_dir, "files/all_samples_tmp.rda"))
}

plot_features <- c("cluster", "nFeature_RNA", "nCount_RNA", "percent.mt",
  "nCount_ADT", "nFeature_ADT", "MKI67", "TOP2A")

# Plot qc metrics
plotDimRed(sample_object = seurat_object_dubrm, col_by = plot_features,
  save_plot = paste0(save_dir, "images/qc_all_samples_dubRm.pdf"))


# Plot UMAPs showing how well clusters overlap
pdf(paste0(save_dir, "images/cluster_by_sample_dubrm.pdf"))

seurat_object_dubrm$Tonsil <- gsub("a|b", "", seurat_object_dubrm$orig.ident)

# Umap coloring by each tonsil
lapply(tonsil_list, function(x) full_umap(seurat_object_dubrm, data_set = x,
  col_by = "cluster", meta_data_col = "Tonsil"))

# Bar plots of percentages of each cluster for each tonsil
stage_list_all_dubrm <- lapply(tonsil_list, function(x)
                         populations_dfs_new(seurat_object_dubrm,
                         x, subsample = TRUE, subsample_by = "Tonsil",
                         meta_data_col = "seurat_clusters"))
stage_list_all_dubrm <- do.call("rbind", stage_list_all_dubrm)

stage_list_all_dubrm$sample <- factor(stage_list_all_dubrm$sample)

population_plots(stage_list_all_dubrm)

dev.off()

# Markers of each cluster
new_list_dubrm <- cell_markers(seurat_object_dubrm,
  sample_name = "all_samples_dubRm", save_dir = save_dir,
  sctransform = FALSE)


# seurat_object_dubrm <- significant_markers(seurat_object_dubrm, sig_logFC = 1,
#   sctransform = FALSE)

saveRDS(seurat_object_dubrm, file = output_file)


# new_list_dubrm <- print_DE(seurat_object_dubrm,
#   output_file = paste0(save_dir, "files/", x, "_DE_genes_dubRm.txt")))

# cluster_markers_dubrm <- cluster_genes(DE_list = new_list_dubrm,
#   seurat_object = seurat_object_dubrm,
#   output_file = paste0(save_dir, "files/", x, "_DE_cluster_markers_dubRm.txt")))