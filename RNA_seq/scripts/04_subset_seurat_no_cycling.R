library(Seurat)
library(dplyr)
library(pheatmap)
library(viridis)
library(harmony)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
function_path <- snakemake@params[[1]]
subset_list <- snakemake@params[[2]]
output_name <- snakemake@params[[3]]
save_dir <- snakemake@params[[4]]
tonsil_list <- snakemake@params[[5]]

seed <- 42

tonsil_list <- names(tonsil_list)

tonsil_list <- setNames(tonsil_list, tonsil_list)

print(subset_list)

if(grepl("_no_cycling", output_name)){
  subset_name <- sub("_no_cycling", "", output_name)
} else {
  subset_name <- output_name
}

subset_cells <- subset_list[[subset_name]]

ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

source(function_path)
set.seed(0)

print("subset list")
print(subset_list)
print("subset cells")
print(subset_name)

print("output name")
print(subset_name)
print(output_name)

data_dir <- getwd()

size_list <- c(BCells = 0.5, GCBcells = 1, ActiveNaiveB = 0.5, TCells = 0.5,
  nonLymphoid = 0.5, NonGCBcells = 0.5)
ndims_list <- c(BCells = 25, GCBcells = 25, ActiveNaiveB = 25, TCells = 20,
  nonLymphoid = 10, NonGCBcells = 25)
resolution_list <- c(BCells = 0.6, GCBcells = 0.6, ActiveNaiveB = 0.6, TCells = 0.6,
  nonLymphoid = 0.4, NonGCBcells = 0.6)

cycling_clusters <- c("10_Dark_GC_B", "11_Dark_GC_B", "12_Dark_GC_B")

# ndims GCB and T cells = 20 for previous run
# T cell resoultion was 0.6 for previous run
# B cell resolution was 0.6 for previous run

# Germinal center looks good. Still working on T cells.

sctransform <- TRUE

vars_to_regress <- NULL

source(function_path)

# print(seed)

seurat_object <- readRDS(input_file)

seurat_object@misc$DE <- NULL

# DefaultAssay(seurat_object) <- "RNA"

# print(seurat_object)

# print("normalize")
# seurat_object <- NormalizeData(seurat_object)

# print(seurat_object)

# DefaultAssay(seurat_object) <- "SCT"


print("subset")
# Changed this line. Cell labels referred to labels decided based on individual tonsils, cell_type referred
# to cell labels decided on all. Cell type was used for the "_new" addition to the output dir
#seurat_subset <- subset(x = seurat_object, subset = cell_labels %in% subset_list)
seurat_subset <- subset(x = seurat_object, subset = cell_type %in% subset_cells)
if(grepl("_no_cycling", output_name)){
  Idents(seurat_subset) <- "num_cell_type" 
  seurat_subset <- subset(x = seurat_subset, idents = cycling_clusters, invert = TRUE)
}
# DefaultAssay(seurat_subset) <- "RNA"

# seurat_subset <- FindVariableFeatures(seurat_subset, nfeatures = 2000)

# variable_features_subset <- VariableFeatures(seurat_subset)

# DefaultAssay(seurat_subset) <- "SCT"

# remove SCT assay
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset[["SCT"]] <- NULL


# # Rerun SCTransform
seurat_subset <- SCTransform(seurat_subset,
      vars.to.regress = vars_to_regress,
      verbose = FALSE)


# Remove the features that are skewing results but have little
# biological significance
variable_features_subset <- VariableFeatures(seurat_subset)
to_remove <- c(variable_features_subset[grepl("IGKC", variable_features_subset)],
               variable_features_subset[grepl("IGLC", variable_features_subset)],
               variable_features_subset[grepl("IGLV", variable_features_subset)],
               variable_features_subset[grepl("HLA", variable_features_subset)],
               variable_features_subset[grepl("IGH", variable_features_subset)])
variable_features_subset <- variable_features_subset[!(variable_features_subset %in% 
                                                   to_remove)]


seurat_subset <- run_harmony(seurat_object_merge = seurat_subset,
  save_dir = save_dir, sample_name = paste0("all_samples_", output_name),
  resolution = resolution_list[[subset_name]],
  plot_ADT = TRUE, ndims = ndims_list[[subset_name]],
  variable_features = variable_features_subset, seed = seed)

# seed = 0 for previous run

# seurat_subset <- harmony_RNA(seurat_object = seurat_subset, save_dir = save_dir,
# 	sample_name = "allSamples", ndims = ndims_list[[output_name]],
#   resolution = resolution_list[[output_name]])

# if(output_name == "GCBcells"){
#   ndims <- 30
#   resolution <- 0.6
#   seurat_subset <- RunPCA(seurat_subset, verbose = FALSE,
#       features = variable_features_subset)

#   seurat_subset <- RunHarmony(seurat_subset, group.by.vars = "orig.ident",
#     assay.use = "SCT")
#   print("UMAP")
#   seurat_subset <- RunUMAP(seurat_subset, reduction = "harmony", dims = 1:ndims,
#     umap.method = "umap-learn", metric = "correlation", seed.use = 12)
#   print("Neighbors")
#   seurat_subset <- FindNeighbors(seurat_subset, reduction = "harmony",
#     dims = 1:ndims)
#   print("clusters")
#   seurat_subset <- FindClusters(seurat_subset, resolution = resolution)

# }

pdf(paste0(save_dir, "images/all_samples_", output_name, "_UMAP.pdf"))

plotDimRed(sample_object = seurat_subset,
	         col_by = c("cluster", "cell_type"), return_plot = TRUE,
           size = size_list[[subset_name]])

lapply(tonsil_list, function(x) full_umap(seurat_subset, data_set = x,
  col_by = "seurat_clusters",
  meta_data_col = "orig.ident", size = size_list[[subset_name]]))

gene_list <- c("POU2F2", "MYC", "CD83", "CD40", "NFKBIA", "EGR1", "EGR2", "EGR3",
  "SLAMF1", "ICAM1", "IL2RG", "BATF", "TCF3", "HMMR",
  "FOXP1", "CXCR4", "PCNA", "TOP2A", "MKI67", "FCRL4",
  "PRDM1", "RORA", "RGS1", "CYTOR", "FOXP3", # T reg
                   "CXCL13", # Tfh CXCL13
                   "CXCR5", "FKBP5", "PDCD1", "BCL6", "PASK", "ST8SIA1", # Tfh
                   "IL7R", "NABP1", "ANK3", "ZFP36L2", # Naive T cells
                   "CCR7", "LEF1", "ITGA6", "KLF2", "SELL", #Tcm
                   "GZMA", "CCL5", "GZMK", "CST7", "CCL4", #CTL
                   "PRF1", "ANXA1", #CTL_PRF1
                   "MKI67", "PCNA", "TOP2A", #Circulating
                   "GNLY", "XCL1", "CTSW", "ID2", # NK
                   "CXCL13", "CXCR5", "CD69", "CCR7",
                   "TXNIP", "MAF", "TRBC1", "TRBC2",
                   "CCR6", "FCRL4", "PTPN1", # FCRL4+MCB
                 "ZBTB20", "CD27", "AHNAK", #MBC
                 "CD72", "FAM129C", "KLF2", #Naive B
                 "CD83", "CCND2", "FCER2", # Activated B
                 "XAF1", "MX1", "IFI44L", #Interferon active B
                 "PCNA", "TUBA1B", "TOP2A", # ciruclating dark GCB
                 "AICDA", "MME", "SUGCT", # Dark zone B
                 "NEIL1", "FGD6", "LMO2", #Light zone b
                 "XBP1", "JCHAIN", "MZB1" #plasmablasts
                 )

gene_list <- gene_list[gene_list %in% rownames(seurat_subset)]

plotDimRed(sample_object = seurat_subset, size = size_list[[subset_name]],
    col_by = gene_list, return_plot = TRUE)

stage_list_all <- lapply(tonsil_list, function(x) populations_dfs_new(seurat_subset,
                         x, subsample = TRUE, subsample_by = "Tonsil",
                         meta_data_col = "seurat_clusters"))
stage_df_all <- do.call("rbind", stage_list_all)

stage_df_all$sample <- factor(stage_df_all$sample)

population_plots(stage_df_all)

dev.off()

saveRDS(seurat_subset, file = output_file)
