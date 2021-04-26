library(Seurat)


# seurat_path_all <- snakemake@input[[1]]
# function_path <- snakemake@params[[1]]
# tonsil_list <- snakemake@params[[2]]
# output_file <- snakemake@output[[1]]
# save_dir <- snakemake@params[[3]]

seurat_path_all <- "output/allSamples_nomt_snakemake/files/allSamples_nomt_dubRm_namedClust_seurat.rda"
seurat_path_bcells <- "output/allSamples_nomt_BCells_snakemake/files/BCells_dubRm_namedClust_seurat.rda"
seurat_path_tcells <- "output/allSamples_nomt_TCells_snakemake/files/TCells_dubRm_namedClust_seurat.rda"

seurat_object_all <- readRDS(seurat_path_all)

seurat_bcells <- readRDS(seurat_path_bcells)

b_cell_meta <- seurat_bcells[[]]

seurat_meta <- seurat_object_all[[]]

seurat_meta_all <- merge(seurat_meta, b_cell_meta, all.x = TRUE, by = "row.names")

rownames(seurat_meta_all) <- seurat_meta_all$Row.names

identical(rownames(seurat_meta_all), rownames(seurat_object_all[[]]))

seurat_object_all$cell_type_bCell <- seurat_meta_all$cell_type_bCell

seurat_meta_all <- NULL
seurat_meta <- NULL
b_cell_meta <- NULL
seurat_bcells <- NULL



seurat_TCells <- readRDS(seurat_path_tcells)

t_cell_meta <- seurat_TCells[[]]

seurat_meta <- seurat_object_all[[]]

seurat_meta_all <- merge(seurat_meta, t_cell_meta, all.x = TRUE, by = "row.names")

rownames(seurat_meta_all) <- seurat_meta_all$Row.names

identical(rownames(seurat_meta_all), rownames(seurat_object_all[[]]))

seurat_object_all$cell_type_tCell <- seurat_meta_all$cell_type_tCell

seurat_meta_all <- NULL
seurat_meta <- NULL
t_cell_meta <- NULL
seurat_TCells <- NULL


seurat_object_meta <- seurat_object_all[[]]

seurat_object_meta$cell_type_high_res <- seurat_object_meta$cell_type_bCell

seurat_object_meta$cell_type_high_res <- as.character(seurat_object_meta$cell_type_high_res)
seurat_object_meta$cell_type_tCell <- as.character(seurat_object_meta$cell_type_tCell)

seurat_object_meta$cell_type_high_res[is.na(seurat_object_meta$cell_type_high_res)] <- 
  seurat_object_meta$cell_type_tCell[is.na(seurat_object_meta$cell_type_high_res)]

seurat_object_meta$cell_type <- as.character(seurat_object_meta$cell_type)
seurat_object_meta$cell_type_high_res[is.na(seurat_object_meta$cell_type_high_res)] <- 
  seurat_object_meta$cell_type[is.na(seurat_object_meta$cell_type_high_res)]

seurat_object_all$cell_type_high_res <- seurat_object_meta$cell_type_high_res

levels <- c(levels(seurat_object_all$cell_type_bCell),
  levels(seurat_object_all$cell_type_tCell), "Monocyte_DC", "Dendritic")

seurat_object_all$cell_type_high_res <- factor(seurat_object_all$cell_type_high_res,
  levels = levels)

saveRDS(seurat_object_all, seurat_path_all)