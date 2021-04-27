library(Seurat)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(Rmagic)
library(wesanderson)
library(openxlsx)
library(cowplot)

# seurat_path_all <- snakemake@input[[1]]
# function_path <- snakemake@params[[1]]
# tonsil_list <- snakemake@params[[2]]
# output_file <- snakemake@output[[1]]
# save_dir <- snakemake@params[[3]]

seurat_path_all <- "output/allSamples_nomt_snakemake/files/allSamples_nomt_dubRm_namedClust_seurat.rda"
seurat_path_bcells <- "output/allSamples_nomt_BCells_snakemake/files/BCells_dubRm_namedClust_seurat.rda"
seurat_path_tcells <- "output/allSamples_nomt_TCells_snakemake/files/TCells_dubRm_namedClust_seurat.rda"

# Doublet files
Tonsil1_path <- "output/mergedTonsil_snakemake/files/Tonsil1_seurat.rda"
Tonsil2_path <- "output/mergedTonsil_snakemake/files/Tonsil2_seurat.rda"
Tonsil3_path <- "output/mergedTonsil_snakemake/files/Tonsil3_seurat.rda"

doublet_path <- "RNA_seq/double_files"

tonsils <- c("Tonsil1a", "Tonsil2a", "Tonsil3a",
             "Tonsil1b", "Tonsil2b", "Tonsil3b")


function_path <- "RNA_seq/scripts/functions.R"
save_dir <- "revision_output/"

tonsil_list <- c("Tonsil1", "Tonsil2", "Tonsil3")


#############
# Functions #
#############

make_barplots <- function(seurat_object, sample_name,
                          meta_data_list = c("seurat_clusters", "cell_type"),
                          color_list = NULL, position = "stack"){
  if(!setequal(names(color_list), meta_data_list)){
    stop("'color_list' must contain the same elements as 'meta_data_list'")
  }
  if(is.null(color_list)){
    color_list <- rep("", times = length(meta_data_list))
    color_list <- setNames(color_list, meta_data_list)
  }
  
  barplots <- lapply(meta_data_list, function(x)
    one_barplot(seurat_object = seurat_object, subsample_list = tonsil_list,
                sample_name = sample_name, subsample_by = "Tonsil", meta_data_col = x,
                color = color_list[[x]], position = position))
  
}

one_barplot <- function(seurat_object, subsample_list, sample_name,
                        subsample_by = "Tonsil",
                        meta_data_col = "seurat_clusters",color = NULL,
                        percent = TRUE, count = FALSE, sep_by = "sample",
                        position = "stack"){
  if(length(color) <= 1){
    color <- NULL
  }
  stage_list_all <- lapply(subsample_list, function(x)
    populations_dfs_new(seurat_object,
                        x, subsample = TRUE, subsample_by = subsample_by,
                        meta_data_col = meta_data_col))
  
  stage_df_all <- do.call("rbind", stage_list_all)
  
  stage_df_all$sample <- factor(stage_df_all$sample)
  print(stage_df_all)
  
  if(sep_by == "sample"){
    if (percent){
      print(population_plots(stage_df_all, color = color,
                             title = paste0(sample_name, "percent"),
                             position = position))
    }
    if (count){
      print(population_plots(stage_df_all, color = color, plot = "count",
                             title = paste0(sample_name, "count"),
                             position = position))
    }    
  } else {
    if (percent){
      print(population_plots_new(stage_df_all, color = color,
                                 title = paste0(sample_name, "percent"),
                                 position = position))
    }
    if (count){
      print(population_plots_new(stage_df_all, color = color, plot = "count",
                                 title = paste0(sample_name, "count"),
                                 position = position))
    }
  }
  return(stage_df_all)
}

make_barplots_ATAC <- function(ATAC_meta, sample_name,
                               meta_data_list = c("Clusters_Harmony", "Clusters2"),
                               color_list = NULL){
  if(!setequal(names(color_list), meta_data_list)){
    stop("'color_list' must contain the same elements as 'meta_data_list'")
  }
  if(is.null(color_list)){
    color_list <- rep("", times = length(meta_data_list))
    color_list <- setNames(color_list, meta_data_list)
  }
  
  barplots <- lapply(meta_data_list, function(x)
    one_barplot_ATAC(ATAC_meta = ATAC_meta, subsample_list = tonsil_list,
                     sample_name = sample_name, subsample_by = "Tonsil", meta_data_col = x,
                     color = color_list[[x]]))
  
}

one_barplot_ATAC <- function(ATAC_meta, subsample_list, sample_name,
                             subsample_by = "Tonsil", meta_data_col = "Clusters_Harmony",
                             color = NULL, percent = TRUE, count = FALSE, sep_by = "sample"){
  if(length(color) <= 1){
    color <- NULL
  }
  stage_list_all <- lapply(subsample_list, function(x)
    populations_dfs_ATAC(ATAC_meta,
                         x, subsample = TRUE, subsample_by = subsample_by,
                         meta_data_col = meta_data_col))
  
  stage_df_all <- do.call("rbind", stage_list_all)
  
  stage_df_all$sample <- factor(stage_df_all$sample)
  
  if(sep_by == "sample"){
    if (percent){
      print(population_plots(stage_df_all, color = color,
                             title = paste0(sample_name, "percent"), position = "dodge"))
    }
    if (count){
      print(population_plots(stage_df_all, color = color, plot = "count",
                             title = paste0(sample_name, "count"), position = "dodge"))
    }    
  } else {
    if (percent){
      print(population_plots_new(stage_df_all, color = color,
                                 title = paste0(sample_name, "percent"), position = "dodge"))
    }
    if (count){
      print(population_plots_new(stage_df_all, color = color, plot = "count",
                                 title = paste0(sample_name, "count"), position = "dodge"))
    }
  }
}

populations_dfs_ATAC <- function(ATAC_meta, sample_name, subsample = FALSE,
                                 subsample_by = "sample",
                                 meta_data_col = "Clusters_Harmony"){
  if (subsample) {
    print(subsample_by)
    ATAC_meta <- ATAC_meta[ATAC_meta[[subsample_by]] == sample_name, ]
  }
  stage_df <- data.frame(table(ATAC_meta[[meta_data_col]]))
  names(stage_df) <- c("cluster", "count")
  stage_df$percent <- stage_df$count / sum(stage_df$count) * 100
  stage_df$sample <- sample_name
  return(stage_df)
}

population_plots_new <- function(stage_df_all, color = NULL, save_plot = NULL,
                                 plot = "percent", title = NULL, position = "stack"){
  if(plot == "percent"){
    plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~cluster,
                                                                    y = ~percent,
                                                                    fill = ~sample))
  } else if (plot == "count") {
    plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~cluster,
                                                                    y = ~count,
                                                                    fill = ~sample))
  } else {
    stop("plot must be either 'percent' or 'count'")
  }
  
  plot_base <- plot_base +
    # ggplot2::theme_classic() + 
    ggplot2::xlab("sample")  +
    ggplot2::geom_bar(stat = "identity", position = position) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  if(!(is.null(color))){
    plot_base <- plot_base + 
      ggplot2::scale_fill_manual(values = color, name = "cluster")
  } else {
    nColors <- length(levels(factor(stage_df_all$cluster)))
    plot_base <- plot_base +
      ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = "cluster")
  }
  if (!is.null(title)){
    plot_base <- plot_base + ggplot2::ggtitle(title)
  }
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = plot_base)
  }
  return(plot_base)
  
}


wes_palette_kw <- function (name, n, type = c("discrete", "continuous")) {
  type <- match.arg(type)
  pal <- wes_palettes[[name]]
  if (is.null(pal)) 
    stop("Palette not found.")
  if (missing(n)) {
    n <- length(pal)
  }
  if (type == "discrete" && n > length(pal)) {
    stop("Number of requested colors greater than what palette can offer")
  }
  out <- switch(type, continuous = (grDevices::colorRampPalette(pal))(n), 
                discrete = pal[1:n])
  # print(out)
  structure(out, class = "palette", name = name)
  return(out)
}

source(function_path)
set.seed(0)

data_dir <- getwd()


cell_types <- c("Naive_B", "Activated_B", "Light_GC_B", "Dark_GC_B",
                "Plasma_B", "CD4_T", "CD8_T", "Monocyte_DC", "Dendritic")

colors_cell_type <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type) <- cell_types

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

############################################
# Plots for all cells in all tonsils (RNA) #
############################################
seurat_object_all <- readRDS(seurat_path_all)
seurat_bcells <- readRDS(seurat_path_bcells)
seurat_tcells <- readRDS(seurat_path_tcells)

# Rename cell types
new_clusters <- c(Naive_B = "Memory_B", Activated_B = "Naive_B",
                  Light_GC_B = "Light_GC_B", Dark_GC_B = "Dark_GC_B",
                  Plasma_B = "Plasmablast", CD4_T = "CD4_T", CD8_T = "CD8_T",
                  Monocyte_DC = "Myeloid_cells", Dendritic = "Dendritic")

seurat_object_all$new_cell_type = new_clusters[seurat_object_all$cell_type]

new_cell_type <- c("Memory_B", "Naive_B", "Light_GC_B", "Dark_GC_B",
                "Plasmablast", "CD4_T", "CD8_T", "Myeloid_cells", "Dendritic")

colors_cell_type_new <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(length(new_cell_type))

names(colors_cell_type_new) <- new_cell_type

# Made sure these matched the manuscript
#plotDimRed(seurat_object_all, col_by = "cell_type", color = colors_cell_type,
#           ggrastr = FALSE, size = 0.5)

#plotDimRed(seurat_object_all, col_by = "new_cell_type", color = colors_cell_type_new,
#           ggrastr = FALSE, size = 0.5)

# Make barplots of percents per tonsil

seurat_object_all$new_cell_type <- factor(seurat_object_all$new_cell_type,
                                      levels = c("Memory_B", "Naive_B",
                                                 "Light_GC_B", "Dark_GC_B",
                                                 "Plasmablast", "CD4_T",
                                                 "CD8_T", "Myeloid_cells",
                                                 "Dendritic"))

pdf(paste0(save_dir, "plots/allSamples_RNA_proportions_barplot_stack.pdf"))

make_barplots(seurat_object = seurat_object_all, sample_name = "allSamples",
              meta_data_list = "new_cell_type",
              color_list = list(new_cell_type = colors_cell_type_new),
              position = "stack")

dev.off()

pdf(paste0(save_dir, "plots/allSamples_RNA_proportions_barplot_dodge.pdf"))

percents_low_res <- make_barplots(seurat_object = seurat_object_all,
                          sample_name = "allSamples",
                          meta_data_list = "new_cell_type",
                          color_list = list(new_cell_type = colors_cell_type_new),
                          position = "dodge")

dev.off()

############
# High res #
############
# Rename clusters and make colors
new_clusters <- c("FCRL4+MCB" = "FCRL4_MCB",
                  "MCB" = "MCB",
                  "Naive_B" = "Naive_B",
                  "Activated_B" = "Activated_B",
                  "Interferon_Active_B" = "IFN_active",
                  "Circulating_Dark_GC_B" = "Cycling_B",
                  "Dark_GC_B" = "Dark_GC_B",
                  "Light_GC_B" = "Light_GC_B",
                  "Plasmablasts" = "Plasmablasts",
                  "Treg" = "Treg",
                  "Tfh" = "TfH",
                  "Tfh_CXCL13" = "TfH_CXCL13",
                  "Naive_Tcells" = "Naive_CD4",
                  "Tcm_CD4" = "Tcm_CD4",
                  "Tcm_CD8" = "Tcm_CD8",
                  "CTL" = "CTL",
                  "CTL_PRF1" = "CTL_PRF1",
                  "Circulating" = "Cycling_T",
                  "NK_cells" = "NK",
                  "Monocyte_DC" = "Myeloid",
                  "Dendritic" = "Dendritic")
seurat_object_all$new_high_res = new_clusters[seurat_object_all$cell_type_high_res]

t_colors <- c("TfH_CXCL13" = "#E41A1C",
              "Naive_CD4" = "#596A98",
              "Tcm_CD4" = "#449B75",
              "Treg" = "#6B886D",
              "CTL" = "#AC5782",
              "TfH" = "#FF7F00",
              "Tcm_CD8" = "#C9992C",
              "Cycling_T" = "#C66764",
              "NK" = "#E485B7",
              "CTL_PRF1" = "#999999")


b_colors <- c("FCRL4_MCB" = "#E41A1C",
              "MCB" = "#377EB8",
              "Naive_B" = "#4DAF4A",
              "Activated_B" = "#984EA3",
              "IFN_active" = "#FF7F00",
              "Cycling_B" = "#FFFF33",
              "Dark_GC_B" = "#A65628",
              "Light_GC_B" = "#F781BF",
              "Plasmablasts" = "#999999")


# B cells
seurat_b_short <- subset(seurat_object_all, cells = colnames(seurat_bcells))
pdf(paste0(save_dir, "plots/b_high_res_RNA_proportions_barplot_stack.pdf"))

make_barplots(seurat_object = seurat_b_short, sample_name = "allSamples",
              meta_data_list = "new_high_res",
              color_list = list(new_high_res = b_colors),
              position = "stack")

dev.off()

pdf(paste0(save_dir, "plots/b_high_res_RNA_proportions_barplot_dodge.pdf"))

percents_b <- make_barplots(seurat_object = seurat_b_short,
                            sample_name = "allSamples",
                            meta_data_list = "new_high_res",
                            color_list = list(new_high_res = b_colors),
                            position = "dodge")
dev.off()

# T cells
seurat_t_short <- subset(seurat_object_all, cells = colnames(seurat_tcells))
pdf(paste0(save_dir, "plots/t_high_res_RNA_proportions_barplot_stack.pdf"))

make_barplots(seurat_object = seurat_t_short, sample_name = "allSamples",
              meta_data_list = "new_high_res",
              color_list = list(new_high_res = t_colors),
              position = "stack")

dev.off()

pdf(paste0(save_dir, "plots/t_high_res_RNA_proportions_barplot_dodge.pdf"))

percents_t <- make_barplots(seurat_object = seurat_t_short,
                            sample_name = "allSamples",
                            meta_data_list = "new_high_res",
                            color_list = list(new_high_res = t_colors),
                            position = "dodge")
dev.off()

# Write to excel file
percents_wb <- createWorkbook()
addWorksheet(percents_wb, "low_resolution")
writeData(percents_wb, "low_resolution", percents_low_res[[1]])

addWorksheet(percents_wb, "b_cell_high_res")
writeData(percents_wb, "b_cell_high_res", percents_b[[1]])

addWorksheet(percents_wb, "t_cell_high_res")
writeData(percents_wb, "t_cell_high_res", percents_t[[1]])

## Save workbook
saveWorkbook(percents_wb,
             file = paste0(save_dir,
                           "files/RNA_percents.xlsx"),
             overwrite = TRUE)

################
# Correlations #
################

# Large cell types
correlations_wb <- createWorkbook()

seurat_object_all$Tonsil_celltype <- paste0(seurat_object_all$Tonsil, "_",
                                           seurat_object_all$new_cell_type)
Idents(seurat_object_all) <- "Tonsil_celltype"

expression <- AverageExpression(seurat_object_all)
expression <- expression$SCT

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(seurat_object_all$new_cell_type), function(x){
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::matches(paste0("Tonsil[0-9]_", as.character(x))))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "files/cell_type_corr.xlsx"),
             overwrite = TRUE)


# High Res cell types
correlations_wb <- createWorkbook()

seurat_object_all$Tonsil_celltype_high_res <- paste0(seurat_object_all$Tonsil, "_",
                                            seurat_object_all$new_high_res)
Idents(seurat_object_all) <- "Tonsil_celltype_high_res"

expression <- AverageExpression(seurat_object_all)
expression <- expression$SCT

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(seurat_object_all$new_high_res), function(x){
  print(x)
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::matches(paste0("Tonsil[0-9]_", as.character(x), "$")))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "files/cell_type_corr_high_res.xlsx"),
             overwrite = TRUE)


####################
# Doublet checking #
####################
# My objects
Tonsil1_obj <- readRDS(Tonsil1_path)
Tonsil2_obj <- readRDS(Tonsil2_path)
Tonsil3_obj <- readRDS(Tonsil3_path)

# Hamish doublet scores
plot_doublets <- function(tonsil, tonsil_obj){
  df_col <- colnames(tonsil_obj[[]])[grepl("DF.classifications",
                                           colnames(tonsil_obj[[]]))]
  a_csv <- read.csv(paste0(doublet_path, "/", tonsil, "a__scrublet_scores.csv"))
  b_csv <- read.csv(paste0(doublet_path, "/", tonsil, "b__scrublet_scores.csv"))
  
  a_csv$barcode <- paste0(tonsil, "a_", a_csv$barcode)
  a_csv$barcode <- sub("-1", "", a_csv$barcode)
  
  b_csv$barcode <- paste0(tonsil, "b_", b_csv$barcode)
  b_csv$barcode <- sub("-1", "", b_csv$barcode)
  
  scrublet_scores <- rbind(a_csv, b_csv)
  scrublet_scores <- scrublet_scores %>%
    dplyr::filter(barcode %in% rownames(tonsil_obj[[]]))
  
  scrublet_scores <- scrublet_scores[order(match(scrublet_scores$barcode,
                                                 rownames(tonsil_obj[[]]))),]
  
  tonsil_obj$doublet <- scrublet_scores$predicted_doublet
  
  db_cols <- c("Doublet" = "red", "Singlet" = "grey",
               "True" = "red", "False" = "grey")
  
  plot1 <- plotDimRed(tonsil_obj, col_by = df_col,
                      color = db_cols, show_legend = FALSE)[[1]] +
    ggplot2::ggtitle("Doublet Finder")
  
  plot2 <- plotDimRed(tonsil_obj, col_by = "doublet",
                      color = db_cols, show_legend = TRUE)[[1]] +
    ggplot2::ggtitle("Scrublet")
  
  doublet_plot <- plot_grid(plot1, plot2,
                            nrow = 1, ncol = 2,
                            align = "hv",
                            axis = "tb",
                            labels = c("A", "B"))
  pdf(paste0(save_dir, "plots/doublets_", tonsil, ".pdf"), width = 11, height = 8)
  print(doublet_plot)
  dev.off()
}

plot_doublets(tonsil = "Tonsil1", tonsil_obj = Tonsil1_obj)
plot_doublets(tonsil = "Tonsil2", tonsil_obj = Tonsil2_obj)
plot_doublets(tonsil = "Tonsil3", tonsil_obj = Tonsil3_obj)
