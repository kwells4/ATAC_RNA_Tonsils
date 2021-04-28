library(ArchR)
library(openxlsx)

ATAC_path <- "/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/ArchR_analysis_scg4/All_Tonsils_high_res_new_chromvar"
save_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/revision_output"
function_path <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/scripts/functions.R"
#############
# Functions #
#############

make_barplots_ATAC <- function(ATAC_meta, sample_name,
	meta_data_list = c("Clusters_Harmony", "Clusters2"),
	color_list = NULL, position = "stack"){
	if(!setequal(names(color_list), meta_data_list)){
		stop("'color_list' must contain the same elements as 'meta_data_list'")
	}
	if(is.null(color_list)){
		color_list <- rep("", times = length(meta_data_list))
		color_list <- setNames(color_list, meta_data_list)
	}

	barplots <- lapply(meta_data_list, function(x)
		one_barplot_ATAC(ATAC_meta = ATAC_meta, subsample_list = tonsil_list,
			sample_name = sample_name, subsample_by = "Sample", meta_data_col = x,
			position = position, color = color_list[[x]]))

}

one_barplot_ATAC <- function(ATAC_meta, subsample_list, sample_name,
	subsample_by = "Sample", meta_data_col = "Clusters_Harmony",
	color = NULL, percent = TRUE, count = FALSE, sep_by = "sample",
  position = "stack"){
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
        title = paste0(sample_name, "percent"), position = position))
    }
    if (count){
      print(population_plots(stage_df_all, color = color, plot = "count",
        title = paste0(sample_name, "count"), position = position))
    }    
  } else {
  	if (percent){
  		print(population_plots_new(stage_df_all, color = color,
  			title = paste0(sample_name, "percent"), position = position))
  	}
  	if (count){
  		print(population_plots_new(stage_df_all, color = color, plot = "count",
  			title = paste0(sample_name, "count"), position = position))
  	}
  }
  return(stage_df_all)
}

populations_dfs_ATAC <- function(ATAC_meta, sample_name, subsample = FALSE,
                                subsample_by = "sample",
                                meta_data_col = "Clusters_Harmony"){
  if (subsample) {
  	print(subsample_by)
    ATAC_meta <- ATAC_meta[ATAC_meta[[subsample_by]] == sample_name, ]
  }
  stage_df <- data.frame(table(ATAC_meta[[meta_data_col]]))
  print(stage_df)
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

# Colors
new_cell_type <- c("Memory_B", "Naive_B", "Light_GC_B", "Dark_GC_B",
                "Plasmablast", "CD4_T", "CD8_T", "Myeloid_cells", "Dendritic")

colors_cell_type_new <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set1"))(length(new_cell_type))

names(colors_cell_type_new) <- new_cell_type

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



###############################
#             ATAC            #
###############################
projTonsils5 <- loadArchRProject(ATAC_path)

# add peak set
peak_set <- readRDS(paste0(ATAC_path, "/files/macs2_peaks.rds"))

projTonsils5 <- addPeakSet(projTonsils5, peakSet = peak_set, force = TRUE)

cell_type_2 <- getCellColData(projTonsils5, c("Clusters2", "RNA_high_res"))
# Rename clusters
# Low res
new_clusters_cell_type <- c(Naive_B = "Memory_B", Activated_B = "Naive_B",
                  Light_GC_B = "Light_GC_B", Dark_GC_B = "Dark_GC_B",
                  Plasma_B = "Plasmablast", CD4_T = "CD4_T", CD8_T = "CD8_T",
                  Monocyte_DC = "Myeloid_cells", Dendritic = "Dendritic")

cell_type_2$new_cell_type <- new_clusters_cell_type[cell_type_2$Clusters2]


# High res
# Rename clusters and make colors
new_clusters_high_res <- c("FCRL4+MCB" = "FCRL4_MCB",
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
cell_type_2$new_high_res <- new_clusters_high_res[cell_type_2$RNA_high_res]

levels(cell_type_2$new_cell_type) <- names(new_clusters_cell_type)
levels(cell_type_2$new_high_res) <- names(new_clusters_high_res)

projTonsils5 <- addCellColData(projTonsils5, data = cell_type_2$new_cell_type,
  name = "new_cell_type", cells = rownames(cell_type_2), force = TRUE)

projTonsils5 <- addCellColData(projTonsils5, data = cell_type_2$new_high_res,
  name = "new_high_res", cells = rownames(cell_type_2), force = TRUE)

############
# Barplots #
############

sample_data <- getCellColData(projTonsils5)

sample_data$Tonsil <- gsub("a|b", "", sample_data$Sample)

tonsil_list <- unique(sample_data$Sample)

pdf(paste0(save_dir, "/plots/allSamples_ATAC_proportions_barplot.pdf"))

percents_low_res <- make_barplots_ATAC(ATAC_meta = sample_data, sample_name = "all_samples",
	color_list = list(new_cell_type = colors_cell_type_new),
  meta_data_list = c("new_cell_type"))

dev.off()

b_meta <- sample_data %>%
  data.frame %>%
  dplyr::filter(new_high_res %in% names(b_colors))

pdf(paste0(save_dir, "/plots/allSamples_ATAC_proportions_high_res_B_barplot.pdf"))
percents_b <- make_barplots_ATAC(ATAC_meta = b_meta, sample_name = "all_samples",
  color_list = list(new_high_res = b_colors),
  meta_data_list = c("new_high_res"))

dev.off()

t_meta <- sample_data %>%
  data.frame %>%
  dplyr::filter(new_high_res %in% names(t_colors))
  
pdf(paste0(save_dir, "/plots/allSamples_ATAC_proportions_high_res_T_barplot.pdf"))
percents_t <- make_barplots_ATAC(ATAC_meta = t_meta, sample_name = "all_samples",
  color_list = list(new_high_res = t_colors),
  meta_data_list = c("new_high_res"))

dev.off()

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
                           "/files/ATAC_percents.xlsx"),
             overwrite = TRUE)


###################################
# DE for cluster 3/5 vs cluster 4 #
###################################

markerTest_PM <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters_Harmony",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4",
  bgdGroups = c("C5", "C3")
)

markerList_PM <- getMarkers(markerTest_PM, cutOff = "FDR <= 0.01 & Log2FC >= 1")

markerTest_GS <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters_Harmony",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4",
  bgdGroups = c("C5", "C3")
)

markerList_GS <- getMarkers(markerTest_GS, cutOff = "FDR <= 0.01 & Log2FC >= 1")


markerTest_GI <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "GeneIntegrationMatrix",
  groupBy = "Clusters_Harmony",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4",
  bgdGroups = c("C5", "C3")
)

markerList_GI <- getMarkers(markerTest_GI, cutOff = "FDR <= 0.01 & Log2FC >= 1")


DE_wb <- createWorkbook()

addWorksheet(DE_wb, "PeakMatrix")
writeData(DE_wb, "PeakMatrix", markerList_PM[[1]])

addWorksheet(DE_wb, "GeneScoreMatrix")
writeData(DE_wb, "GeneScoreMatrix", markerList_GS[[1]])

addWorksheet(DE_wb, "GeneIntegrationMatrix")
writeData(DE_wb, "GeneIntegrationMatrix", markerList_GI[[1]])

saveWorkbook(DE_wb,
             file = paste0(save_dir,
                           "/files/Harmony_DE.xlsx"),
             overwrite = TRUE)

################
# Correlations #
################

# Large cell types
sample_data$new_sample <- sub("a|b", "", sample_data$Sample)
sample_data$Tonsil_celltype <- paste0(sample_data$new_sample, "_",
                                      sample_data$new_cell_type)

projTonsils5 <- addCellColData(projTonsils5, data = sample_data$Tonsil_celltype,
  name = "Tonsil_celltype", cells = rownames(sample_data), force = TRUE)

correlations_wb <- createWorkbook()
# Get average expression
average_PM <- getGroupSE(
  ArchRProj = projTonsils5,
  useMatrix = "PeakMatrix",
  groupBy = "Tonsil_celltype",
)

expression <- assays(average_PM)$PeakMatrix

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(sample_data$new_cell_type), function(x){
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::contains(as.character(x)))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

correlations_all <- stats::cor(expression, use = "pairwise.complete.obs",
                           method = "spearman")
addWorksheet(correlations_wb, "all_correlations")
writeData(correlations_wb, "all_correlations", correlations_all,
          rowNames = TRUE)


## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "/files/ATAC_PM_cell_type_corr.xlsx"),
             overwrite = TRUE)


correlations_wb <- createWorkbook()

average_GS <- getGroupSE(
  ArchRProj = projTonsils5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Tonsil_celltype",
)

expression <- assays(average_GS)$GeneScoreMatrix

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(sample_data$new_cell_type), function(x){
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::contains(as.character(x)))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

correlations_all <- stats::cor(expression, use = "pairwise.complete.obs",
                           method = "spearman")
addWorksheet(correlations_wb, "all_correlations")
writeData(correlations_wb, "all_correlations", correlations_all,
          rowNames = TRUE)


## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "/files/ATAC_GS_cell_type_corr.xlsx"),
             overwrite = TRUE)

# High Res
sample_data$Tonsil_celltype_high_res <- paste0(sample_data$new_sample, "_",
                                               sample_data$new_high_res)

projTonsils5 <- addCellColData(projTonsils5, data = sample_data$Tonsil_celltype_high_res,
  name = "Tonsil_celltype_high_res", cells = rownames(sample_data), force = TRUE)

correlations_wb <- createWorkbook()

average_PM <- getGroupSE(
  ArchRProj = projTonsils5,
  useMatrix = "PeakMatrix",
  groupBy = "Tonsil_celltype_high_res",
)

expression <- assays(average_PM)$PeakMatrix

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(sample_data$new_high_res), function(x){
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::contains(as.character(x)))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

correlations_all <- stats::cor(expression, use = "pairwise.complete.obs",
                           method = "spearman")
addWorksheet(correlations_wb, "all_correlations")
writeData(correlations_wb, "all_correlations", correlations_all,
          rowNames = TRUE)

## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "/files/ATAC_PM_cell_type_high_res_corr.xlsx"),
             overwrite = TRUE)


correlations_wb <- createWorkbook()

average_GS <- getGroupSE(
  ArchRProj = projTonsils5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Tonsil_celltype_high_res",
)

expression <- assays(average_GS)$GeneScoreMatrix

expression <- expression %>%
  dplyr::na_if(0)

lapply(unique(sample_data$new_high_res), function(x){
  expression_short <- expression %>%
    data.frame %>%
    dplyr::select(dplyr::contains(as.character(x)))
  correlations <- stats::cor(expression_short, use = "pairwise.complete.obs",
                             method = "spearman")
  addWorksheet(correlations_wb, as.character(x))
  writeData(correlations_wb, as.character(x), correlations,
            rowNames = TRUE)
  return(correlations)
})

correlations_all <- stats::cor(expression, use = "pairwise.complete.obs",
                           method = "spearman")
addWorksheet(correlations_wb, "all_correlations")
writeData(correlations_wb, "all_correlations", correlations_all,
          rowNames = TRUE)

## Save workbook to working directory
saveWorkbook(correlations_wb,
             file = paste0(save_dir,
                           "/files/ATAC_GS_cell_type_high_res_corr.xlsx"),
             overwrite = TRUE)