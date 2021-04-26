library(Seurat)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggrastr)
library(ggpubr)
library(gridExtra)
library(Rmagic)
library(wesanderson)

# seurat_path_all <- snakemake@input[[1]]
# function_path <- snakemake@params[[1]]
# tonsil_list <- snakemake@params[[2]]
# output_file <- snakemake@output[[1]]
# save_dir <- snakemake@params[[3]]

seurat_path_all <- "output/allSamples_nomt_snakemake/files/allSamples_nomt_dubRm_namedClust_seurat.rda"
seurat_path_gc <- "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_dubRm_namedClust_seurat.rda"
seurat_path_bcells <- "output/allSamples_nomt_BCells_snakemake/files/BCells_dubRm_namedClust_seurat.rda"
seurat_path_nonLymphoid <- "output/allSamples_nomt_nonLymphoid_snakemake/files/nonLymphoid_dubRm_seurat.rda"
seurat_path_tcells <- "output/allSamples_nomt_TCells_snakemake/files/TCells_dubRm_namedClust_seurat.rda"
ATAC_path <- "/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/ArchR_analysis_scg4/All_Tonsils_high_res_new_chromvar"
seurat_path_gc_new <- "output/gc_seeds/allSamples_nomt_GCBcells_snakemake_17/files/GCBcells_dubRm_seurat.rda"

function_path <- "scripts/functions.R"
save_dir <- "final_plots/"

tonsil_list <- c("Tonsil1", "Tonsil2", "Tonsil3")

all_markers <- c("TOP2A", "PCNA", "MKI67", "ICAM1", "SEMA7A", "BACH2", "MS4A1",
	             "STX7", "HMMR", "PLK1", "XBP1", "PRDM1", "CD27", "AIM2", "SCIMP",
	             "CXCR6", "FCRL4", "NEAT1", "CD72", "TCL1A", "NFKB1", "CCND2",
	             "IRF4", "CD83", "EGR3", "TCF4", "CCR7", "CXCL13", "FKBP5",
	             "LEF1", "KLF2", "KLF6", "RORA", "IL7R", "NKG7", "ANXA1", "TXNIP",
	             "GZMA", "CCL4", "CLDND1", "BANK1", "ANXA2", "XCL1", "CD81", "CLEC4C",
	             "CD14", "FCGR2A", "IFNGR1", "CCR6", "PSME2", "CD69", "CXCR4", "CXCR5",
	             "VIM", "CD44", "PLAC8", "EZR", "DUSP2", "LY86", "JCHAIN", "MZB1", "MX1",
	             "CD52", "FYB1")

markerGenes <- c(
    "SEMA7A", "BACH2", #Germinal center dz
    "MS4A1", "STX7",  # germinal center lz
    "XBP1", "PRDM1", # plasmablasts
    "CD27", "AIM2", # MBC
    # FCRL4 MBC
    "CD72", "TCL1A", # Naive B cells
    # preGC/Acrivated B
    "LEF1", "KLF2", "KLF6", "NKG7",
    "BANK1", 
    "CD81", "CD8A", "CD4", "CD3D", "CD3G" # T cell # DC
    )


#############
# Functions #
#############

make_barplots <- function(seurat_object, sample_name,
	meta_data_list = c("seurat_clusters", "cell_type"),
	color_list = NULL){
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
			color = color_list[[x]]))

}

one_barplot <- function(seurat_object, subsample_list, sample_name,
	subsample_by = "Tonsil", meta_data_col = "seurat_clusters",
	color = NULL, percent = TRUE, count = FALSE, sep_by = "sample"){
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

colors_cell_type_new <- grDevices::colorRampPalette(
	RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type_new) <- cell_types

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

############################################
# Plots for all cells in all tonsils (RNA) #
############################################
seurat_object_all <- readRDS(seurat_path_all)
adt_rna_list <- c("IGHM", "IGHD", "CD44", "MME", "MS4A1",  "CD38",  "CD27",
  "CD3G", "CD8A", "CD4")

# Run magic
seurat_object_all <- magic(seurat_object_all, genes = adt_rna_list)


pdf(paste0(save_dir, "allSamples_gene_imputed.pdf"))
DefaultAssay(seurat_object_all) <- "MAGIC_SCT"
pal <- wes_palette_kw("Zissou1", 21, type = "continuous")
plots <- plotDimRed(seurat_object_all, col_by = adt_rna_list, ggrastr = TRUE,
  wesanderson = TRUE, color = pal)
plots
# plots <- lapply(1:length(plots), function(x){
#   plot <- plots[[x]]
#   plot <- plot + scale_color_viridis(option = "plasma")
#   print(plot)
#   return(plot)
#   })

dev.off()

# Run magic on ADTs
DefaultAssay(seurat_object_all) <- "ADT"
seurat_object_all <- magic(seurat_object_all)


pdf(paste0(save_dir, "allSamples_ADT_imputed.pdf"))
DefaultAssay(seurat_object_all) <- "MAGIC_ADT"
#pal <- wes_palette_kw("Zissou1", 21, type = "continuous")
plots <- plotDimRed(seurat_object_all, col_by = rownames(seurat_object_all),
  ggrastr = TRUE)
#plots
# plots <- lapply(1:length(plots), function(x){
#   plot <- plots[[x]]
#   plot <- plot + scale_color_viridis()
#   print(plot)
#   return(plot)
#   })

plots

dev.off()

DefaultAssay(seurat_object_all) <- "SCT"

# Colors for cluster
nclusters <- length(levels(seurat_object_all$seurat_clusters))
colors_cluster <- grDevices::colorRampPalette(
	RColorBrewer::brewer.pal(9, "Set3"))(nclusters)

names(colors_cluster) <- 0:(nclusters - 1)

# Figure 1
# UMAP of all samples, colored by cell type
pdf(paste0(save_dir, "allSamples_RNA_UMAP.pdf"),
	height = 12, width = 10)

plotDimRed(seurat_object_all, col_by = "cell_type", color = colors_cell_type_new,
	ggrastr = TRUE, size = 0.5)

dev.off()

# UMAP of all samples, colored by cluster
pdf(paste0(save_dir, "allSamples_RNA_UMAP_cluster.pdf"),
	height = 12, width = 12)

plotDimRed(seurat_object_all, col_by = "num_cell_type", ggrastr = TRUE,
	size = 0.5, color = colors_cluster)


dev.off()

# UMAP of all samples, colored by ADT
pdf(paste0(save_dir, "allSamples_RNA_UMAP_ADT.pdf"),
	height = 8, width = 25)
adt_list <- rownames(seurat_object_all[["ADT"]])

adt_list <- adt_list[!(adt_list %in% 
	c("CD184-(CXCR4)-", "CD185-(CXCR5)-", "IgG-Fc-"))]

adt_list <- c("IgM-", "IgD-", "CD44-", "CD10-", "CD20-",  "CD38-",  "CD27-",
	"CD3-", "CD8a-", "CD4-")

plot_list <- plotDimRed(seurat_object_all, col_by = adt_list, ggrastr = TRUE,
	size = 0.5)

plot_list2 <- lapply(1:length(plot_list), function(x){
	plot <- plot_list[[x]]
	plot <- plot + ggplot2::ggtitle(adt_list[[x]])
	return(plot)
	})

ggarrange(plotlist = plot_list2, ncol = 5, nrow = 2)

dev.off()


adt_list <- c("IgM-", "IgD-", "CD44-", "CD10-", "CD20-", "CD184-(CXCR4)-",
              "CD185-(CXCR5)-",  "CD38-",  "CD27-", "CD3-", "CD8a-", "CD4-")

adt_list_fetch <- paste0("adt_", adt_list)
adt_rna_list <- c("IGHM", "IGHD", "CD44", "MME", "MS4A1", "CXCR4", "CXCR5",
                  "CD38",  "CD27", "CD3G", "CD8A", "CD4")

# Heatmap for both

plot_data_adt <- t(FetchData(object = seurat_object_all, vars = adt_list_fetch))
plot_data_rna <- t(FetchData(object = seurat_object_all, vars = adt_rna_list))
annotations <- FetchData(object = seurat_object_all, vars = "cell_type")

adt_list_all <- c("IGHM", "adt_IgM-", "IGHD", "adt_IgD-", "CD44", "adt_CD44-",
                  "MME", "adt_CD10-", "MS4A1", "adt_CD20-", "CXCR4",
                  "adt_CD184-(CXCR4)-", "CXCR5", "adt_CD185-(CXCR5)-", "CD38",
                  "adt_CD38-", "CD27", "adt_CD27-", "CD3G", "adt_CD3-", "CD8A",
                  "adt_CD8a-", "CD4", "adt_CD4-")

plot_data_rna <- t(scale(t(plot_data_rna)))

plot_data_adt <- t(scale(t(plot_data_adt)))

cluster_color <- colors_cell_type_new

if(all.equal(colnames(plot_data_rna), colnames(plot_data_adt))){
  plot_data_all <- rbind(plot_data_rna, plot_data_adt)
} else {
  print("ADT and RNA have different cells!")
}

annotations$cell_type <- factor(annotations$cell_type,
  levels = c("Naive_B", "Activated_B", "Light_GC_B", "Dark_GC_B",
             "Plasma_B", "CD4_T", "CD8_T", "Monocyte_DC", "Dendritic"))
annotations_s <- annotations[order(annotations$cell_type), , drop = FALSE]
plot_data_s <- plot_data_all[order(match(rownames(plot_data_all), adt_list_all)),
  order(match(colnames(plot_data_all), rownames(annotations_s)))]
# plot_data_s <- plot_data[order(match(rownames(plot_data),
#                                      rownames(cluster_rows_df))),
#   order(match(colnames(plot_data), rownames(annotations_s)))]

# plot_data_s2 <- plot_data[order(match(rownames(plot_data),
#                                       rownames(cluster_rows_df))),
#   order(match(colnames(plot_data), rownames(cluster_columns_df)))]

# cluster_columns_df$cell_cluster <- factor(cluster_columns_df$cell_cluster)
# cluster_rows_df$gene_cluster <- factor(cluster_rows_df$gene_cluster)
pdf(paste0(save_dir, "/allSamples_RNA_ADT_gene_expression_heatmap.pdf"),
  height = 5, width = 10)


all_colors <- list(cell_type = cluster_color)
  
  # This comes from the seurat DoHeatmap function
plot_data_s <- MinMax(plot_data_s, -2.5, 2.5)
pheatmap(mat = plot_data_s,
               show_colnames = FALSE,
               annotation_col = annotations_s,
               #annotation_row = cluster_rows_df,
               annotation_colors = all_colors,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               color = viridis(10),
               main = "ADT RNA Expression")

dev.off()

plot_data_adt_s <- plot_data_adt[order(match(rownames(plot_data_adt), adt_list_fetch)),
  order(match(colnames(plot_data_adt), rownames(annotations_s)))]

plot_data_rna_s <- plot_data_rna[order(match(rownames(plot_data_rna), adt_rna_list)),
  order(match(colnames(plot_data_rna), rownames(annotations_s)))]

pdf(paste0(save_dir, "/allSamples_RNA_ADT_gene_expression_separate_heatmap.pdf"),
  height = 5, width = 10)
  
  # This comes from the seurat DoHeatmap function
plot_data_adt_s <- MinMax(plot_data_adt_s, -2.5, 2.5)
pheatmap(mat = plot_data_adt_s,
               show_colnames = FALSE,
               annotation_col = annotations_s,
               #annotation_row = cluster_rows_df,
               annotation_colors = all_colors,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               color = inferno(10),
               main = "ADT expression")

plot_data_rna_s <- MinMax(plot_data_rna_s, -2.5, 2.5)
pheatmap(mat = plot_data_rna_s,
               show_colnames = FALSE,
               annotation_col = annotations_s,
               #annotation_row = cluster_rows_df,
               annotation_colors = all_colors,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               color = viridis(10),
               main = "RNA expression")

dev.off()

# Dotplot of RNA markers
cell_groups <- c(Activated_B = "B_cell",
                 Naive_B = "B_cell",
                 Dark_GC_B = "B_cell",
                 Light_GC_B = "B_cell",
                 Plasma_B = "B_cell",
                 CD4_T = "T_cell",
                 CD8_T = "T_cell",
                 Monocyte_DC = "non_lymphoid",
                 Dendritic = "non_lymphoid"
                 )

seurat_object_all$cell_type <- factor(seurat_object_all$cell_type,
    levels = names(cell_groups))
slot <- "scale.data"
Idents(seurat_object_all) <- "cell_type"
sample_markers_cell_type <- FindAllMarkers(seurat_object_all, only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = 0.25, slot = slot)

sample_markers_ct_short <- sample_markers_cell_type[sample_markers_cell_type$pct.2 < 0.5,]


top10_cell_type <- sample_markers_ct_short %>% group_by(cluster) %>% 
    top_n(n = 5, wt = avg_diff)

zohar_markers <- c("CLEC4C", "CD14", "FCGR2A", "IFNGR1", "CD94", "GZMB", #Dendritic
	"LYZ", "VIM", "CST3", #monocytes
	"KLRG1", "GZMA", "IL32", "NKG7", "CCL5", #CD8T
	"IL6R", "CD28", "CD3D", #CD4T
	"PRDM1", "XBP1", "JCHAIN", "MZB1", #plasma B
	"LPP", "NEIL1", "FCRL3", "CD40", "MME", "HMMR", "LMO2", #Light GC
	"TUBA1B", "HMGB2", "STMN1", "TOP2A", #Dark GC
	"FCER2", "CD72", "TCL1A", "BANK1", # Naive B
	"FCRL4", "CCR6", "PLAC8", "TNFRSF13B", "CD83" # Activated B
	)

pdf(paste0(save_dir, "allSamples_RNA_markers_dotplot.pdf"),
	height = 20, width = 10)

plot1 <- DotPlot(object = seurat_object_all, features = zohar_markers,
    group.by = "cell_type",
    cols = c("lightgrey", "blue"))
plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + coord_flip()

dev.off()

# Make barplots of percents per tonsil

seurat_object_all$cell_type <- factor(seurat_object_all$cell_type,
  levels = c("Naive_B", "Activated_B", "Light_GC_B", "Dark_GC_B",
             "Plasma_B", "CD4_T", "CD8_T", "Monocyte_DC", "Dendritic"))

pdf(paste0(save_dir, "allSamples_RNA_proportions_barplot.pdf"))

make_barplots(seurat_object = seurat_object_all, sample_name = "allSamples",
	color_list = list(cell_type = colors_cell_type_new, 
		seurat_clusters = colors_cluster))

dev.off()
##############################################################
# Add in percents of higher resolution cell types from B and T cell subsets
pdf(paste0(save_dir, "allSamples_RNA_UMAP_b_cell_cluster.pdf"),
  height = 12, width = 12)

plotDimRed(seurat_object_all, col_by = "cell_type_bCell", ggrastr = TRUE,
  size = 0.5)


dev.off()

pdf(paste0(save_dir, "allSamples_RNA_proportions_barplot_b_cells.pdf"))


one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
      sample_name = "allSamples", subsample_by = "Tonsil",
      meta_data_col = "cell_type_bCell")


dev.off()


pdf(paste0(save_dir, "allSamples_RNA_UMAP_T_cell_cluster.pdf"),
  height = 12, width = 12)

plotDimRed(seurat_object_all, col_by = "cell_type_tCell", ggrastr = TRUE,
  size = 0.5)


dev.off()

pdf(paste0(save_dir, "allSamples_RNA_proportions_barplot_t_cells.pdf"))


one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
      sample_name = "allSamples", subsample_by = "Tonsil",
      meta_data_col = "cell_type_tCell")


dev.off()

colors_tonsil <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(8, "Dark2"))(3)

pdf(paste0(save_dir, "allSamples_RNA_proporation_barplot_all_cells.pdf"),
  height = 12, width = 12)

one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
      sample_name = "allSamples", subsample_by = "Tonsil",
      meta_data_col = "cell_type_high_res", sep_by = "cluster",
      color = colors_tonsil)

dev.off()

pdf(paste0(save_dir, "allSamples_RNA_proportion_barplot_large_clusters.pdf"))

one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
      sample_name = "allSamples", subsample_by = "Tonsil",
      meta_data_col = "cell_type", sep_by = "cluster",
      color = colors_tonsil)

dev.off()

pdf(paste0(save_dir, "allSamples_RNA_proportion_barplot_high_res_RNA_clusters.pdf"))
one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil",
  meta_data_col = "cell_type_high_res",
  color = colors_tonsil, sep_by = "cluster")

one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil",
  meta_data_col = "cell_type_high_res",
  sep_by = "sample")
dev.off()

pdf(paste0(save_dir, "allSamples_RNA_proportion_barplot_large_clusters.pdf"))
one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil",
  meta_data_col = "cell_type",
  color = colors_tonsil, sep_by = "cluster")

one_barplot(seurat_object = seurat_object_all, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil",
  meta_data_col = "cell_type",
  color = colors_cell_type_new, sep_by = "sample")

dev.off()



# Figure 2
# # Markers on UMAP dimensional reduction
# markerGenes <- c("MS4A1", "PRDM1", "TCL1A", "CD8A", "CD4", "CD3D", "POU2F2",
#   "RORA", "TCF4", "IRF4")

# pdf(paste0(save_dir, "all_sample_marker_genes.pdf"), height = 8, width = 8)

# plotDimRed(seurat_object_all, col_by = "cell_type", color = colors_cell_type_new,
#   ggrastr = TRUE, size = 0.5, show_legend = FALSE)

# plotDimRed(seurat_object_all, col_by = markerGenes, ggrastr = TRUE)

# dev.off()


###############################
# Plots for all B cells (RNA) #
###############################
seurat_bcells <- readRDS(seurat_path_bcells)

save_dir_bcell <- "output/allSamples_nomt_BCells_snakemake/"


# plotDimRed(seurat_bcells, col_by = all_markers,
# 	save_plot = paste0(save_dir_bcell, "images/UMAP_markers.pdf"))

# new_markers <- c("XBP1", "JCHAIN", "MZB1")

# pdf(paste0(save_dir_bcell, "images/dotplot_markers.pdf"),
# 	height = 10, width = 8)

# plot1 <- DotPlot(object = seurat_bcells, features = new_markers,
#     group.by = "seurat_clusters",
#     cols = c("lightgrey", "blue"))
# plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

# dev.off()


pdf(paste0(save_dir, "allSamples_BCells_RNA_UMAP.pdf"),
	height = 12, width = 10)

plotDimRed(seurat_bcells, col_by = "cell_type", color = colors_cell_type_new,
	ggrastr = TRUE, size = 0.5)

dev.off()

pdf(paste0(save_dir, "allSamples_BCells_RNA_UMAP_cluster.pdf"),
	height = 12, width = 12)

plotDimRed(seurat_bcells, col_by = "new_cluster", ggrastr = TRUE,
	size = 0.5)

dev.off()

pdf(paste0(save_dir, "allSamples_BCells_RNA_UMAP_BCell_cell_type.pdf"),
	height = 12, width = 12)

plotDimRed(seurat_bcells, col_by = "cell_type_bCell", ggrastr = TRUE,
	size = 0.5)

dev.off()


# slot <- "scale.data"
# sample_markers_cell_type <- FindAllMarkers(seurat_bcells, only.pos = TRUE,
#     min.pct = 0.25, logfc.threshold = 0.25, slot = slot)

# sample_markers_ct_short <- sample_markers_cell_type[sample_markers_cell_type$pct.2 < 0.5,]


# top10_cell_type <- sample_markers_ct_short %>% group_by(cluster) %>% 
#     top_n(n = 10, wt = avg_diff)

# top_10 <- sample_markers_cell_type %>% group_by(cluster) %>%
# 	top_n(n = 10, wt = avg_diff)

# genes <- unique(top10_cell_type$gene)

zohar_markers <- c("CCR6", "FCRL4", "PTPN1", # FCRL4+MCB
	               "ZBTB20", "CD27", "AHNAK", "TNFRSF13B", #MBC
	               "CD72", "FAM129C", "KLF2", #Naive B
	               "CD83", "CCND2", "FCER2", # Activated B
	               "XAF1", "MX1", "IFI44L", #Interferon active B
	               "PCNA", "TUBA1B", "TOP2A", # ciruclating dark GCB
	               "AICDA", "MME", "SUGCT", # Dark zone B
	               "NEIL1", "FGD6", "LMO2", #Light zone b
	               "XBP1", "JCHAIN", "MZB1" #plasmablasts
	               )

# zohar_markers <- rev(zohar_markers)

pdf(paste0(save_dir, "allSamples_BCells_RNA_markers_dotplot.pdf"),
	height = 7.5, width = 5)

plot1 <- DotPlot(object = seurat_bcells, features = zohar_markers,
    group.by = "cell_type_bCell",
    cols = c("lightgrey", "red"))
plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + coord_flip()

dev.off()


###########################
# Plots for T cells (RNA) #
###########################

seurat_TCells <- readRDS(seurat_path_tcells)

t_cell_dir <- "output/allSamples_nomt_TCells_snakemake"
# Merge PRF1 CTL
# levels <- levels(seurat_TCells$cell_type_tCell)
# levels <- levels[-8]
# seurat_TCells$cell_type_tCell <- sub("CTL_PRF1", "CTL", 
#   seurat_TCells$cell_type_tCell)
# seurat_TCells$cell_type_tCell <- factor(seurat_TCells$cell_type_tCell,
#   levels = levels)

pdf(paste0(save_dir, "allSamples_TCells_RNA_UMAP.pdf"),
	height = 12, width = 10)

plotDimRed(seurat_TCells, col_by = "cell_type", color = colors_cell_type_new,
	ggrastr = TRUE, size = 1)

dev.off()

pdf(paste0(save_dir, "allSamples_TCells_RNA_UMAP_cluster.pdf"),
	height = 12, width = 12)

plotDimRed(seurat_TCells, col_by = "new_cluster", ggrastr = TRUE,
	size = 1)

dev.off()

pdf(paste0(save_dir, "allSamples_TCells_RNA_UMAP_named_cluster.pdf"),
  height = 12, width = 12)

plotDimRed(seurat_TCells, col_by = "cell_type_tCell", ggrastr = TRUE,
  size = 1)

dev.off()

# slot <- "scale.data"
# sample_markers_cell_type <- FindAllMarkers(seurat_TCells, only.pos = TRUE,
#     min.pct = 0.25, logfc.threshold = 0.25, slot = slot)

# sample_markers_ct_short <- sample_markers_cell_type[sample_markers_cell_type$pct.2 < 0.5,]


# top10_cell_type <- sample_markers_ct_short %>% group_by(cluster) %>% 
#     top_n(n = 10, wt = avg_diff)

# genes <- top10_cell_type$gene

zohar_markers <- c("PRDM1", "RORA", "RGS1", "CYTOR", "FOXP3", # T reg
                   "CXCL13", # Tfh CXCL13
                   "CXCR5", "FKBP5", "PDCD1", "BCL6", "PASK", "ST8SIA1", # Tfh
                   "IL7R", "NABP1", "ANK3", "ZFP36L2", # Naive T cells
                   "CCR7", "LEF1", "ITGA6", "KLF2", "SELL", #Tcm
                   "GZMA", "CCL5", "GZMK", "CST7", "CCL4", #CTL
                   "PRF1", "ANXA1", #CTL_PRF1
                   "MKI67", "PCNA", "TOP2A", #Circulating
                   "GNLY", "XCL1", "CTSW", "ID2" # NK
                   )

pdf(paste0(save_dir, "allSamples_TCells_RNA_markers_dotplot.pdf"),
	height = 8, width = 5)

plot1 <- DotPlot(object = seurat_TCells, features = zohar_markers,
    group.by = "cell_type_tCell",
    cols = c("lightgrey", "red"))
plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + coord_flip()

dev.off()


###############################
#             ATAC            #
###############################

# UMAP of markers of interest
library(ArchR)

projTonsils5 <- loadArchRProject(ATAC_path)

# add peak set
peak_set <- readRDS(paste0(ATAC_path, "/files/macs2_peaks.rds"))

projTonsils5 <- addPeakSet(projTonsils5, peakSet = peak_set, force = TRUE)

cell_type_2 <- getCellColData(projTonsils5, c("Clusters2", "RNA_high_res"))


cell_type_2$new_clusters2 <- cell_type_2$Clusters2

cell_type_2$new_clusters2[cell_type_2$RNA_high_res == "Light_GC_B"] <- "Light_GC_B"

cell_type_2$new_clusters2[cell_type_2$RNA_high_res == "Circulating_Dark_GC_B"] <- "Dark_GC_B"

levels(cell_type_2$new_clusters2) <- cell_types

projTonsils5 <- addCellColData(projTonsils5, data = cell_type_2$new_clusters2,
  name = "new_clusters2", cells = rownames(cell_type_2), force = TRUE)

p1 <- plotEmbedding(projTonsils5, colorBy = "cellColData", name = "new_clusters2",
    pal = colors_cell_type_new, embedding = "UMAP_Harmony")

plotPDF(p1, name = "Plot-UMAP-New-Clusters.pdf", ArchRProj = projTonsils5,
  addDOC = FALSE, width = 5, height = 5)


sample_data <- getCellColData(projTonsils5)

sample_data$Tonsil <- gsub("a|b", "", sample_data$Sample)

nclusters_ATAC <- length(unique(sample_data$Clusters_Harmony))
colors_cluster_ATAC <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(9, "Set3"))(nclusters_ATAC)

names(colors_cluster_ATAC) <- paste0("C", 1:(nclusters_ATAC))
sample_data$new_clusters2 <- factor(sample_data$new_clusters2,
    levels = names(cell_groups))


pdf("/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/final_plots/allSamples_ATAC_proportions_barplot.pdf")

make_barplots_ATAC(ATAC_meta = sample_data, sample_name = "all_samples",
	color_list = list(Clusters2 = colors_cell_type_new, 
		Clusters_Harmony = colors_cluster_ATAC))
dev.off()


pdf(paste0(save_dir, "allSamples_ATAC_proportion_barplot_high_res_RNA_clusters.pdf"))
one_barplot_ATAC(ATAC_meta = sample_data, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil", meta_data_col = "RNA_high_res",
  color = colors_tonsil, sep_by = "cluster")

one_barplot_ATAC(ATAC_meta = sample_data, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil", meta_data_col = "RNA_high_res",
  sep_by = "sample")

dev.off()

pdf(paste0(save_dir, "allSamples_ATAC_proportion_barplot_RNA_clusters.pdf"))
one_barplot_ATAC(ATAC_meta = sample_data, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil", meta_data_col = "new_clusters2",
  color = colors_tonsil, sep_by = "cluster")

one_barplot_ATAC(ATAC_meta = sample_data, subsample_list = tonsil_list,
  sample_name = "allSamples", subsample_by = "Tonsil", meta_data_col = "new_clusters2",
  color = colors_cell_type_new, sep_by = "sample")
dev.off()

# Footprints 
# footprint_markers <- c("DNMT1", "SMAD5", "KLF5", "KLF7", "ZFX", "SP2", "EGR1",
# 	"SP1", "KLF6", "CTCFL", "ZBTB7A", "HMGA1", "KLF13", "KLF10", "CEBPB",
# 	"BCL11A", "SMAD4", "HES4", "E2F7", "HMGA2", "MEIS2", "E2F3", )

footprint_markers <- c("FOXO1")

motifs <- footprint_markers
markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
  useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs


#################

# Get correlations with the RNA expression preditions

# Identify devient motifs
seGroupMotif <- getGroupSE(ArchRProj = projTonsils5,
  useMatrix = "MotifMatrix", groupBy = "Clusters2")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
    ArchRProj = projTonsils5,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGIM_MM <- correlateMatrices(
    ArchRProj = projTonsils5,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name,
  rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name,
  rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

TF_regulator <- sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

TF_regulator_ord <- c("POU2F2", "MEF2C", "PAX5", "SPIB", "BCL11A",
  "STAT1", "CTCF", "HMGA1", "TCF3", "SPI1", "POU2F1", "PKNOX1",
  "LMO2", "EBF1", "IRF4", "NFIA", "CEBPA", "EOMES", "CEBPD", "FOSL2",
  "CEBPB", "RUNX1")

pdf(paste0(save_dir, "allSamples_TFRegulators_from_ATAC_dotplot.pdf"),
  height = 10, width = 8)

plot1 <- DotPlot(object = seurat_object_all, features = TF_regulator_ord,
    group.by = "cell_type_high_res",
    cols = c("lightgrey", "blue"))
plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + coord_flip()

dev.off()


############################



p <- plotEmbedding(
    ArchRProj = projTonsils5, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils5)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-GC-Genes-Deviations-W-Imputation.pdf", 
    ArchRProj = projTonsils5, 
    addDOC = FALSE, width = 5, height = 5)
motifPositions <- getPositions(projTonsils5)

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
  value = TRUE)))
markerMotifs

seFoot_bcell <- getFootprints(
  ArchRProj = projTonsils5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2",
  useGroups = c("Dark_GC_B", "Light_GC_B", "Activated_B", "Naive_B", "Plasma_B")
)



p2 <- plotFootprints(
  seFoot = seFoot_bcell,
  ArchRProj = projTonsils5, 
  normMethod = "Subtract",
  plotName = "FOXO1-footprints-subtract-bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)


markerTest <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Dark_GC_B",
  bgdGroups = "Light_GC_B"
)

markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerListLight <- getMarkers(markerTest, cutOff = "FDR <= 0.01 & Log2FC <= -1.25")
markerTestLIght <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Light_GC_B",
  bgdGroups = "Dark_GC_B"
)

markerListLight <- getMarkers(markerTestLIght, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")


markerMotifTest <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "MotifMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Dark_GC_B",
  bgdGroups = "Light_GC_B"
)

markerMotifList <- getMarkers(markerMotifTest, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerMotifListLight <- getMarkers(markerMotifTest, cutOff = "FDR <= 0.01 & Log2FC <= -1.25")


markerMotifTestLight <- getMarkerFeatures(
  ArchRProj = projTonsils5, 
  useMatrix = "MotifMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Light_GC_B",
  bgdGroups = "Dark_GC_B"
)

markerMotifListLight <- getMarkers(markerMotifTestLight, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")


# > markerMotifListLight[["Light_GC_B"]]$name
#  [1] "SP4_180"             "KLF14_251"           "KLF11_847"          
#  [4] "KLF2_846"            "SP8_226"             "SP7_241"            
#  [7] "SP9_283"             "SP6_275"             "SP5_279"            
# [10] "SP3_247"             "KLF16_205"           "KLF13_238"          
# [13] "KLF10_826"           "KLF7_189"            "NRF1_805"           
# [16] "ZFP42_848"           "KLF4_208"            "NHLH2_80"           
# [19] "ZNF143_231"          "SIX5_540"            "SP1_267"            
# [22] "TFAP4_23"            "ARNTL_50"            "BHLHE41_44"         
# [25] "USF1_60"             "TFE3_18"             "ARNT_57"            
# [28] "TFEB_32"             "TCFL5_25"            "BHLHE40_51"         
# [31] "TFEC_27"             "ARNTL2_16"           "MAFA_146"           
# [34] "USF2_26"             "MLXIPL_15"           "E2F3_313"           
# [37] "CEBPB_140"           "MAX_14"              "ZNF76_164"          
# [40] "HES4_95"             "E2F7_316"            "ZNF148_222"         
# [43] "E2F5_315"            "KLF9_192"            "MYCL1_37"           
# [46] "SRY_771"             "E2F1_312"            "HLTF_802"           
# [49] "MEIS2_471"           "MXI1_40"             "TFCP2_392"          
# [52] "RFX6_729"            "SMAD4_739"           "MYPOP_652"          
# [55] "HOXC6_561"           "ZIC2_162"            "BARHL1_456"         
# [58] "MSX2_449"            "BSX_559"             "MYBL2_647"          
# [61] "MAZ_178"             "HOXD11_464"          "HOXA13_418"         
# [64] "NKX11_597"           "HMX2_558"            "TGIF1_541"          
# [67] "NOBOX_420"           "ENSG00000250542_156" "VSX1_408"           
# [70] "MZF1_171"            "PITX3_426"           "OTX1_437"           
# [73] "HMGA2_13"            "ESRRB_669"           "GCM2_389"           
# [76] "CIC_749"             "MITF_91"             "HOXC11_453"         
# > markerMotifList[["Dark_GC_B"]]$name
#  [1] "PITX2_504"   "HMGA1_12"    "DNMT1_301"   "ZNF524_243"  "CBFB_801"   
#  [6] "MBD2_644"    "HEYL_65"     "TBP_793"     "ARID3C_10"   "MAFB_150"   
# [11] "NR2E3_657"   "ZSCAN10_206" "IRF4_632"    "MAFG_148"    "HOXC10_548" 
# [16] "BCL11A_194"  "BCL11B_825"  "ARID5B_7"    "NR3C1_666"  