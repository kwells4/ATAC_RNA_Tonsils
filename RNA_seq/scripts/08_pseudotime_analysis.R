library(Rmagic)
library(Seurat)
library(slingshot)
library(RColorBrewer)
library(ArchR)
library(viridis)
library(gam)



save_path <- "output/allSamples_nomt_BCells_snakemake/"

function_path <- "scripts/functions.R"
ATAC_path <- "/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/ArchR_analysis_scg4/B_cells"

source(function_path)

seurat_object <- readRDS(paste0(save_path,
    "files/BCells_dubRm_namedClust_seurat.rda"))

projTonsils5 <- loadArchRProject(ATAC_path)

ATAC_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/ArchR_analysis_scg4/"

pseudotime_scores <- "slingshot"

sce.pca <- readRDS(paste0(save_path,
  "files/slingshot_object_harmony_cell_type_bCell.rds"))

cell_types <- levels(seurat_object$cell_type_bCell)

colors_cell_type_new <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type_new) <- cell_types

add_pseudotime <- FALSE

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

# Decide how to group pseudotime scoes (this groups by every two blocks of psudotime)
groupEvery <- 2
smoothWindow <- 7
labelMarkers <- c("BHLHE40", "PSME2", "APEX2", "MIR155HG", "CD69",
  "JUN", "FOS", "DUSP2", "NME2", "BATF")

labelMarkers_chr <- c("chr3:BHLHE40", "chr14:PSME2", "chrX:APEX2",
  "chr21:MIR155HG", "chr12:CD69", "chr1:JUN", "chr14:FOS", "chr2:DUSP2",
  "chr17:NME2", "chr14:BATF")

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("Starting analysis with  Slingshot")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


########################
# Slingshot pseudotime #
########################

##################### Add pseudotime data to seurat object #######################
seurat_pseudotime <- seurat_object

if(add_pseudotime){
  seurat_pseudotime$curve1 <- NULL
  pseudotime <- data.frame(slingshot::slingPseudotime(sce.pca))

  pseudotime$curve1_notScaled <- pseudotime$curve1
  pseudotime$curve2_notScaled <- pseudotime$curve2
  pseudotime$curve1 <- NULL
  pseudotime$curve2 <- NULL
    # scale to be between 0 and 100 (necessary for ArchR)
  actual_min <- min(pseudotime$curve1_notScaled[!is.na(pseudotime$curve1_notScaled)])
  actual_max <- max(pseudotime$curve1_notScaled[!is.na(pseudotime$curve1_notScaled)])
  desired_min <- 0
  desired_max <- 100
  pseudotime$curve1 <- ((pseudotime$curve1_notScaled - actual_min) / 
    (actual_max - actual_min)) * (desired_max - desired_min) + desired_min

  actual_min2 <- min(pseudotime$curve2_notScaled[!is.na(pseudotime$curve2_notScaled)])
  actual_max2 <- max(pseudotime$curve2_notScaled[!is.na(pseudotime$curve2_notScaled)])
  pseudotime$curve2 <- ((pseudotime$curve2_notScaled - actual_min2) / 
    (actual_max2 - actual_min2)) * (desired_max - desired_min) + desired_min

  seurat_meta <- seurat_pseudotime[[]]

  seurat_meta <- merge(seurat_meta, pseudotime, by = "row.names", all.x = TRUE)

  rownames(seurat_meta) <- seurat_meta$Row.names

  seurat_meta <- seurat_meta[match(colnames(seurat_pseudotime), rownames(seurat_meta)), ]
  identical(rownames(seurat_meta), colnames(seurat_pseudotime))

  seurat_pseudotime$curve1_notScaled <- seurat_meta$curve1_notScaled
  seurat_pseudotime$slingshot_curve1 <- seurat_meta$curve1

  seurat_pseudotime$curve2_notScaled <- seurat_meta$curve2_notScaled
  seurat_pseudotime$slingshot_curve2 <- seurat_meta$curve2

  seurat_pseudotime$in_trajectory <- "not_in_trajectory"

  seurat_pseudotime$in_trajectory[!is.na(seurat_pseudotime$slingshot_curve1)] <- "in_trajectory"

  seurat_pseudotime$in_trajectory_curve2 <- "not_in_trajectory"

  seurat_pseudotime$in_trajectory_curve2[!is.na(seurat_pseudotime$slingshot_curve2)] <- "in_trajectory"



  # Plot UMAP of RNA
  pseudotime_plot <- plotDimRed(seurat_pseudotime, col_by = "slingshot_curve1",
    highlight_group = TRUE,
    meta_data_col = "in_trajectory", group = "in_trajectory",
    ggrastr = TRUE)
  pseudotime_plot <- pseudotime_plot[[1]]
  # Add lines to plot
  # pull out coordinates
  plot_coord <- data.frame(Embeddings(object = seurat_pseudotime, reduction = "umap"))
  plot_coord$cell_type <- seurat_pseudotime$cell_type_bCell

  # Find centroid of clusters
  gg <- merge(plot_coord, aggregate(cbind(mean.x=UMAP_1, 
    mean.y=UMAP_2) ~cell_type, plot_coord, mean), by="cell_type")

  # remove UMAP columns
  gg <- gg[ , -c(2:3)]

  # Keep only one row for each centroid
  gg <- dplyr::distinct(gg)

  # Keep only clusters that are in the analysis
  gg <- gg[gg$cell_type %in% sce.pca@lineages$Lineage1, ]

  # Order based on pseudotime order
  gg$cell_type <- factor(gg$cell_type, levels = sce.pca@lineages$Lineage1)
  gg <- gg[order(gg$cell_type), ]

  pdf(paste0(save_path, "images/pseudotime/slingshot/pseudotime_umap.pdf"))

  print(pseudotime_plot + scale_color_viridis(option = "plasma") +
    ggplot2::ggtitle("slingshot pseudotime"))
  print(pseudotime_plot + scale_color_viridis(option = "plasma") +
    ggplot2::ggtitle("slingshotpseudotime") + 
    ggplot2::geom_path(data = gg, aes(mean.x, mean.y),
      size = 2, arrow = arrow(angle = 15, type = "closed")))
  dev.off()

  # get color scale:
  full_plot <- pseudotime_plot + scale_color_viridis(option = "plasma")
  pseudotime_data <- full_plot$data
  color_data <- ggplot_build(full_plot)
  color_data <- color_data$data[[1]]$colour
  pseudotime_data$color <- color_data
  saveRDS(pseudotime_data, paste0(save_path, "files/slingshot_files/color_mapping.rda"))
  
    # Plot UMAP of RNA
  pseudotime_plot2 <- plotDimRed(seurat_pseudotime, col_by = "slingshot_curve2",
    highlight_group = TRUE,
    meta_data_col = "in_trajectory_curve2", group = "in_trajectory",
    ggrastr = TRUE)
  pseudotime_plot2 <- pseudotime_plot2[[1]]
  # Add lines to plot
  # pull out coordinates
  plot_coord2 <- data.frame(Embeddings(object = seurat_pseudotime, reduction = "umap"))
  plot_coord2$cell_type <- seurat_pseudotime$cell_type_bCell

  # Find centroid of clusters
  gg2 <- merge(plot_coord2, aggregate(cbind(mean.x=UMAP_1, 
    mean.y=UMAP_2) ~cell_type, plot_coord2, mean), by="cell_type")

  # remove UMAP columns
  gg2 <- gg2[ , -c(2:3)]

  # Keep only one row for each centroid
  gg2 <- dplyr::distinct(gg2)

  # Keep only clusters that are in the analysis
  gg2 <- gg2[gg2$cell_type %in% sce.pca@lineages$Lineage2, ]

  # Order based on pseudotime order
  gg2$cell_type <- factor(gg$cell_type, levels = sce.pca@lineages$Lineage2)
  gg2 <- gg2[order(gg$cell_type), ]

  pdf(paste0(save_path, "images/pseudotime/slingshot/pseudotime_umap_curve-2.pdf"))

  print(pseudotime_plot2 + scale_color_viridis(option = "plasma") +
    ggplot2::ggtitle("slingshot pseudotime curve 2"))
  print(pseudotime_plot + scale_color_viridis(option = "plasma") +
    ggplot2::ggtitle("slingshot pseudotime") + 
    ggplot2::geom_path(data = gg2, aes(mean.x, mean.y),
      size = 2, arrow = arrow(angle = 15, type = "closed")))
  dev.off()

  #################### Add pseudotime info to ATAC object ###################
  ATAC_sample_data <- getCellColData(projTonsils5)

  # Get pseudotime
  pseudotime$RNA_cell <- rownames(pseudotime)
  ATAC_sample_data$RNA_cell <- ATAC_sample_data$predictedCell

  ATAC_sample_data$ATAC_cell <- rownames(ATAC_sample_data)

  ATAC_pseudotime <- merge(ATAC_sample_data, pseudotime, by = "RNA_cell",
    all.x = TRUE)
  rownames(ATAC_pseudotime) <- ATAC_pseudotime$ATAC_cell

  ATAC_pseudotime <- ATAC_pseudotime[order(match(rownames(ATAC_pseudotime),
    rownames(ATAC_sample_data))), ]

  if(identical(rownames(ATAC_pseudotime), rownames(projTonsils5@cellColData))){
    projTonsils5@cellColData[, "RNA_pseudotime1_not_scaled"] <- ATAC_pseudotime[["curve1_notScaled"]]
    projTonsils5@cellColData[, "RNA_pseudotime1"] <- ATAC_pseudotime[["curve1"]]
    projTonsils5@cellColData[, "RNA_pseudotime2"] <- ATAC_pseudotime[["curve2"]]
  } else {
    print("rownames don't match!!!")
  }

  # Make UMAP with ArchR

  p <- plotTrajectory(projTonsils5, trajectory = "RNA_pseudotime1",
      colorBy = "cellColData", name = "RNA_pseudotime1",
      embedding = "UMAP_Harmony")

  plotPDF(p, name = "Plot-RNA-Traj-UMAP.pdf", ArchRProj = projTonsils5,
      addDOC = FALSE, width = 5, height = 5)


  p <- plotTrajectory(projTonsils5, trajectory = "RNA_pseudotime2",
      colorBy = "cellColData", name = "RNA_pseudotime2",
      embedding = "UMAP_Harmony")

  plotPDF(p, name = "Plot-RNA-Traj-UMAP-Curve2.pdf", ArchRProj = projTonsils5,
      addDOC = FALSE, width = 5, height = 5)


  saveRDS(seurat_pseudotime, paste0(save_path,
      "files/BCells_dubRm_namedClust_seurat.rda"))

  # saveRDS(seurat_object, paste0(save_path,
  #     "files/BCells_dubRm_namedClust_seurat_old.rda"))

  saveArchRProject(ArchRProj = projTonsils5, outputDirectory = ATAC_path)
   
}



######################### Run Magic on RNA ##########################
# print("!!!!!!!!!!!!!!!!!!")
# print("pseudotime of ADTs")
# print("!!!!!!!!!!!!!!!!!!")
seurat_stable <- seurat_pseudotime

# adt_list <- c("IgM-", "IgD-", "CD44-", "CD10-", "CD20-", "CD184-(CXCR4)-",
#               "CD185-(CXCR5)-",  "CD38-",  "CD27-", "CD3-", "CD8a-", "CD4-")

# adt_rna_list <- c("IGHM", "IGHD", "CD44", "MME", "MS4A1", "CXCR4", "CXCR5",
#                   "CD38",  "CD27", "CD3G", "CD8A", "CD4")

# # IGM and IGD weren't in the ATAC lists
# adt_list <- c("CD44-", "CD10-", "CD20-", "CD184-(CXCR4)-",
#               "CD185-(CXCR5)-",  "CD38-",  "CD27-", "CD3-", "CD8a-", "CD4-")

# adt_rna_list <- c("CD44", "MME", "MS4A1", "CXCR4", "CXCR5",
#                   "CD38",  "CD27", "CD3G", "CD8A", "CD4")
# # Run magic
# seurat_pseudotime <- magic(seurat_stable, genes = adt_rna_list)

# DefaultAssay(seurat_pseudotime) <- "ADT"
# seurat_pseudotime <- magic(seurat_pseudotime)

# DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"


########################## Pseudotime for ADTs ####################
# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_gene_expression.pdf"),
#   height = 9, width = 7)
# plot_list <- lapply(1:length(adt_rna_list), function(x){
#   print(x)
#   print(adt_rna_list[x])
#   print(adt_list[x])

#   # Plot RNA
#   DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"
#   pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#       y_val = adt_rna_list[[x]],
#       col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#       color = "viridis", ggrastr = FALSE)
#   pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#   # Plot ADT
#   DefaultAssay(seurat_pseudotime) <- "MAGIC_ADT"
#   pseudotime_plot_ADT <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#       y_val = adt_list[[x]],
#       col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#       color = "viridis", ggrastr = FALSE)
#   pseudotime_plot_ADT <- pseudotime_plot_ADT + ggplot2::ggtitle("ADT expression")

#   # Plot ATAC
#   pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#       archr_matrix = "GeneScoreMatrix", gene = adt_rna_list[[x]],
#       trajectory = "RNA_pseudotime1", color = "viridis",
#       col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#       impute = TRUE, ggrastr = FALSE)
#   pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#   # make grid of three plots object
#   gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#       pseudotime_plot_ADT, pseudotime_plot_ATAC), nrow = 3)
#   })

# dev.off()


########################## Pseudotime for genes of interest ####################
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("pseudotime of genes of interest")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

# Run magic

# seurat_pseudotime <- magic(seurat_stable, genes = labelMarkers)
# DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"

# motifs <- labelMarkers
# markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
#   useMatrix = "MotifMatrix")

# markerMotifs <- grep("z:", markerMotifs, value = TRUE)

# markerGenes <- sub("z:", "", markerMotifs)
# markerGenes <- sub("_.*", "", markerGenes)
# labelMarkers_short <- labelMarkers[labelMarkers %in% markerGenes]

# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_Hamish_genes.pdf"),
#   height = 9, width = 7)
# plot_list <- lapply(1:length(labelMarkers_short), function(x){
#   print(x)
#   print(labelMarkers_short[x])
#   # Figure this out...
#   gene_motif <- markerMotifs[grepl(paste0("z:", labelMarkers_short[x], "_"),
#     markerMotifs)]
#   print(gene_motif)


#   # Plot RNA
#   pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#     y_val = labelMarkers_short[[x]],
#     col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#     color = "viridis")
#   pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#   # Plot ATAC GeneScoreMatrix
#   pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#     archr_matrix = "GeneScoreMatrix", gene = labelMarkers_short[[x]],
#     trajectory = "RNA_pseudotime1", color = "viridis",
#     col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#     impute = TRUE)
#   pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#   # Plot ATAC Motif Matrix
#   pseudotime_plot_ATAC_dev <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#     archr_matrix = "MotifMatrix", gene = gene_motif,
#     trajectory = "RNA_pseudotime1", color = "viridis",
#     col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#     impute = TRUE, useSeqnames = "z")
#   pseudotime_plot_ATAC_dev <- pseudotime_plot_ATAC_dev + ggplot2::ggtitle("ATAC deviation z scores")


#   # make grid of three plots object
#   gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#     pseudotime_plot_ATAC, pseudotime_plot_ATAC_dev), nrow = 3)
#   })

# dev.off()

# labelMarkers_not_tf <- labelMarkers[!(labelMarkers %in% markerGenes)]

# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_Hamish_genes_not_tf.pdf"),
#   height = 6, width = 7)
# plot_list <- lapply(1:length(labelMarkers_not_tf), function(x){
#   print(x)
#   print(labelMarkers_not_tf[x])
#   # Figure this out...

#   # Plot RNA
#   pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#     y_val = labelMarkers_not_tf[[x]],
#     col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#     color = "viridis")
#   pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#   # Plot ATAC GeneScoreMatrix
#   pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#     archr_matrix = "GeneScoreMatrix", gene = labelMarkers_not_tf[[x]],
#     trajectory = "RNA_pseudotime1", color = "viridis",
#     col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#     impute = TRUE)
#   pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#   # make grid of three plots object
#   gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#     pseudotime_plot_ATAC), nrow = 2)
#   })

# dev.off()


################### Pseudotime for txn regulators ################
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("pseudotime on txn regualtors")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

# Motif matrix correlation
txn_regulators_MM <- c("BATF", "BCL11A", "CEBPA-DT", "CEBPB", "CEBPD", "EBF1",
                       "EOMES", "ERF", "ETS1-AS1", "FOSL2", "IRF1", "IRF4",
                       "KLF11", "LEF1", "LMO2", "MEF2A", "MEF2B", "MEF2C",
                       "MYF6", "NFIA", "NFYA", "PAX5", "POU2F1", "POU2F2",
                       "POU3F3", "RUNX2", "SP4", "SPI1", "SPIB", "STAT1",
                       "TAL2", "TCF15", "TWIST1")  

# Integration matrix correlation
txn_regulators_IM <- c("BATF", "BCL11A", "CEBPA", "CEBPB", "CEBPD", "CTCFL", "EBF1",
                       "ELK1", "EOMES", "ERF", "HMGA1", "IRF2", "IRF9", "KLF2",
                       "LEF1", "LMO2", "MEF2B", "MEF2C", "NFE2", "NFIA", "NFIB",
                       "PAX5", "PKNOX1", "POU2F1", "POU2F2", "SPI1", "SPIB",
                       "STAT1", "STAT2", "TEAD1") 



# All
txn_regulators <- unique(txn_regulators_IM, txn_regulators_MM)


# seurat_pseudotime <- magic(seurat_stable, genes = txn_regulators)
# DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"

# motifs <- txn_regulators
# markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
#   useMatrix = "MotifMatrix")

# markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# # remove markers that weren't correctly grabbed
# # markerMotifs <- markerMotifs[!grepl("CTCFL", markerMotifs)]
# # markerMotifs <- markerMotifs[!grepl("SREBF1", markerMotifs)]
# # markerMotifs <- markerMotifs[!grepl("PRDM16", markerMotifs)]

# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_tf_regulator_expression.pdf"),
#   height = 9, width = 7)
# plot_list <- lapply(1:length(txn_regulators), function(x){
#   print(x)
#   print(txn_regulators[x])
#   # Figure this out...
#   gene_motif <- markerMotifs[grepl(paste0("z:", txn_regulators[x], "_"),
#     markerMotifs)]
#   print(gene_motif)


#   # Plot RNA
#   pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#     y_val = txn_regulators[[x]],
#     col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#     color = "viridis")
#   pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#   # Plot ATAC GeneScoreMatrix
#   pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#     archr_matrix = "GeneScoreMatrix", gene = txn_regulators[[x]],
#     trajectory = "RNA_pseudotime1", color = "viridis",
#     col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#     impute = TRUE)
#   pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#   # Plot ATAC Motif Matrix
#   pseudotime_plot_ATAC_dev <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#     archr_matrix = "MotifMatrix", gene = gene_motif,
#     trajectory = "RNA_pseudotime1", color = "viridis",
#     col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#     impute = TRUE, useSeqnames = "z")
#   pseudotime_plot_ATAC_dev <- pseudotime_plot_ATAC_dev + ggplot2::ggtitle("ATAC deviation z scores")


#   # make grid of three plots object
#   gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#     pseudotime_plot_ATAC, pseudotime_plot_ATAC_dev), nrow = 3)
#   })

# dev.off()


############### Txn regulators from B cells alone ##################
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("pseudotime on txn regulators from B cells")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


# seGroupMotif <- getGroupSE(ArchRProj = projTonsils5,
#   useMatrix = "MotifMatrix", groupBy = "Clusters2")

# seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

# rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
#   rowMaxs(assay(seZ) - assay(seZ)[,x])
# }) %>% Reduce("cbind", .) %>% rowMaxs

# corGSM_MM <- correlateMatrices(
#     ArchRProj = projTonsils5,
#     useMatrix1 = "GeneScoreMatrix",
#     useMatrix2 = "MotifMatrix",
#     reducedDims = "Harmony"
# )

# corGSM_MM

# corGIM_MM <- correlateMatrices(
#     ArchRProj = projTonsils5,
#     useMatrix1 = "GeneIntegrationMatrix",
#     useMatrix2 = "MotifMatrix",
#     reducedDims = "Harmony"
# )

# corGIM_MM

# corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name,
#   rowData(seZ)$name), "maxDelta"]
# corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name,
#   rowData(seZ)$name), "maxDelta"]

# corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
# corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
# corGSM_MM$TFRegulator <- "NO"
# corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
# tfReg_GSM_MM <- sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


# corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
# corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
# corGIM_MM$TFRegulator <- "NO"
# corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
# tfReg_GIM_MM <- sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

# txn_regulators_bcell <- unique(tfReg_GSM_MM, tfReg_GIM_MM)

# txn_regulators_bcell <- txn_regulators_bcell[txn_regulators_bcell %in% 
#   rownames(seurat_pseudotime)]


# seurat_pseudotime <- magic(seurat_stable, genes = txn_regulators_bcell)
# DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"

# motifs <- txn_regulators_bcell
# markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
#   useMatrix = "MotifMatrix")

# markerMotifs <- grep("z:", markerMotifs, value = TRUE)


# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_b_cell_tf_regulator_expression.pdf"),
#   height = 9, width = 7)
# plot_list <- lapply(1:length(txn_regulators_bcell), function(x){
#   print(x)
#   print(txn_regulators_bcell[x])
#   # Figure this out...
#   gene_motif <- markerMotifs[grepl(paste0("z:", txn_regulators_bcell[x], "_")
#     , markerMotifs)]
#   print(gene_motif)
#   if(length(gene_motif) > 0){

#     # Plot RNA
#     pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#       y_val = txn_regulators_bcell[[x]],
#       col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#       color = "viridis")
#     pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#     # Plot ATAC GeneScoreMatrix
#     pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#       archr_matrix = "GeneScoreMatrix", gene = txn_regulators_bcell[[x]],
#       trajectory = "RNA_pseudotime1", color = "viridis",
#       col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#       impute = TRUE)
#     pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#     # Plot ATAC Motif Matrix
#     pseudotime_plot_ATAC_dev <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#       archr_matrix = "MotifMatrix", gene = gene_motif,
#       trajectory = "RNA_pseudotime1", color = "viridis",
#       col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#       impute = TRUE, useSeqnames = "z")
#     pseudotime_plot_ATAC_dev <- pseudotime_plot_ATAC_dev + ggplot2::ggtitle("ATAC deviation z scores")


#     # make grid of three plots object
#     gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#       pseudotime_plot_ATAC, pseudotime_plot_ATAC_dev), nrow = 3)
#   }
# })

# dev.off()

############### Heatmaps of genes that correlate with pseudotime ################
# Keep working to see if it will work with Jeff's pipeline.
# ATAC

# trajectory <- c("Naive_B", "Activated_B", "Light_GC_B", "Plasmablasts")

# projTonsils5 <- addTrajectory(
#     ArchRProj = projTonsils5, 
#     name = "ATAC_trajectory", 
#     groupBy = "Clusters2",
#     trajectory = trajectory, 
#     embedding = "UMAP_Harmony", 
#     force = TRUE
# )
print("!!!!!!!!!!!!!!!!!!!")
print("pseudotime heatmaps")
print("!!!!!!!!!!!!!!!!!!!")

###########################

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("Pseudotime based on even splits")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

# Make a trajectory based on even splits
seurat_pseudotime <- seurat_stable

ATAC_sample_data <- getCellColData(projTonsils5)

# Get pseudotime
# ATAC_sample_data$RNA_cells <- ATAC_sample_data$predictedCell

# ATAC_sample_data$ATAC_cells <- rownames(ATAC_sample_data)

# ATAC_pseudotime <- merge(ATAC_sample_data, pseudotime, by = "RNA_cell",
#   all.x = TRUE)


trajectory_ATAC <- getCellColData(projTonsils5, c("RNA_pseudotime1",
  "predictedCell"))



trajectory_ATAC <- trajectory_ATAC[!is.na(trajectory_ATAC[,1]),,drop=FALSE]

# trajectory_ATAC$pseudotime_merge <- trajectory_ATAC$RNA_pseudotime1

trajectory_ATAC$RNA_cells <- trajectory_ATAC$predictedCell

trajectory_RNA <- seurat_pseudotime[["slingshot_curve1"]]

trajectory_RNA <- trajectory_RNA[!is.na(trajectory_RNA[,1]), , drop=FALSE]
# trajectory_RNA$pseudotime_merge <- trajectory_RNA$slingshot_curve1
trajectory_RNA$RNA_cells <- rownames(trajectory_RNA)
trajectory_ATAC$ATAC_cells <- rownames(trajectory_ATAC)

# trajectory_all <- merge(trajectory_ATAC, trajectory_RNA, by = "pseudotime_merge",
#   all.x = TRUE, all.y = TRUE)

trajectory_all <- merge(trajectory_ATAC, trajectory_RNA, by = "RNA_cells",
  all.x = TRUE, all.y = TRUE)

trajectory_all <- trajectory_all[order(trajectory_all$slingshot_curve1), ]

# Find a split that leaves at least ncells for both RNA and ATAC
ncells <- 100
line_number <- 1
new_group <- TRUE
cell_group <- 0
groupListATAC <- list()
groupListRNA <- list()
pseudotimeList <- list()
while(line_number <= nrow(trajectory_all)){
  trajectory_sub <- trajectory_all[line_number, ]
  if (new_group){
    RNA_cells <- c()
    ATAC_cells <- c()
    pseudotime_start <- trajectory_sub$slingshot_curve1
    cell_group <- cell_group + 1
  }
  # Only add if it exists in the trajectory and is not already present in the list
  if (!is.na(trajectory_sub$ATAC_cells)){
    if(length(groupListATAC) > 0){
      if(!trajectory_sub$ATAC_cells %in% unlist(groupListATAC)){
        ATAC_cells <- c(ATAC_cells, trajectory_sub$ATAC_cells)
      }
    } else {
      ATAC_cells <- c(ATAC_cells, trajectory_sub$ATAC_cells)
    } 
  }
  if (!is.na(trajectory_sub$RNA_cells)){
    if(length(groupListRNA) > 0){
      if(!trajectory_sub$RNA_cells %in% unlist(groupListRNA)){
        RNA_cells <- c(RNA_cells, trajectory_sub$RNA_cells)
      }
    } else {
      RNA_cells <- c(RNA_cells, trajectory_sub$RNA_cells)
    }  
  }
  if (length(unique(RNA_cells)) >= ncells & length(unique(ATAC_cells)) >= ncells) {
    pseudotime_end <- trajectory_sub$slingshot_curve1
    print(line_number)
    pseudotimeList[[cell_group]] <- c(pseudotime_start, pseudotime_end)
    groupListRNA[[cell_group]] <- unique(RNA_cells)
    groupListATAC[[cell_group]] <- unique(ATAC_cells)
    new_group <- TRUE
  } else {
    new_group <- FALSE
  }
  line_number <- line_number + 1
}

if(!new_group){
  if(length(unique(RNA_cells)) > ncells/2 & length(unique(ATAC_cells)) > ncells/2){
    pseudotimeList[[cell_group]] <- c(pseudotime_start, trajectory_sub$slingshot_curve1)
    groupListRNA[[cell_group]] <- unique(RNA_cells)
    groupListATAC[[cell_group]] <- unique(ATAC_cells)
  } else {
    pseudotimeList[[cell_group - 1]][2] <- trajectory_sub$slingshot_curve1
    groupListRNA[[cell_group - 1]] <- c(groupListRNA[[cell_group - 1]],
      unique(RNA_cells))
    groupListATAC[[cell_group - 1]] <- c(groupListATAC[[cell_group - 1]],
      unique(ATAC_cells))
  }
}

if(!length(unlist(groupListRNA)) == nrow(trajectory_RNA)){
  stop("Cell lists don't match!")
}

if(!length(unlist(groupListATAC)) == nrow(trajectory_ATAC)){
  stop("Cell lists don't match!")
}
# # refine list to include all cells at the same point in pseudotime
# pseudotimeList <- lapply(seq_along(pseudotimeList), function(x){
#   if(x == 1){
#     return(pseudotimeList[[x]])
#   } else {
#     return(c(pseudotimeList[[x-1]][2], pseudotimeList[[x]][2]))
#   }
#   })

groupNames <- lapply(seq_along(pseudotimeList), function(x){
  return(paste0("T.", pseudotimeList[[x]][1], "_", pseudotimeList[[x]][2]))
  })

# groupListATAC <- lapply(seq_along(pseudotimeList), function(x){
#   if(x == 1){
#     cells <- rownames(trajectory_ATAC)[which(trajectory_ATAC[,1] >= 
#       pseudotimeList[[x]][1] &
#       trajectory_ATAC[,1] <= pseudotimeList[[x]][2])]
#   } else {
#     cells <- rownames(trajectory_ATAC)[which(trajectory_ATAC[,1] >
#       pseudotimeList[[x]][1] &
#       trajectory_ATAC[,1] <= pseudotimeList[[x]][2])]
#   }
#   print(length(cells))
#   return(cells)
# })

names(groupListATAC) <- groupNames

# Find same breaks in RNA
# groupListRNA <- lapply(seq_along(pseudotimeList), function(x){
#   if(x == 1){
#     cells <- rownames(trajectory_RNA)[which(trajectory_RNA[,1] >=
#       pseudotimeList[[x]][1] &
#       trajectory_RNA[,1] <= pseudotimeList[[x]][2])]
#   } else {
#     cells <- rownames(trajectory_RNA)[which(trajectory_RNA[,1] >
#       pseudotimeList[[x]][1] &
#       trajectory_RNA[,1] <= pseudotimeList[[x]][2])]
#   }
#   print(length(cells))
#   return(cells)
# })
names(groupListRNA) <- groupNames

##################################
# repeat trajectories using these group lists
trajMM_new  <- getTrajectory_Kristen(ArchRProj = projTonsils5,
  name = "RNA_pseudotime1", useMatrix = "MotifMatrix",
  log2Norm = FALSE, groupEvery = groupEvery,
  smoothWindow = smoothWindow, groupList = groupListATAC)
# This failed because there were some time points that there are no ATAC cells
# groupEvery = 2 fixes this because more point of pseudotime are combined
# If this smooths it too much, I can rewrite the function to remove those
# points of pseudotime
trajGSM_new <- getTrajectory_Kristen(ArchRProj = projTonsils5,
  name = "RNA_pseudotime1", useMatrix = "GeneScoreMatrix",
  log2Norm = TRUE, groupEvery = groupEvery, smoothWindow = smoothWindow,
  groupList = groupListATAC)

# Gene integration matrix
trajGIM_new <- getTrajectory_Kristen(ArchRProj = projTonsils5,
  name = "RNA_pseudotime1", useMatrix = "GeneIntegrationMatrix",
  log2Norm = FALSE, groupEvery = groupEvery, smoothWindow = smoothWindow,
  groupList = groupListATAC)

# Peak matrix
trajPM_new  <- getTrajectory_Kristen(ArchRProj = projTonsils5,
  name = "RNA_pseudotime1", useMatrix = "PeakMatrix",
  log2Norm = TRUE, groupEvery = groupEvery, smoothWindow = smoothWindow,
  groupList = groupListATAC)

# RNA (I wrote my own function to be modeled after Jeffs, can now be plotted
# with his functions)
trajRNA_new <- getTrajectoryRNA(seurat_object = seurat_pseudotime, name = "slingshot_curve1",
  log2Norm = TRUE, groupEvery = groupEvery, smoothWindow = smoothWindow,
  groupList = groupListRNA)



p1_new <- plotTrajectoryHeatmap(trajMM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, useSeqnames = "z", labelMarkers = labelMarkers_chr)
mat1_new <- plotTrajectoryHeatmap(trajMM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, useSeqnames = "z", returnMat = TRUE)

# Gene score matrix
p2_new <- plotTrajectoryHeatmap(trajGSM_new,
  pal = paletteContinuous(set = "horizonExtra"),
  labelTop = 50, labelMarkers = labelMarkers_chr)
mat2_new <- plotTrajectoryHeatmap(trajGSM_new,
  pal = paletteContinuous(set = "horizonExtra"),
  labelTop = 50, returnMat = TRUE)

# Gene integration matrix
p3_new <- plotTrajectoryHeatmap(trajGIM_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, labelMarkers = labelMarkers_chr)
mat3_new <- plotTrajectoryHeatmap(trajGIM_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, returnMat = TRUE)

# Peak matrix
p4_new <- plotTrajectoryHeatmap(trajPM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50)
mat4_new <- plotTrajectoryHeatmap(trajPM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, returnMat = TRUE)

# RNA (I wrote my own function to be modeled after Jeffs, can now be plotted
  # with his functions)
pRNA_new <- plotTrajectoryHeatmap(trajRNA_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, labelMarkers = labelMarkers)
matRNA_new <- plotTrajectoryHeatmap(trajRNA_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, returnMat = TRUE)

# Also get peak matrix based on nearest gene:
output_dir <- getOutputDirectory(projTonsils5)
peak_set <- readRDS(paste0(output_dir, "/files/macs2_peaks.rds"))
peak_set$comb_names <- paste0(seqnames(peak_set), ":", ranges(peak_set))
peak_set$comb_names <- sub("-", "_", peak_set$comb_names)

peak_set_ord_new <- peak_set[match(rownames(rowData(trajPM_new)), peak_set$comb_names), ]
identical(peak_set_ord_new$comb_names, rownames(rowData(trajPM_new)))
identical(peak_set_ord_new$comb_names, rownames(assays(trajPM_new)$smoothMat))

trajPM_nearest_gene_new <- trajPM_new

smoothMat_new <- assays(trajPM_nearest_gene_new)$smoothMat

rownames(smoothMat_new) <- paste0(rownames(smoothMat_new), ":",
  peak_set_ord_new$nearestGene)
rownames(trajPM_nearest_gene_new) <- rownames(smoothMat_new)

p5_new <- plotTrajectoryHeatmap(trajPM_nearest_gene_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50)
mat5_new <- plotTrajectoryHeatmap(trajPM_nearest_gene_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, returnMat = TRUE)

write.csv(rownames(mat1_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/MM_gene_order_smooth100cells.csv"))
write.csv(rownames(mat2_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/GSM_gene_order_smooth100cells.csv"))
write.csv(rownames(mat3_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/GIM_gene_order_smooth100cells.csv"))
write.csv(rownames(mat4_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/PM_gene_order_smooth100cells.csv"))
write.csv(rownames(matRNA_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/RNA_gene_order_smooth100cells.csv"))
write.csv(rownames(mat5_new), paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/PM_gene_order_nearest_gene_smooth100cells.csv"))



pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/pseudotime_heatmaps_smooth100cells.pdf"),
  height = 9, width = 7)
p1_new

p2_new

p3_new

p4_new

p5_new

pRNA_new

dev.off()

# ATAC
trajectory_ATAC <- getCellColData(projTonsils5, c("RNA_pseudotime1", "Clusters2"))
trajectory_ATAC <- trajectory_ATAC[!is.na(trajectory_ATAC[,1]),,drop=FALSE]
trajectory_ATAC$Clusters2[trajectory_ATAC$Clusters2 == "Circulating_Dark_GC_B"] <-
  "GC_B"
trajectory_ATAC$Clusters2[trajectory_ATAC$Clusters2 == "Light_GC_B"] <-
  "GC_B"

max_ident_ATAC <- lapply(seq_along(groupListATAC), function(x){
  cells <- groupListATAC[[x]]
  trajectory_set <- trajectory_ATAC[cells, ]
  if(nrow(trajectory_set) > 0){
    return(names(which.max(table(trajectory_set$Clusters2))))
  } else {
    return("not_mapped")
  }
})

names(max_ident_ATAC) <- names(groupListATAC)

max_ident_ATAC <- unlist(max_ident_ATAC)

max_ident_ATAC <- max_ident_ATAC[max_ident_ATAC != "not_mapped"]
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# Fix the numbers here when you have time.
# Make DF of cell percent for each
trajectory_ATAC$Clusters2 <- factor(trajectory_ATAC$Clusters2,
  levels = levels(seurat_pseudotime$cell_type_bCell))
cell_types_ATAC <- lapply(seq_along(groupListATAC), function(x){
  cells <- groupListATAC[[x]]
  trajectory_set <- trajectory_ATAC[cells, ]
  if(nrow(trajectory_set) > 0){
    df <- data.frame(table(trajectory_set$Clusters2))
    names(df) <- c("cell_type", "frequency")
    df$percent <- df$frequency/sum(df$frequency) * 100
    df$time <- x*2
    return(df)
  } else {
    return(NULL)
  }
})

cell_types_ATAC <- do.call(rbind, cell_types_ATAC)


# RNA
trajectory_RNA <- seurat_pseudotime[[c("slingshot_curve1", "cell_type_bCell")]]
trajectory_RNA <- trajectory_RNA[!is.na(trajectory_RNA[,1]),,drop=FALSE]

max_ident_RNA <- lapply(seq_along(groupListRNA), function(x){
  cells <- groupListRNA[[x]]
  trajectory_set <- trajectory_RNA[cells, ]
  if(nrow(trajectory_set) > 0){
    return(names(which.max(table(trajectory_set$cell_type_bCell))))
  } else {
    return("not_mapped")
  }
})

names(max_ident_RNA) <- names(groupListRNA)

max_ident_RNA <- unlist(max_ident_RNA)


# Make DF of cell percent for each
print("!!!!!!!!!!!!!!!!")
print("percent of cells")
print("!!!!!!!!!!!!!!!!")

trajectory_RNA$cell_type_bCell <- factor(trajectory_RNA$cell_type_bCell,
  levels = levels(seurat_pseudotime$cell_type_bCell))
cell_types_RNA <- lapply(seq_along(groupListRNA), function(x){
  cells <- groupListRNA[[x]]
  trajectory_set <- trajectory_RNA[cells, ]
  if(nrow(trajectory_set) > 0){
    df <- data.frame(table(trajectory_set$cell_type_bCell))
    names(df) <- c("cell_type", "frequency")
    df$percent <- df$frequency/sum(df$frequency) * 100
    df$time <- x*2
    return(df)
  } else {
    return(NULL)
  }
})

cell_types_RNA <- do.call(rbind, cell_types_RNA)

colors_cell_type_new["GC_B"] <- "#E75480"

# Make a plot with both
pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/population_percents_smooth100cells.pdf"),
  height = 4, width = 8)

ggplot2::ggplot(cell_types_RNA, aes(x = time, y = percent, fill = cell_type)) +
  ggplot2::geom_area() +
  ggplot2::scale_fill_manual(values = colors_cell_type_new) +
  ggplot2::xlab("pseudotime") +
  ggplot2::ylab("Percent of population") +
  ggplot2::ggtitle("RNA cell types over pseudotime")


ggplot2::ggplot(cell_types_ATAC, aes(x = time, y = percent, fill = cell_type)) +
  ggplot2::geom_area() +
  ggplot2::scale_fill_manual(values = colors_cell_type_new) +
  ggplot2::xlab("pseudotime") +
  ggplot2::ylab("Percent of population") +
  ggplot2::ggtitle("ATAC cell types over pseudotime")

dev.off()


colors_cell_type_new["not_mapped"] <- "#DFDFDF"
colorMap <- list(cell_type = colors_cell_type_new)
attr(colorMap[[1]], "discrete") <- TRUE


# # Remake heatmaps with the cell type labeled
# p1 <- draw_heatmap(max_ident = max_ident_ATAC, mat = mat1,
#   pal = paletteContinuous(set = "solarExtra"), seTrajectory = trajMM)

# # Gene score matrix
# p2 <- draw_heatmap(max_ident = max_ident_ATAC, mat = mat2,
#   pal = paletteContinuous(set = "horizonExtra"), seTrajectory = trajGSM)

# # Gene integration matrix
# p3 <- draw_heatmap(max_ident = max_ident_ATAC, mat = mat3,
#   pal = paletteContinuous(set = "blueYellow"), seTrajectory = trajGIM)

# # Peak matrix
# p4 <- draw_heatmap(max_ident = max_ident_ATAC, mat = mat4,
#   pal = paletteContinuous(set = "solarExtra"), seTrajectory = trajPM)

# # RNA (I wrote my own function to be modeled after Jeffs, can now be plotted
#   # with his functions)
# pRNA <- draw_heatmap(max_ident = max_ident_RNA, mat = matRNA,
#   pal = paletteContinuous(set = "blueYellow"), seTrajectory = trajRNA)


p1_names_new <- plotTrajectoryHeatmap_kristen(trajMM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, useSeqnames = "z", max_ident= max_ident_ATAC,
  colorMap = colorMap, labelMarkers = labelMarkers_chr)


# Gene score matrix
p2_names_new <- plotTrajectoryHeatmap_kristen(trajGSM_new,
  pal = paletteContinuous(set = "horizonExtra"),
  labelTop = 50, max_ident= max_ident_ATAC, colorMap = colorMap,
  labelMarkers = labelMarkers_chr)

# Gene integration matrix
p3_names_new <- plotTrajectoryHeatmap_kristen(trajGIM_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, max_ident= max_ident_ATAC, colorMap = colorMap,
  labelMarkers = labelMarkers_chr)

# Peak matrix
p4_names_new <- plotTrajectoryHeatmap_kristen(trajPM_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, max_ident= max_ident_ATAC, colorMap = colorMap)

p5_names_new <- plotTrajectoryHeatmap_kristen(trajPM_nearest_gene_new,
  pal = paletteContinuous(set = "solarExtra"),
  labelTop = 50, max_ident= max_ident_ATAC, colorMap = colorMap)

# RNA (I wrote my own function to be modeled after Jeffs, can now be plotted
  # with his functions)
pRNA_names_new <- plotTrajectoryHeatmap_kristen(trajRNA_new,
  pal = paletteContinuous(set = "blueYellow"),
  labelTop = 50, max_ident= max_ident_RNA, colorMap = colorMap,
  labelMarkers = labelMarkers)

pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/pseudotime_heatmaps_cell_names_smooth100cells.pdf"),
  height = 9, width = 12)
p1_names_new

p2_names_new

p3_names_new

p4_names_new

p5_names_new

pRNA_names_new

dev.off()

print("!!!!!!!!!!!!!!!")
print("correlate genes")
print("!!!!!!!!!!!!!!!")

# print(trajGSM_new)
# print(trajMM_new)
# # Find accessibility that correlates with gene score
# corGSM_MM_new <- correlateTrajectories(trajGSM_new, trajMM_new,
#   varCutOff1 = 0.75, varCutOff2 = 0.75, corCutOff = 0.25)
# print("correlation")
# print(nrow(corGSM_MM_new[[1]]) > 1)
# # Subset to only the features that correlated
# if(nrow(corGSM_MM_new[[1]]) > 1){
#   trajGSM2_sub <- trajGSM_new[corGSM_MM_new[[1]]$name1, ]
#   trajMM2_sub <- trajMM_new[corGSM_MM_new[[1]]$name2, ]

#   # Order the heatmaps
#   trajCombined_new <- trajGSM2_sub
#   assay(trajCombined_new) <- t(apply(assay(trajGSM2_sub), 1, scale)) +
#     t(apply(assay(trajMM2_sub), 1, scale))

#   combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new, returnMat = TRUE,
#     varCutOff = 0)

#   rowOrder <- match(rownames(combinedMat_new), rownames(trajGSM2_sub))

#   ht1 <- plotTrajectoryHeatmap_kristen(trajGSM2_sub, 
#     pal = paletteContinuous(set = "horizonExtra"), 
#     varCutOff = 0, rowOrder = rowOrder,
#     max_ident= max_ident_ATAC, colorMap = colorMap)
#   ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
#     pal = paletteContinuous(set = "solarExtra"),
#     varCutOff = 0, rowOrder = rowOrder,
#     max_ident= max_ident_ATAC, colorMap = colorMap)
#   pdf(paste0(save_path, "images/pseudotime/",
#     pseudotime_scores, "/ATAC_correlate_GSM_MM.pdf"), height = 11, width = 16)
#   print(ht1 + ht2)
#   dev.off()

#   # repeate with RNA heatmap GSM genes
#   corGSM_MM_RNA_new <- corGSM_MM_new[[1]][corGSM_MM_new[[1]]$matchname1 %in% 
#     rownames(trajRNA_new),]
#   trajRNA2_sub <- trajRNA_new[corGSM_MM_RNA_new$matchname1, ]
#   trajMM2_sub <- trajMM_new[corGSM_MM_RNA_new$name2, ]


#   # Order the heatmaps
#   trajCombined_new <- trajRNA2_sub
#   assay(trajCombined_new) <- t(apply(assay(trajRNA2_sub), 1, scale)) +
#     t(apply(assay(trajMM2_sub), 1, scale))

#   combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new, returnMat = TRUE,
#     varCutOff = 0)

#   rowOrder <- match(rownames(combinedMat_new), rownames(trajRNA2_sub))

#   ht1 <- plotTrajectoryHeatmap_kristen(trajRNA2_sub,
#     pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0,
#     rowOrder = rowOrder, max_ident= max_ident_RNA, colorMap = colorMap)
#   ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
#     pal = paletteContinuous(set = "solarExtra"),
#     varCutOff = 0, rowOrder = rowOrder,
#     max_ident= max_ident_ATAC, colorMap = colorMap)

#   pdf(paste0(save_path, "images/pseudotime/",
#     pseudotime_scores, "/ATAC_correlate_GSM_MM_with_RNA_smooth100cells.pdf"),
#   height = 11, width = 16)
#   print(ht1 + ht2)
#   dev.off()

# }


# Repeat with gene integration
var_cutoff <- 0.7
corGIM_MM_new <- correlateTrajectories(trajGIM_new, trajMM_new,
  varCutOff1 = var_cutoff, varCutOff2 = var_cutoff, corCutOff = 0.25)
print("correlation")
print(nrow(corGIM_MM_new[[1]]) > 1)
if(nrow(corGIM_MM_new[[1]]) > 1) {

  trajGIM2_sub <- trajGIM_new[corGIM_MM_new[[1]]$name1, ]
  trajMM2_sub <- trajMM_new[corGIM_MM_new[[1]]$name2, ]

  trajCombined_new <- trajGIM2_sub
  assay(trajCombined_new) <- t(apply(assay(trajGIM2_sub), 1, scale)) + 
    t(apply(assay(trajMM2_sub), 1, scale))

  combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new,
    returnMat = TRUE, varCutOff = 0)

  rowOrder <- match(rownames(combinedMat_new), rownames(trajGIM2_sub))
  ht1 <- plotTrajectoryHeatmap_kristen(trajGIM2_sub, 
    pal = paletteContinuous(set = "horizonExtra"), 
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap)
  ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap)

  pdf(paste0(save_path, "images/pseudotime/",
    pseudotime_scores, "/ATAC_correlate_GIM_MM_smooth100cells_var",
    var_cutoff, ".pdf"), height = 11, width = 16)
  print(ht1 + ht2)
  dev.off()



  # repeat with RNA heatmap GIM genes
  corGIM_MM_RNA_new <- corGIM_MM_new[[1]][corGIM_MM_new[[1]]$matchname1 %in% rownames(trajRNA_new),]
  trajRNA2_sub <- trajRNA_new[corGIM_MM_RNA_new$matchname1, ]
  trajMM2_sub <- trajMM_new[corGIM_MM_RNA_new$name2, ]

  trajCombined_new <- trajRNA2_sub
  assay(trajCombined_new) <- t(apply(assay(trajRNA2_sub), 1, scale)) + 
    t(apply(assay(trajMM2_sub), 1, scale))

  combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new, returnMat = TRUE,
    varCutOff = 0)

  rowOrder <- match(rownames(combinedMat_new), rownames(trajRNA2_sub))
  ht1 <- plotTrajectoryHeatmap_kristen(trajRNA2_sub,
    pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0,
    rowOrder = rowOrder, max_ident= max_ident_RNA, colorMap = colorMap)
  ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap)

  pdf(paste0(save_path, "images/pseudotime/",
    pseudotime_scores, "/ATAC_correlate_GIM_MM_with_RNA_smooth100cells_var",
    var_cutoff, ".pdf"), height = 11, width = 16)
  print(ht1 + ht2)
  dev.off()
} 

var_cutoff <- 0.25
corGIM_MM_new <- correlateTrajectories(trajGIM_new, trajMM_new,
  varCutOff1 = var_cutoff, varCutOff2 = var_cutoff, corCutOff = 0.25)
print("correlation")
print(nrow(corGIM_MM_new[[1]]) > 1)
if(nrow(corGIM_MM_new[[1]]) > 1) {

  cor_df <- corGIM_MM_new[[1]][grep("z:", corGIM_MM_new[[1]]$name2), ]
  trajGIM2_sub <- trajGIM_new[cor_df$name1, ]
  trajMM2_sub <- trajMM_new[cor_df$name2, ]

  trajCombined_new <- trajGIM2_sub
  assay(trajCombined_new) <- t(apply(assay(trajGIM2_sub), 1, scale)) + 
    t(apply(assay(trajMM2_sub), 1, scale))

  combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new,
    returnMat = TRUE, varCutOff = 0)

  rowOrder <- match(rownames(combinedMat_new), rownames(trajGIM2_sub))
  ht1 <- plotTrajectoryHeatmap_kristen(trajGIM2_sub, 
    pal = paletteContinuous(set = "horizonExtra"), 
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap,
    labelRows = TRUE, labelTop = 100)
  ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap,
    labelRows = TRUE, labelTop = 100)

  pdf(paste0(save_path, "images/pseudotime/",
    pseudotime_scores, "/ATAC_correlate_GIM_MM_smooth100cells_var",
    var_cutoff, ".pdf"), height = 11, width = 16)
  print(ht1 + ht2)
  dev.off()



  # repeat with RNA heatmap GIM genes
  cor_df_RNA <- cor_df[cor_df$matchname1 %in% rownames(trajRNA_new), ]

  trajRNA2_sub <- trajRNA_new[cor_df_RNA$matchname1, ]
  trajMM2_sub <- trajMM_new[cor_df_RNA$name2, ]

  trajCombined_new <- trajRNA2_sub
  assay(trajCombined_new) <- t(apply(assay(trajRNA2_sub), 1, scale)) + 
    t(apply(assay(trajMM2_sub), 1, scale))

  combinedMat_new <- plotTrajectoryHeatmap(trajCombined_new, returnMat = TRUE,
    varCutOff = 0)

  rowOrder <- match(rownames(combinedMat_new), rownames(trajRNA2_sub))
  ht1 <- plotTrajectoryHeatmap_kristen(trajRNA2_sub,
    pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0,
    rowOrder = rowOrder, max_ident= max_ident_RNA, colorMap = colorMap,
    labelRows = TRUE, labelTop = 100)
  ht2 <- plotTrajectoryHeatmap_kristen(trajMM2_sub,
    pal = paletteContinuous(set = "solarExtra"),
    varCutOff = 0, rowOrder = rowOrder,
    max_ident= max_ident_ATAC, colorMap = colorMap,
    labelRows = TRUE, labelTop = 100)

  pdf(paste0(save_path, "images/pseudotime/",
    pseudotime_scores, "/ATAC_correlate_GIM_MM_with_RNA_smooth100cells_var",
    var_cutoff, ".pdf"), height = 11, width = 16)
  print(ht1 + ht2)
  dev.off()
} 


saveRDS(groupListATAC, paste0(save_path, "files/slingshot_files/slingshot_ATAC_groups.rda"))
saveRDS(groupListRNA, paste0(save_path, "files/slingshot_files/slingshot_RNA_groups.rda"))
saveRDS(max_ident_RNA, paste0(save_path, "files/slingshot_files/slingshot_cell_groups_RNA.rda"))
saveRDS(max_ident_ATAC, paste0(save_path, "files/slingshot_files/slingshot_cell_groups_ATAC.rda"))
saveRDS(colorMap, paste0(save_path, "files/slingshot_files/slingshot_colorMap.rda"))
saveRDS(trajMM_new, paste0(save_path, "files/slingshot_files/slingshot_Motif_matrix.rda"))
saveRDS(trajGIM_new, paste0(save_path, "files/slingshot_files/slingshot_Gene_integration_matrix.rda"))
saveRDS(trajGSM_new, paste0(save_path, "files/slingshot_files/slingshot_Gene_score_matrix.rda"))
saveRDS(trajRNA_new, paste0(save_path, "files/slingshot_files/slingshot_RNA_matrix.rda"))
saveRDS(trajPM_new, paste0(save_path, "files/slingshot_files/slingshot_Peak_matrix.rda"))


print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("pseudotime of correlated genes")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

corr_genes <- c(corGIM_MM_new[[1]]$matchname1, corGSM_MM_new[[1]]$matchname1)

corr_genes <- corr_genes[corr_genes %in% rownames(seurat_pseudotime)]

corr_genes_pseudotime2 <- corr_genes

seurat_pseudotime <- seurat_stable
seurat_pseudotime <- Rmagic::magic(seurat_pseudotime, genes = corr_genes)
DefaultAssay(seurat_pseudotime) <- "MAGIC_SCT"

motifs <- corr_genes
markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
  useMatrix = "MotifMatrix")

markerMotifs <- grep("z:", markerMotifs, value = TRUE)




# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_correlated_genes_smooth100cells.pdf"),
#   height = 9, width = 7)
# plot_list <- lapply(1:length(corr_genes), function(x){
#   print(x)
#   print(corr_genes[x])
#   # Figure this out...
#   gene_motif <- markerMotifs[grepl(paste0("z:", corr_genes[x], "_")
#     , markerMotifs)]
#   print(gene_motif)
#   if(length(gene_motif) > 0){

#     # Plot RNA
#     pseudotime_plot_gene <- plot_sling_pseudotime(seurat_object = seurat_pseudotime,
#       y_val = corr_genes[[x]],
#       col_by = "slingshot_curve1", pseudotime_curve = "slingshot_curve1",
#       color = "viridis")
#     pseudotime_plot_gene <- pseudotime_plot_gene + ggplot2::ggtitle("RNA expression")

#     # Plot ATAC GeneScoreMatrix
#     pseudotime_plot_ATAC <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#       archr_matrix = "GeneScoreMatrix", gene = corr_genes[[x]],
#       trajectory = "RNA_pseudotime1", color = "viridis",
#       col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#       impute = TRUE)
#     pseudotime_plot_ATAC <- pseudotime_plot_ATAC + ggplot2::ggtitle("ATAC gene scores")

#     # Plot ATAC Motif Matrix
#     pseudotime_plot_ATAC_dev <- ATAC_plot_pseudotime(ArchRProj = projTonsils5,
#       archr_matrix = "MotifMatrix", gene = gene_motif,
#       trajectory = "RNA_pseudotime1", color = "viridis",
#       col_by = "RNA_pseudotime1", imputeWeights = getImputeWeights(projTonsils5),
#       impute = TRUE, useSeqnames = "z")
#     pseudotime_plot_ATAC_dev <- pseudotime_plot_ATAC_dev + ggplot2::ggtitle("ATAC deviation z scores")


#     # make grid of three plots object
#     gridExtra::grid.arrange(grobs = list(pseudotime_plot_gene,
#       pseudotime_plot_ATAC, pseudotime_plot_ATAC_dev), nrow = 3)
#   }
# })

# dev.off()

##########################################
# TF footprints


plot_footprints <- function(ArchRProj, gene_list, groupBy, gene_list_name,
  normMethod = "Subtract"){
  motifPositions <- getPositions(ArchRProj)
  markerMotifs <- unlist(lapply(gene_list, function(x) grep(paste0("^", x, "_"),
    names(motifPositions), value = TRUE)))
  # Add group coverage
  seFoot <- getFootprints(
    ArchRProj = ArchRProj, 
    positions = motifPositions[markerMotifs], 
    groupBy = groupBy
  )

  plot <- plotFootprints(
    seFoot = seFoot,
    ArchRProj = ArchRProj, 
    normMethod = normMethod,
    plotName = "Footprints-Subtract-Bias",
    addDOC = FALSE,
    smoothWindow = 5,
    plot = FALSE,
    force = TRUE
  )
  return(plot)
}

# Split up by quantiles and
# 1. find correlations across quantiles
# 2. plot footprints along quantiles
# Can think of two ways to do this. Evenly across psuedotime or evenly across cell
# numbers
# First I will do numbers, then I will do pseudotime
nquantiles <- 10

group_in_quantile <- round(length(groupListATAC)/nquantiles)

groupListATAC2 <- lapply(1:nquantiles, function(x){
  quantile_list <- lapply(1:group_in_quantile, function(y){
    index <- group_in_quantile*(x-1) + y
    print(index)
    if(index > length(groupListATAC)){
      return(NULL)
    } else {
      cells <- groupListATAC[[index]]
      cells <- data.frame(cells = cells)
      return(cells) 
    }

  })
  quantile_list_group <- do.call(rbind, quantile_list)
  start_index <- group_in_quantile*(x-1) + 1
  starting_pseudotime <- pseudotimeList[[start_index]]
  end_index <- start_index + group_in_quantile - 1
  if(end_index > length(groupListATAC)){
      ending_pseudotime <- pseudotimeList[[length(groupListATAC)]]
  } else {
    ending_pseudotime <- pseudotimeList[[end_index]]
  }
  print(starting_pseudotime)
  print(ending_pseudotime)
  quantile_list_group$time <- paste0("T.", starting_pseudotime[1], "_",
    ending_pseudotime[2])
  return(quantile_list_group)
})


cell_quantile <- do.call(rbind, groupListATAC2)


projTonsils5 <- addCellColData(projTonsils5, data = cell_quantile$time,
  name = "pseudotime_quantile", cells = as.character(cell_quantile$cells),
  force = TRUE)

projTonsils <- subsetArchRProject(ArchRProj = projTonsils5,
  cells = as.character(cell_quantile$cells),
  outputDirectory = paste0(ATAC_dir, pseudotime_scores))

# seGroupMotif <- getGroupSE(ArchRProj = projTonsils,
#   useMatrix = "MotifMatrix", groupBy = "pseudotime_quantile")

# seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

# rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
#   rowMaxs(assay(seZ) - assay(seZ)[,x])
# }) %>% Reduce("cbind", .) %>% rowMaxs


# quantiles <- unique(cell_quantile$time)

# correlations <- lapply(seq_along(quantiles), function(x){
#   group_correlation <- findCorrelations(ArchRProj = projTonsils,
#     atac_cluster = "pseudotime_quantile",
#     useMatrix = "GeneScoreMatrix",
#     cell_type = quantiles[x],
#     group_name = quantiles[x], seZ = seZ)

#   })

# geom_tile_plot_gsm <- plot_correlations_geomtile(ArchRProj = projTonsils,
#   atac_cluster = "pseudotime_quantile",
#   cell_type_levels = unique(cell_quantile$time),
#   useMatrix = "GeneScoreMatrix")

# geom_tile_plot_gim <- plot_correlations_geomtile(ArchRProj = projTonsils,
#   atac_cluster = "pseudotime_quantile",
#   cell_type_levels = unique(cell_quantile$time),
#   useMatrix = "GeneIntegrationMatrix")


# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/TF_regulators_over_pseudotime_number_cells.pdf"),
#   height = 22, width = 10)

# geom_tile_plot_gsm + ggtitle("Gene Score Matrix vs Motif Matrix") +
#   xlab("pseudotime quantile (split by number of cells)")

# geom_tile_plot_gim + ggtitle("Gene Integration Matrix vs Motif Matrix") +
#   xlab("pseudotime quantile (split by number of cells)")

# dev.off()

projTonsils <- addGroupCoverages(ArchRProj = projTonsils,
  groupBy = "pseudotime_quantile")
plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators,
  groupBy = "pseudotime_quantile")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/tf_regulator_footprint_number_cells.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()

plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators_bcell,
  groupBy = "pseudotime_quantile")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/bCell_tf_regulator_footprint_number_cells.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()


plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators,
  groupBy = "Clusters2")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/tf_regulator_footprint_cell_type.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()

plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators_bcell,
  groupBy = "Clusters2")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/bCell_tf_regulator_footprint_cell_type.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()

# First do across all cell types
# Transcription regulators
txn_regulators
marker_txn_regulators <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_tf_regulator_expression.pdf"),
#   height = 9, width = 7)

# Transcription regulators B cell

# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_b_cell_tf_regulator_expression.pdf"),
#   height = 9, width = 7)
txn_regulators_bcell

# correlated genes pseudotime1
# pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
#   "/pseudotime_correlated_genes.pdf"),
#   height = 9, width = 7)
corr_genes_pseudotime1

# corr_genes_pseudtoime 2
pdf(paste0(save_path, "images/pseudotime/", pseudotime_scores,
  "/pseudotime_heatmaps_cell_names_smooth100cells.pdf"),
  height = 9, width = 12)
corr_genes_pseudotime2

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
seFoot <- getFootprints(
  ArchRProj = projTonsils5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)

# Repeat with 20 quantiles
nquantiles <- 20

group_in_quantile <- round(length(groupListATAC)/nquantiles)

groupListATAC2 <- lapply(1:nquantiles, function(x){
  quantile_list <- lapply(1:group_in_quantile, function(y){
    index <- group_in_quantile*(x-1) + y
    print(index)
    if(index > length(groupListATAC)){
      return(NULL)
    } else {
      cells <- groupListATAC[[index]]
      cells <- data.frame(cells = cells)
      return(cells) 
    }

  })
  quantile_list_group <- do.call(rbind, quantile_list)
  start_index <- group_in_quantile*(x-1) + 1
  starting_pseudotime <- pseudotimeList[[start_index]]
  end_index <- start_index + group_in_quantile - 1
  if(end_index > length(groupListATAC)){
      ending_pseudotime <- pseudotimeList[[length(groupListATAC)]]
  } else {
    ending_pseudotime <- pseudotimeList[[end_index]]
  }
  print(starting_pseudotime)
  print(ending_pseudotime)
  quantile_list_group$time <- paste0("T.", starting_pseudotime[1], "_",
    ending_pseudotime[2])
  return(quantile_list_group)
})


cell_quantile <- do.call(rbind, groupListATAC2)


projTonsils <- addCellColData(projTonsils, data = cell_quantile$time,
  name = "pseudotime_quantile_20", cells = as.character(cell_quantile$cells),
  force = TRUE)

projTonsils <- addGroupCoverages(ArchRProj = projTonsils,
  groupBy = "pseudotime_quantile_20")
plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators,
  groupBy = "pseudotime_quantile_20")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/tf_regulator_footprint_number_cells_20.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()

plot_list <- plot_footprints(ArchRProj = projTonsils, gene_list = txn_regulators_bcell,
  groupBy = "pseudotime_quantile_20")

pdf(paste0(save_path ,"images/pseudotime/", pseudotime_scores,
  "/bCell_tf_regulator_footprint_number_cells_20.pdf"), height = 6, width = 4)
plot_list2 <- lapply(seq_along(plot_list), function(x){
  if(x!=1){
    grid::grid.newpage()
  }
  grid::grid.draw(plot_list[[x]])
  })

dev.off()


# Save seurat object
