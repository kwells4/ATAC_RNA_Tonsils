create_seurat_object <- function(sample, data_dir, ADT = TRUE, min_features = 200,
  min_cells = 3){
  sample_data <- Read10X(data.dir = paste0(data_dir, sample,
    "_count/outs/filtered_feature_bc_matrix/"))
  if (ADT){
    sample_object <- CreateSeuratObject(counts = sample_data[["Gene Expression"]],
      project = sample, min.cells = min_cells, min.features = min_features)
    sample_object[["ADT"]] <- CreateAssayObject(counts =
      sample_data[["Antibody Capture"]][ , Cells(sample_object)])
  } else {
    sample_object <- CreateSeuratObject(counts = sample_data,
      project = sample, min.cells = min_cells,
      min.features = min_features)
  }
  return(sample_object)
}

create_seurat_object_snakemake <- function(sample, count_path, ADT = TRUE,
  min_features = 200, min_cells = 3){
  print(count_path)
  print(sample)
  sample_data <- Read10X(data.dir = paste0(count_path, "/filtered_feature_bc_matrix/"))
  if (ADT){
    sample_object <- CreateSeuratObject(counts = sample_data[["Gene Expression"]],
      project = sample, min.cells = min_cells, min.features = min_features)
    sample_object[["ADT"]] <- CreateAssayObject(counts =
      sample_data[["Antibody Capture"]][ , Cells(sample_object)])
  } else {
    sample_object <- CreateSeuratObject(counts = sample_data,
      project = sample, min.cells = min_cells,
      min.features = min_features)
  }
  return(sample_object)
}

initial_processing <- function(sample_object, sample, data_dir, save_dir,
                               project_name = "scRNA_seq",
                               ADT = TRUE, min_cells = 3,
                               min_features = 200, selection_method = "vst",
                               nfeatures = 2000, percent_mito = 10,
                               ngene = c(200, 7500), nADT = 4000,
                               sctransform = FALSE, vars_to_regress = NULL){
  print(sample)
  pdf(paste0(save_dir, "images/qualityPlots_", sample, ".pdf"))
  sample_data <- Read10X(data.dir = paste0(data_dir, sample,
    "_count/outs/filtered_feature_bc_matrix/"))
  if (ADT){
    sample_object <- CreateSeuratObject(counts = sample_data[["Gene Expression"]],
      project = project_name, min.cells = min_cells, min.features = min_features)
    sample_object[["ADT"]] <- CreateAssayObject(
      counts = sample_data[["Antibody Capture"]][ , Cells(sample_object)])
  } else {
    sample_object <- CreateSeuratObject(counts = sample_data,
      project = project_name, min.cells = min_cells,
      min.features = min_features)
  }
  
  sample_object[["percent.mt"]] <- PercentageFeatureSet(sample_object,
    pattern = "^MT-")
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
  plot(FeatureScatter(sample_object, feature1 = "percent.mt", "ACTB"))
  if (ADT){
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2] & nCount_ADT < nADT)
  } else {
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2])
  }
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
 if(sctransform){
    sample_object <- SCTransform(sample_object, vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    sample_object <- NormalizeData(sample_object, verbose = FALSE)
    sample_object <- FindVariableFeatures(sample_object, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    sample_object <- ScaleData(sample_object, vars.to.regress = vars_to_regress)
  }
  if (ADT) {
    sample_object <- NormalizeData(sample_object, assay = "ADT",
      normalization.method = "CLR")
    plot(VlnPlot(sample_object, features = c("nCount_ADT", "nFeature_ADT")))
    #DefaultAssay(object = sample_object) <- "ADT"

  }

  
  dev.off()
  return(sample_object)

}

initial_processing_snakemake <- function(sample_object, sample, save_dir,
                               project_name = "scRNA_seq",
                               ADT = TRUE, selection_method = "vst",
                               nfeatures = 2000, percent_mito = 10,
                               ngene = c(200, 7500), nADT = 4000,
                               sctransform = FALSE, vars_to_regress = NULL){
  print(sample)
  pdf(paste0(save_dir, "images/qualityPlots_", sample, ".pdf"))
  sample_object[["percent.mt"]] <- PercentageFeatureSet(sample_object,
  	pattern = "^MT-")
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  	ncol = 3))
  plot(FeatureScatter(sample_object, feature1 = "percent.mt", "ACTB"))
  if (ADT){
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2] & nCount_ADT < nADT)
  } else {
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2])
  }
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  	ncol = 3))
 if(sctransform){
    sample_object <- SCTransform(sample_object, vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    sample_object <- NormalizeData(sample_object, verbose = FALSE)
    sample_object <- FindVariableFeatures(sample_object, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    sample_object <- ScaleData(sample_object, vars.to.regress = vars_to_regress)
  }
  if (ADT) {
    sample_object <- NormalizeData(sample_object, assay = "ADT",
    	normalization.method = "CLR")
  	plot(VlnPlot(sample_object, features = c("nCount_ADT", "nFeature_ADT")))
  	#DefaultAssay(object = sample_object) <- "ADT"

  }

  
  dev.off()
  return(sample_object)

}

initial_processing_merge <- function(object1, object2, save_dir,
                                     cell_ids, project_name = "scRNA_seq",
                                     ADT = TRUE, min_cells = 3,
                                     selection_method = "vst",
                                     nfeatures = 2000, percent_mito = 10,
                                     ngene = c(200, 7500), nADT = 4000,
                                     sctransform = FALSE, vars_to_regress = NULL){
  pdf(paste0(save_dir, "images/qualityPlots_", sample, ".pdf"))
  sample_object <- merge(x = object1, y = object2,
    add.cell.ids = cell_ids, project = project_name)
  sample_object[["percent.mt"]] <- PercentageFeatureSet(sample_object,
    pattern = "^MT-")
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
  plot(FeatureScatter(sample_object, feature1 = "percent.mt", "ACTB"))
  if (ADT){
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2]& nCount_ADT < nADT)
  } else {
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2])
  }
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
  if(sctransform){
    sample_object <- SCTransform(sample_object, vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    sample_object <- NormalizeData(sample_object, verbose = FALSE)
    sample_object <- FindVariableFeatures(sample_object, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    sample_object <- ScaleData(sample_object, vars.to.regress = vars_to_regress)
  }
  if (ADT) {
    sample_object <- NormalizeData(sample_object, assay = "ADT",
      normalization.method = "CLR")
    plot(VlnPlot(sample_object, features = c("nCount_ADT", "nFeature_ADT")))
    #DefaultAssay(object = sample_object) <- "ADT"

  }

  
  dev.off()
  return(sample_object)

}

initial_processing_merge_two <- function(sample1, sample2, seurat_list,
                                     merge_name, project_name = "scRNA_seq",
                                     ADT = TRUE, selection_method = "vst",
                                     nfeatures = 2000, percent_mito = 10,
                                     ngene = c(200, 7500), nADT = 4000,
                                     sctransform = FALSE, vars_to_regress = NULL){
  pdf(paste0(save_dir, "images/qualityPlots_", merge_name, ".pdf"))
  sample_object_1 <- seurat_list[[sample1]]
  sample_object_2 <- seurat_list[[sample2]]
  sample_object <- merge(sample_object_1, sample_object_2,
    add.cell.ids = c(sample1, sample2), project = project_name)
  sample_object[["percent.mt"]] <- PercentageFeatureSet(sample_object,
    pattern = "^MT-")
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
  plot(FeatureScatter(sample_object, feature1 = "percent.mt", "ACTB"))
  if (ADT){
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2]& nCount_ADT < nADT)
  } else {
    sample_object <- subset(x = sample_object, subset = percent.mt < percent_mito &
      nFeature_RNA > ngene[1] & nFeature_RNA < ngene[2])
  }
  plot(VlnPlot(sample_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3))
  if(sctransform){
    sample_object <- SCTransform(sample_object, vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    sample_object <- NormalizeData(sample_object, verbose = FALSE)
    sample_object <- FindVariableFeatures(sample_object, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    sample_object <- ScaleData(sample_object, vars.to.regress = vars_to_regress)
  }
  if (ADT) {
    sample_object <- NormalizeData(sample_object, assay = "ADT",
      normalization.method = "CLR")
    plot(VlnPlot(sample_object, features = c("nCount_ADT", "nFeature_ADT")))
    #DefaultAssay(object = sample_object) <- "ADT"

  }

  
  dev.off()
  return(sample_object)

}

initial_processing_from_matrix <- function(sample, data_dir, rna_file, ADT = TRUE,
  min_cells = 3, min_features = 200, selection_method = "vst", nfeatures = 2000,
  percent_mito = 10, ngene = c(200, 7500), nADT = 4000, sctransform = FALSE,
  vars_to_regress = NULL){
  rna_matrix <- readRDS(paste0(data_dir, rna_file))
  sample_object <- CreateSeuratObject(counts = rna_matrix, project = sample,
    min.cells = min_cells, min.features = min_features)
  if(ADT){
    adt_file <- gsub("RNA", "ADT", rna_file)
    adt_matrix <- readRDS(paste0(data_dir, adt_file))
    # Dumb workaround to get rid of "NA" rows
    adt_matrix[is.na(adt_matrix)] <- 0
    adt_matrix <- adt_matrix[rowSums(adt_matrix) > 0, ]
    if(all_equal(colnames(adt_matrix), colnames(rna_matrix))){
      sample_object[["ADT"]] <- CreateAssayObject(counts = adt_matrix)
    } else {
      sample_object[["ADT"]] <- CreateAssayObject(
        counts = adt_matrix[ , Cells(sample_object)])
    }
  }
  if (ADT){
    # sample_object <- subset(x = sample_object, nFeature_RNA > ngene[1] & 
    #   nFeature_RNA < ngene[2] & nCount_ADT < nADT)
    sample_object <- subset(x = sample_object, subset = nFeature_RNA > ngene[1] &
      nFeature_RNA < ngene[2])
  } else {
    sample_object <- subset(x = sample_object, subset = nFeature_RNA > ngene[1] &
      nFeature_RNA < ngene[2])
  }
  if(sctransform){
    sample_object <- SCTransform(sample_object, vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    sample_object <- NormalizeData(sample_object, verbose = FALSE)
    sample_object <- FindVariableFeatures(sample_object, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    sample_object <- ScaleData(sample_object, vars.to.regress = vars_to_regress)
  }
  if (ADT) {
    sample_object <- NormalizeData(sample_object, assay = "ADT",
      normalization.method = "CLR")
    sample_object <- ScaleData(sample_object, assay = "ADT")

  }
  return(sample_object)
}

integrate_samples <- function(object_list, sctransform = TRUE){
  if(sctransform) {
    object_features <- SelectIntegrationFeatures(object.list = object_list,
      nfeatures = 3000)
    object_list <- PrepSCTIntegration(object.list = object_list,
      anchor.feature = object_features)
    object_anchors <- FindIntegrationAnchors(object.list = object_list,
      normalization.method = "SCT", anchor.features = object_features)
    object_integrated <- IntegrateData(anchorset = object_anchors, 
      normalization.method = "SCT")
  } else {
    print("Not written to handle non sctransform data yet")
  }
  return(object_integrated)
}


PCA_dimRed <- function(sample_object, sample_name, save_dir){
  print(sample_name)
  pdf(paste0(save_dir, "images/PCA_jackstraw_", sample_name, ".pdf"))
  sample_object <- RunPCA(sample_object,
    features = VariableFeatures(object = sample_object))
  plot(VizDimLoadings(sample_object, dims = 1:2, reduction = "pca"))
  plot(DimPlot(sample_object, reduction = "pca"))
  plot(DimPlot(sample_object, group.by = "orig.ident", reduction = "pca"))
  plot(FeaturePlot(sample_object, features = "percent.mt", reduction = "pca") +
    ggplot2::scale_color_gradient(low = "#00AFBB", high = "#FC4E07"))
  plot(FeaturePlot(sample_object, features = "nFeature_RNA", reduction = "pca") +
    ggplot2::scale_color_gradient(low = "#00AFBB", high = "#FC4E07"))
  plot(FeaturePlot(sample_object, features = "nCount_RNA", reduction = "pca") +
    ggplot2::scale_color_gradient(low = "#00AFBB", high = "#FC4E07"))
  # Currently crashing, check tomorrow
  #plot(DimHeatmap(sample_object, dims = 1:12, cells = 500, balanced = TRUE))
  # When submitting this, uncomment
  sample_object <- JackStraw(sample_object, num.replicate = 100)
  sample_object <- ScoreJackStraw(sample_object, dims = 1:20)
  plot(JackStrawPlot(sample_object, dims = 1:20))
  plot(ElbowPlot(sample_object))
  dev.off()
  return(sample_object)
}

group_cells <- function(sample_object, sample_name, save_dir, nPCs = 10,
  resolution = 0.8, plot_ADT = FALSE, ...){
  print(sample_name)
  save_plot <- paste0(save_dir, "images/UMAP_", sample_name, ".pdf")
  sample_object <- FindNeighbors(sample_object, dims = 1:nPCs)
  sample_object <- FindClusters(sample_object, resolution = resolution)
  sample_object <- RunUMAP(sample_object, umap.method = "umap-learn",
    metric = "correlation", dims = 1:nPCs)
  col_by_list <- c("cluster", "orig.ident")
  if (plot_ADT){
  	ADT_genes <- rownames(GetAssayData(object = sample_object, slot = "data",
      assay = "ADT"))
  	#ADT_genes <- paste0("adt_", ADT_genes)
    col_by_list <- c(col_by_list, ADT_genes)
  }
  plotDimRed(sample_object = sample_object, save_plot = save_plot,
    col_by = col_by_list, return_plot = FALSE, ...)

  return(sample_object)
}

# group_cells <- function(sample_object, sample_name, save_dir, nPCs = 10,
#   resolution = 0.8, plot_ADT = FALSE){
#   print(sample_name)
#   pdf(paste0(save_dir, "images/UMAP_", sample_name, ".pdf"))
#   sample_object <- FindNeighbors(sample_object, dims = 1:nPCs)
#   sample_object <- FindClusters(sample_object, resolution = resolution)
#   sample_object <- RunUMAP(sample_object, umap.method = "umap-learn",
#   metric = "correlation", dims = 1:nPCs)
#   print(DimPlot(sample_object, reduction = "umap"))
#   plot(DimPlot(sample_object, group.by = "orig.ident", reduction = "umap"))
#   if (plot_ADT){
#     DefaultAssay(object = sample_object) <- "ADT"
#     ADT_genes <- rownames(GetAssayData(object = sample_object, slot = "data"))
#     DefaultAssay(object = sample_object) <- "RNA"
#     ADT_genes <- paste0("adt_", ADT_genes)
#     lapply(ADT_genes, function(x) {
#       print(FeaturePlot(sample_object, features = x, reduction = "umap") +
#         ggplot2::scale_color_gradient(low = "#00AFBB", high = "#FC4E07"))
#       })
#   }
#   dev.off()
#   return(sample_object)
# }

plot_umap <- function(sample_object, sample_name, save_dir, feature_list) {
  pdf(paste0(save_dir, "images/umap_", sample_name, "_qc.pdf"))
  lapply(feature_list, function(x) {
    plot(FeaturePlot(sample_object, features = x, reduction = "umap") +
    ggplot2::scale_color_gradient(low = "#00AFBB", high = "#FC4E07"))
  })
  dev.off()
}

cell_markers <- function(sample_object, sample_name, save_dir, sctransform = TRUE,
  seed = 0, cell_k = NULL, gene_k = NULL, cluster_color = NULL){
  print(sample_name)
  if(sctransform){
    slot <- "scale.data"
  } else {
    slot <- "data"
  }
  pdf(paste0(save_dir, "images/markers_", sample_name, ".pdf"),
  	width = 10, height = 20)
  # idents <- levels(Idents(sample_object))
  # idents <- setNames(idents, idents)
  # sample_markers <- lapply(idents, function(x) FindMarkers(object = sample_object,
  #   ident.1 = x, slot = slot, only.pos = TRUE, reduction = reduction,
  #   min.pct = 0.25, logfc.threshold = 0.25))
  # sample_markers <- unlist(sample_markers)
  sample_markers <- FindAllMarkers(sample_object, only.pos = TRUE,
  	min.pct = 0.25, logfc.threshold = 0.25, slot = slot)
  print(head(sample_markers))
  if(sctransform){
    top10 <- sample_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
  } else {
    top10 <- sample_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  }

  sample_object <- GetResidual(sample_object, features = top10$gene, assay = "SCT")
  
  plot(DoHeatmap(sample_object, features = top10$gene, label = FALSE) + NoLegend())
  unique_genes <- unique(top10$gene)
  genes <- unique_genes[unique_genes %in% rownames(GetAssayData(sample_object,
  	slot = "scale.data"))]
  plot_data <- t(FetchData(object = sample_object, vars = genes, slot = "scale.data"))
  annotations <- FetchData(object = sample_object, vars = "ident")

  if (is.null(cluster_color)){
    nColors <- length(levels(sample_object))
    cluster_color <- grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(9, "Set1"))(nColors)
    names(cluster_color) <- levels(Idents(sample_object))
  }
  
  if (is.null(cell_k)){
    cell_k <- length(levels(sample_object))
  }

  if (is.null(gene_k)){
    gene_k <- 10
  }
  # cluster_rows <- tglkmeans::TGL_kmeans(df = plot_data, k = gene_k, seed = seed,
  #   id_column = FALSE)
  # cluster_rows_df <- data.frame(gene_cluster = cluster_rows$cluster)
  # rownames(cluster_rows_df) <- rownames(plot_data)
  # cluster_rows_df <- cluster_rows_df[order(cluster_rows_df$gene_cluster), , drop = FALSE]

  # cluster_columns <- tglkmeans::TGL_kmeans(df = t(plot_data), k = cell_k,
  #   seed = seed, id_column = FALSE)
  # cluster_columns_df <- data.frame(cell_cluster = cluster_columns$cluster)
  # rownames(cluster_columns_df) <- colnames(plot_data)
  # cluster_columns_df <- cluster_columns_df[order(cluster_columns_df$cell_cluster), ,
  # drop = FALSE] 

  annotations_s <- annotations[order(annotations$ident), , drop = FALSE]
  plot_data_s <- plot_data[ , order(match(colnames(plot_data), rownames(annotations_s)))]
  # plot_data_s <- plot_data[order(match(rownames(plot_data),
  #                                      rownames(cluster_rows_df))),
  #   order(match(colnames(plot_data), rownames(annotations_s)))]

  # plot_data_s2 <- plot_data[order(match(rownames(plot_data),
  #                                       rownames(cluster_rows_df))),
  #   order(match(colnames(plot_data), rownames(cluster_columns_df)))]

  # cluster_columns_df$cell_cluster <- factor(cluster_columns_df$cell_cluster)
  # cluster_rows_df$gene_cluster <- factor(cluster_rows_df$gene_cluster)


  all_colors <- list(ident = cluster_color)
  
  # This comes from the seurat DoHeatmap function
  plot_data_s <- MinMax(plot_data_s, -2.5, 2.5)
  pheatmap(mat = plot_data_s,
  	            show_colnames = FALSE,
  	            annotation_col = annotations_s,
                #annotation_row = cluster_rows_df,
                annotation_colors = all_colors,
  	            cluster_cols = FALSE,
                cluster_rows = TRUE,
  	            color = viridis(10),
                main = "gene clustering")
  plot_data_s <- MinMax(plot_data_s, -2.5, 2.5)
  pheatmap(mat = plot_data_s,
  	            show_colnames = FALSE,
                cluster_cols = TRUE,
                cluster_rows = TRUE,
  	            #annotation_col = cluster_columns_df,
                #annotation_row = cluster_rows_df,
  	            color = viridis(10),
                main = "gene and cell clustering")
  dev.off()
  return(sample_markers)
}

significant_markers <- function(sample_object,
                                cluster_list = NULL, ...){
  print(sample_object)
  # Do DE of all clusters if no list is given
  if (is.null(cluster_list)){
    clusters <- Idents(sample_object)
    cluster_vals <- levels(factor(clusters))
  } else {
    cluster_vals <- cluster_list
  }
  
  # Determine all pair-wise combinations of clusters given
  combinations <- utils::combn(cluster_vals, 2)

  # Run DE between all pairwise clusters and return a list of all DE genes
  de_list_all <- sapply(1:ncol(combinations), function(x)
    DE_two_clusters(sample_object, combinations[1, x], combinations[2, x], ...))
  
  # Add the DE lists to the seurat object (append to an existing list or create a
    # new slot)
  if (is.null(sample_object@misc$DE)){
    sample_object@misc$DE <- de_list_all
  } else {
    previously_covered <- names(de_list_all) %in% names(sample_object@misc$DE)
    previously_covered_DE <- names(de_list_all)[previously_covered]
    if(length(previously_covered_DE) > 0) {
      print(paste0("replacing ", previously_covered_DE, " in DE slot"))
      sample_object@misc$DE <- modifyList(sample_object@misc$DE, de_list_all)
    } else {
      sample_object@misc$DE <- c(sample_object@misc$DE, de_list_all)
    }
  }

  return(sample_object)
}

DE_two_clusters <- function(seurat_object,
                            ident_1, ident_2,
                            adj_p_val        = 0.05,
                            sig_logFC        = 1.0,
                            num_genes        = 5,
                            test_use         = "wilcox",
                            logfc_threshold  = 0.01,
                            reduction        = NULL,
                            sctransform      = TRUE){
  
  # if(sctransform){
  #   slot <- "scale.data"
  # } else {
  #   slot <- "data"
  # }
  slot <- "data"
  # Find markers between two clusters
  cluster_markers <- find_markers(seurat_object, ident_1 = ident_1,
    ident_2 = ident_2, adj_p_val = adj_p_val, sig_logFC = sig_logFC,
    test_use = test_use, logfc_threshold = logfc_threshold,
    reduction = reduction, slot = slot)
    
  # Name the comparison
  slot_name <- paste0(as.character(ident_1), "vs", as.character(ident_2))
    
      
  # Make a dataframe consisting of significant and high fold change genes
  DE_sig <- cluster_markers[cluster_markers$group ==
                                "Significant&FoldChange", ]
 
  DE_sig$genes <- rownames(DE_sig)
    
  # Add newest DE to list
  de_list <- list(DE_sig)
  names(de_list) <- slot_name
  return(de_list)
}

find_markers <- function(seurat_object, ident_1, ident_2, logfc_threshold,
                         adj_p_val = 0.05, sig_logFC = 0.25, test_use = "wilcox",
                         reduction = NULL, slot = "data"){
  print(logfc_threshold)
  print(ident_1)
  print(ident_2)
  # Run FindMarkers from seurat on the two clusters that are called. Use
  # A random seed for consistent results.
  cluster_markers <- Seurat::FindMarkers(object          = seurat_object,
                                         ident.1         = ident_1,
                                         ident.2         = ident_2,
                                         random.seed     = 0,
                                         logfc.threshold = logfc_threshold,
                                         test.use        = test_use,
                                         slot            = slot,
                                         reduction       = reduction)
   
  # Add a "group" column to dataframe. Make default "NotSignificant"
  cluster_markers["group"] <- "NotSignificant"
  
  # Find rows that have a low enough adjusted p_value to be significant
  # But which has a low log fold change. Change these values to "significant"
  # in the group column
  cluster_markers[which(cluster_markers["p_val_adj"] < adj_p_val &
                        abs(cluster_markers["avg_logFC"]) <=
                          sig_logFC), "group"] <- "Significant"
  
  # Find rows that have a high p_value but large log fold change. Call these
  # "FoldChange"
  cluster_markers[which(cluster_markers["p_val_adj"] >= adj_p_val &
                        abs(cluster_markers["avg_logFC"]) >= 
                          sig_logFC), "group"] <- "FoldChange"
  
  # Find rows with small p_value and high fold change. Call these
  # "Significant&FoldChange"
  cluster_markers[which(cluster_markers["p_val_adj"] < adj_p_val &
                        abs(cluster_markers["avg_logFC"]) > 
                          sig_logFC), "group"] <- "Significant&FoldChange"

  cluster_markers$cluster_name <- "cluster_name"
  cluster_markers[which(cluster_markers["avg_logFC"] > 0),
    "cluster_name"] <- ident_1
  
  cluster_markers[which(cluster_markers["avg_logFC"] < 0),
    "cluster_name"] <- ident_2

  return(cluster_markers)
}

plotDimRed <- function(sample_object, col_by, save_plot = NULL,
                       plot_type = "umap",
                       dims_use = NULL, highlight_group = FALSE,
                       group = NULL, meta_data_col = "orig.ident",
                       return_plot = TRUE, ...) {
  print(col_by)
  print(save_plot)
  plot_list <- lapply(col_by, function(x) {
    plotDimRedSingle(seurat_object = sample_object, col_by = x, plot_type = plot_type,
               dims_use = dims_use, highlight_group = highlight_group,
               group = group, meta_data_col = meta_data_col, ...)
  })
  if (!is.null(save_plot)){
    pdf(save_plot)
    print(plot_list)
    dev.off()
  }
  return(plot_list)
}

plotDimRedSingle <- function(seurat_object, col_by, plot_type = "umap",
                            dims_use = NULL, highlight_group = FALSE,
                            group = NULL, meta_data_col = "orig.ident", ...) {
  # Determine where in Seurat object to find variable to color by
  if (col_by == "cluster" | col_by == "Cluster"){
    col_by_data <- as.data.frame(Idents(object = seurat_object))
  }else if (col_by %in% rownames(seurat_object) |
    col_by %in% colnames(seurat_object[[]])){
    col_by_data <- FetchData(object = seurat_object, vars = col_by)
  }else if (col_by %in% rownames(seurat_object[["ADT"]])){
    col_by_data <- FetchData(object = seurat_object, vars = paste0("adt_", col_by))
  }else {
    stop("col_by must be a gene, metric from meta data or 'cluster'")
  }

  # Make the name in the data frame the same regardless of what it was originally
  names(col_by_data) <- "colour_metric"

  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  # Make a data frame based on the cell embeddings from the plot type of choice
  if (plot_type %in% names(seurat_object)){
    plot_coord <- Embeddings(object = seurat_object, reduction = plot_type)
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[ , colnames(plot_coord) %in% plot_cols]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")
    rownames(plot_df) <- plot_df$Row.names
    plot_df$Row.names <- NULL
  } else {
    stop("plot type must be a dimensional reduction in Seurat object")
  }
  # Add in group information if highlighting one group.
  if (highlight_group){
    if (is.null(group)){
      stop("if highlight_group is true, group must be a value from the meta_data
            column specified")
    }
    if (!identical(rownames(seurat_object[[]]), rownames(plot_df))) {
      print("must reorder cells")
      plot_df <- plot_df[match(rownames(seurat_object[[]]),
                                       rownames(plot_df)), , drop = FALSE]
    }
    plot_df[[meta_data_col]] <- seurat_object[[meta_data_col]]
    if (is.factor(plot_df$all)){
      plot_df$all <- factor(plot_df$all,
        levels = c("all_samples", levels(plot_df$all)))
    }
    plot_df$all[!(plot_df[[meta_data_col]] %in% group)] <- "all_samples"

    # Plot as descrete
    if (!is.numeric(plot_df$colour_metric)){
      return_plot <- groupDiscretePlots(group, plot_df, axis_names = axis_names,
                                        col_by = col_by, ...)
    # Plot as continuouts
    }else{
      return_plot <- groupContinuousPlots(group, plot_df, axis_names = axis_names,
                                          col_by = col_by, ...)
    }
  }
  # Plot as discrete
  if (!is.numeric(plot_df$colour_metric)){
    return_plot <- discretePlots(plot_df, axis_names = axis_names,
                                 col_by = col_by, ...)

  # Plot as continuous
  }else{
    return_plot <- continuousPlots(plot_df, axis_names = axis_names,
                                   col_by = col_by, ...)
  }
  return(return_plot)
}



discretePlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                          color = NULL, save_plot = NULL, show_legend = TRUE,
                          size = 0.25, ggrastr = FALSE,
                          raster_width = NULL){
  base_plot <- ggplot2::ggplot(data = plot_df,
                               ggplot2::aes_(~dim1, ~dim2))

  # Add colors based on metric chosen
  if (ggrastr){
    base_plot <- base_plot +
        ggrastr::geom_point_rast(ggplot2::aes_(colour = ~colour_metric),
                            show.legend = show_legend, size = size,
                            raster.width = raster_width) +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  } else {
    base_plot <- base_plot +
        ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                            show.legend = show_legend, size = size) +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  }
  nColors <- length(levels(factor(plot_df$colour_metric)))

  # Color based on RColorBrewer if own palette isn't chosen
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_manual(values = color, name = col_by)
  }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}



continuousPlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                            color = NULL, save_plot = NULL, show_legend = TRUE,
                            size = 0.25, ggrastr = FALSE,
                            raster_width = NULL, wesanderson = FALSE){
  if (is.null(color)) {
    low <- "#00AFBB"
    high <- "#FC4E07"
  } else {
    low <- color[1]
    high <- color[2]
  }

  base_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes_(~dim1, ~dim2))

  if (ggrastr){
    base_plot <- base_plot +
      ggrastr::geom_point_rast(ggplot2::aes_(colour = ~colour_metric),
                          show.legend = show_legend, size = size,
                          raster.width = raster_width) 
  } else {
    base_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                          show.legend = show_legend, size = size)
      
  }

  if(wesanderson){
    base_plot <- base_plot + scale_color_gradientn(colours = color) 
      # scale_x_discrete(expand = c(0, 0)) +
      # scale_y_discrete(expand = c(0, 0)) + 
      # coord_equal() 
  } else {
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low,
                                    high = high, name = col_by)
  }
  base_plot <- base_plot + ggplot2::ggtitle(col_by)
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}


groupDiscretePlots <- function(group, plot_df, col_by, axis_names = c("dim1", "dim2"),
                                color = NULL, save_plot = NULL, show_legend = TRUE,
                                size = 0.25, ggrastr = FALSE,
                                raster_width = NULL) {
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1,
                                                             ~dim2))
  
  if (ggrastr) {
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                           ggplot2::aes_(~dim1, ~dim2), 
                                           color = "#DCDCDC",
                                           size = size,
                                           show.legend = FALSE,
                                           raster.width = raster_width)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                           ggplot2::aes_(~dim1, ~dim2,
                                                         color = ~all),
                                           size = size,
                                           show.legend = show_legend,
                                           raster.width = raster_width)
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                           ggplot2::aes_(~dim1, ~dim2), 
                                           color = "#DCDCDC",
                                           size = size,
                                           show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                           ggplot2::aes_(~dim1, ~dim2,
                                                         color = ~all),
                                           size = size,
                                           show.legend = show_legend)
  }
  
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(group, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2]) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  if (is.null(color)) {
    nColors <- length(levels(factor(plot2$all)))
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
   } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
   }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}

groupContinuousPlots <- function(group, plot_df, col_by, color = NULL,
                                 limits = NULL, axis_names = c("dim1", "dim2"),
                                 save_plot = NULL, show_legend = TRUE,
                                 size = 0.25, ggrastr = FALSE,
                                 raster_width = NULL) {
  plot_name_comb <- paste(group, collapse = "_")
  if (is.null(color)) {
    low <- "#00AFBB"
    high <- "#FC4E07"
  }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]

  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1, ~dim2))
  
  if (ggrastr) {
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                           ggplot2::aes_(~dim1, ~dim2), 
                                           color = "#DCDCDC",
                                           size = size,
                                           show.legend = FALSE,
                                           raster.width = raster_width)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                           ggplot2::aes_(~dim1, ~dim2,
                                                         color = ~colour_metric),
                                           size = size,
                                           show.legend = show_legend,
                                           raster.width = raster_width)
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                           ggplot2::aes_(~dim1, ~dim2), 
                                           color = "#DCDCDC",
                                           size = size,
                                           show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                           ggplot2::aes_(~dim1, ~dim2,
                                                         color = ~colour_metric),
                                           size = size,
                                           show.legend = show_legend)
  }

  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(is.null(limits)){
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by, limits = limits)
  }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}


populations_dfs_new <- function(seurat_object, sample_name, subsample = FALSE,
                                subsample_by = "orig.ident",
                                meta_data_col = "seurat_clusters"){
  if (subsample) {
    cells_use <- rownames(seurat_object[[]])[
      seurat_object[[subsample_by]] == sample_name]
    seurat_object <- subset(seurat_object, cells = cells_use)
  }
  stage_df <- data.frame(table(seurat_object[[meta_data_col]]))
  names(stage_df) <- c("cluster", "count")
  # stage_df$percent <- stage_df$count / sum(stage_df$count) * 100
  stage_df$percent <- stage_df$count / nrow(seurat_object[[]]) * 100
  stage_df$sample <- sample_name
  return(stage_df)
}

population_plots <- function(stage_df_all, color = NULL, save_plot = NULL,
  plot = "percent", title = NULL, position = "stack", sep_by = "sample",
  fill = "cluster"){
  if(plot == "percent"){
    plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~sample,
                                                                    y = ~percent,
                                                                    fill = ~cluster))
  } else if (plot == "count") {
    plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~sample,
                                                                    y = ~count,
                                                                    fill = ~cluster))
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

full_umap <- function(seurat_object, data_set, col_by, plot_type = "umap",
                      dims_use = NULL, meta_data_col = "exp", ...) {
  # Determine where in Seurat object to find variable to color by
  if (col_by == "cluster" | col_by == "Cluster"){
    col_by_data <- as.data.frame(Idents(object = seurat_object))
  }else if (col_by %in% rownames(seurat_object) |
    col_by %in% colnames(seurat_object[[]])){
    col_by_data <- FetchData(object = seurat_object, vars = col_by)
  }else if (col_by %in% rownames(seurat_object[["ADT"]])){
    col_by_data <- FetchData(object = seurat_object, vars = paste0("adt_", col_by))
  }else {
    stop("col_by must be a gene, metric from meta data or 'cluster'")
  }


  # Make the name in the data frame the same regardless of what it was originally
  names(col_by_data) <- "colour_metric"
  
  col_by_data$all <- col_by_data$colour_metric
  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  if (!identical(rownames(seurat_object[[]]), rownames(col_by_data))) {
    print("must reorder cells")
    col_by_data <- col_by_data[match(rownames(seurat_object[[]]),
                                     rownames(col_by_data)), , drop = FALSE]
  }
  meta_data_all <- seurat_object[[]]
  col_by_data[[meta_data_col]] <- meta_data_all[[meta_data_col]]
  if (is.factor(col_by_data$all)){
    col_by_data$all <- factor(col_by_data$all,
      levels = c("all_samples", levels(col_by_data$all)))
  }

  col_by_data$all[!(col_by_data[[meta_data_col]] %in% data_set)] <- "all_samples"
  if (plot_type %in% names(seurat_object)){
    plot_coord <- Embeddings(object = seurat_object, reduction = plot_type)
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[ , colnames(plot_coord) %in% plot_cols]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")

  } else {
    stop("plot type must be a dimensional reduction in dr slot")
  }

   # Plot as discrete
  if (!is.numeric(col_by_data$colour_metric)){
    return_plot <- full_discrete_plots(data_set, plot_df, axis_names = axis_names,
      col_by = col_by, ...)
  # Plot as continuous
  }else{
    return_plot <- full_continuous_plots(data_set, plot_df, col_by = col_by, ...)
  }
  return(return_plot)

}

full_discrete_plots <- function(data_set, plot_df, col_by, axis_names = c("dim1", "dim2"),
                                color = NULL, save_plot = NULL, show_legend = TRUE,
                                size = 0.25) {
  # if (!(is.null(save_plot))){
  #   extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
  #   if (extension == "pdf"){
  #     pdf(save_plot)
  #   } else if (extension == "png") {
  #     png(save_plot)
  #   } else {
  #     print("save plot must be .png or .pdf")
  #   }
  # }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1,
                                                             ~dim2))
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                         ggplot2::aes_(~dim1, ~dim2), 
                                         color = "#DCDCDC",
                                         size = size,
                                         show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                         ggplot2::aes_(~dim1, ~dim2,
                                                       color = ~all),
                                         size = size,
                                         show.legend = show_legend) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(data_set, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2]) 
  if (is.null(color)) {
    nColors <- length(levels(factor(plot2$all)))
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
   } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
   }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot)
    # print(base_plot)
    # dev.off()
  }
  return(base_plot)
}

full_continuous_plots <- function(data_set, plot_df, col_by, color = NULL,
                                  limits = NULL, axis_names = c("dim1", "dim2"),
                                  save_plot = NULL, show_legend = TRUE,
                                  size = 0.25) {
  # if (!(is.null(save_plot))){
  #   extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
  #   if (extension == "pdf"){
  #     pdf(save_plot)
  #   } else if (extension == "png") {
  #     png(save_plot)
  #   } else {
  #     print("save plot must be .png or .pdf")
  #   }
  # }
  plot_name_comb <- paste(data_set, collapse = "_")
  if (is.null(color)) {
    low <- "#00AFBB"
    high <- "#FC4E07"
  }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]

  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1, ~dim2))
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                         ggplot2::aes_(~dim1, ~dim2), 
                                         color = "#DCDCDC",
                                         size = size,
                                         show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                         ggplot2::aes_(~dim1, ~dim2,
                                                       color = ~colour_metric),
                                         size = size,
                                         show.legend = show_legend)
  
  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(is.null(limits)){
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by, limits = limits)
  }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot)
    # print(base_plot)
    # dev.off()
  }
  return(base_plot)
  
}

get_gene_lists <- function(comparison, DE_genes, output_file = NULL){
  cluster_1 <- gsub("v.*", "", comparison)
  cluster1_genes <- rownames(DE_genes[DE_genes$avg_logFC > 0,])
  cluster_2 <- gsub(".*v", "", comparison)
  cluster2_gene <- rownames(DE_genes[DE_genes$avg_logFC < 0, ])
  gene_list <- list(cluster1_genes, cluster2_gene)
  names(gene_list) <- c(cluster_1, cluster_2)
  if (!is.null(output_file)) {
    write.table(comparison, output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table(cluster_1, output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table(cluster1_genes, output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table("\n", output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table(cluster_2, output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table(cluster2_gene, output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)
    write.table("\n", output_file, append = TRUE, col.names = FALSE,
      row.names = FALSE)

  }
  return(gene_list)
}

print_DE <- function(seurat_object, output_file = NULL){
  comparisons <- names(seurat_object@misc$DE)
  comparisons <- setNames(comparisons, comparisons)

  DE_list <- lapply(comparisons, function(x) get_gene_lists(x,
    seurat_object@misc$DE[[x]],
    output_file = output_file))

  return(DE_list)
}


cluster_genes <- function(DE_list, seurat_object, output_file = NULL){
  clusters <- levels(Idents(seurat_object))
  clusters <- setNames(clusters, clusters)
  all_markers <- lapply(clusters, function(x) one_cluster(DE_list, x,
    output_file = output_file))
  return(all_markers)
}

one_cluster <- function(DE_list, cluster, output_file = NULL){
  beg <- paste0("^", cluster, "v")
    end <- paste0("v", cluster, "$")
    cluster_comparisons <- c(grep(beg, names(DE_list)),
      grep(end, names(DE_list)))
    cluster_genes <- lapply(cluster_comparisons, function(x){
      DE_list[[x]][[cluster]]
      })
    cluster_genes_unique <- unique(unlist(cluster_genes))
    if (!is.null(output_file)) {
      write.table(cluster, output_file, append = TRUE, col.names = FALSE,
        row.names = FALSE)
      write.table(cluster_genes_unique, output_file, append = TRUE,
        col.names = FALSE, row.names = FALSE)
      write.table("\n", output_file, append = TRUE, col.names = FALSE,
          row.names = FALSE)

    }
    return(cluster_genes_unique)
}


doublet_finder <- function(seurat_object, sample, nPCs = 10, sct = FALSE,
  expted_doublets = 0.039) {
  # Taken from doublet finder tutorial
  ## pK Identification (no ground-truth) --------------------------------------
  sweep.res.list_seurat <- paramSweep_v3(seurat_object, PCs = 1:nPCs, sct = sct)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = FALSE)
  bcmvn_seurat <- find.pK(sweep.stats_seurat)
  pk <- bcmvn_seurat[bcmvn_seurat$MeanBC == max(bcmvn_seurat$MeanBC), ]$pK
  pk <- as.numeric(as.character(pk))

  ## Homotypic Doublet Proportion Estimate ----------------------------------------
  annotations <- seurat_object$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(expted_doublets*length(colnames(x = seurat_object)))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies -------------------
  # Pk value comes from the bcmnv_seurat. Want the peak "MeanBC" value
  seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:nPCs, pN = 0.25, pK = pk,
    nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  col_by <- colnames(seurat_object[[]])[grepl("DF", colnames(seurat_object[[]]))]
  plotDimRed(sample_object = seurat_object,
    save_plot = paste0(save_dir, "images/umap_", sample, "_doubletFinder.pdf"),
    col_by = col_by, return_plot = FALSE)
  return(seurat_object)
}

remove_doublets <- function(seurat_object){
  doublet_col <- colnames(seurat_object[[]])[grepl("DF", colnames(seurat_object[[]]))]
  seurat_object$doublet_col <- seurat_object[[doublet_col]]
  seurat_object <- subset(x = seurat_object, subset = doublet_col == "Singlet")
  return(seurat_object)
}



featDistPlot <- function(seurat_object, geneset, cell_cycle = FALSE,
                         plot_type = "violin",
                         color = NULL, sep_by = "cluster", save_plot = NULL,
                         nrow = NULL, ncol = NULL){
  geneset <- setNames(geneset, geneset)
  print(geneset)
  if (plot_type == "jitter") {
    # Make jitter plots colored by cell cycle stage
    if(cell_cycle){
      gene_list_cycle <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        col_by = "cycle_phase", color = c("black", "red", "purple")))
          
      # Arrange all plots into one figure
      save_plot <- gridExtra::grid.arrange(grobs = gene_list_cycle,
        nrow = length(geneset))
    } else {
      # Make a jitter plot based on expression of each gene given in the gene
      # set color by stage
      
      gene_list_stage <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        color = color))

      # Make a plot consisting of all plots made above
      save_plot <- gridExtra::grid.arrange(grobs = gene_list_stage,
        nrow = length(geneset))
    }
 
  }
  if (plot_type == "violin" || plot_type == "both") {
    if (plot_type == "both"){
      plot_jitter <- TRUE
    } else {
      plot_jitter <- FALSE
    }
    gene_list_stage <- lapply(geneset, function(x) violinPlot(
      seurat_object = seurat_object, y_val = x, x_val = sep_by,
      color = color, plot_jitter = plot_jitter))

    save_plot <- gridExtra::grid.arrange(grobs = gene_list_stage,
        nrow = length(geneset))
        
  }
  if (!(is.null(save_plot))){
    pdf(save_plot)
    plot(save_plot)
    dev.off()
  }
  return(save_plot)
}

jitterPlot <- function(seurat_object, y_val, x_val,
                        col_by = NULL, color = NULL) {
  plot_data <- plotDF(seurat_object, y_val, x_val,
                            col_by)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               color = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_point() + ggplot2::geom_jitter(shape = 16) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_color_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_color_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}  


violinPlot <- function(seurat_object, y_val, x_val,
                        col_by = NULL, color = NULL,
                        plot_jitter = FALSE) {
  plot_data <- plotDF(seurat_object, y_val, x_val,
                            col_by)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               fill = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  if (plot_jitter) {
    plot_base <- plot_base + ggplot2::geom_jitter(shape = 16)
  }
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_fill_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_fill_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}

plotDF <- function(seurat_object, y_val, x_val,
                          col_by = NULL) {
  # Add y_value to a data frame used for plotting. This value can be a gene
  # or a value from meta data like nGene
    # Determine where in Seurat object to find variable to color by
  if (y_val %in% rownames(seurat_object) |
    y_val %in% colnames(seurat_object[[]])){
    plot_data <- FetchData(object = seurat_object, vars = y_val)
  }else if (y_val %in% rownames(seurat_object[["ADT"]])){
    plot_data <- FetchData(object = seurat_object, vars = paste0("adt_", y_val))
  }else {
    stop("y_val must be a gene, metric from meta data")
  }
  # Name the column
  names(plot_data) <- "y_value"
  
  # Add a column contining the x_value. This should be something discrete
  # Like timepoint or cluster
  if (x_val %in% colnames(seurat_object[[]])) {
    # Should be able to fix this. Take out as and do x_val = then don't need names()
    x_plot_data <- as.data.frame(seurat_object[[]][, x_val, drop = FALSE])
    #x_plot_data <- data.frame("x_value" = seurat_object@meta.data[, x_val,
     #                                                             drop = FALSE])
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else if (x_val == "cluster") {
    x_plot_data <- as.data.frame(Idents(seurat_object))
    #x_plot_data <- data.frame("x_value" = seurat_object@ident)
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  
  # Name the appropriate column of the plotting data
  names(plot_data)[2] <- "x_value"
  plot_data <- plot_data[match(colnames(seurat_object),
                               rownames(plot_data)), ]
  
  # Determine how to color the plot. Default is the x_value but can be any
  # discrete value.
  if (is.null(col_by)) {
    plot_data$col_by <- plot_data$x_value
  } else if (col_by %in% colnames(seurat_object[[]])) {
    plot_data$col_by <- seurat_object[[]][ , col_by]
  } else if (col_by == "cluster") {
    plot_data$col_by <- Idents(seurat_object)
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  return(plot_data)
}

harmony_RNA <- function(seurat_object, save_dir, sample_name, seurat_list = NULL,
                        nfeatures = 2000, vars_to_regress = NULL,
                        selection_method = "vst",  resolution = 0.8,
                        sctransform = TRUE, plot_ADT = TRUE,
                        ndims = 30, ...){
  save_plot <- paste0(save_dir, "images/UMAP_", sample_name, ".pdf")
  # If the seurat object was subset, it may have already been merged.
  if(!is.null(seurat_list)){
    print("merge)")
    seurat_object_merge <- merge(x = seurat_object, y = seurat_list)
  } else {
    seurat_object_merge <- seurat_object
  }
  if(sctransform){
    print("transform")
    seurat_object_merge <- SCTransform(seurat_object_merge,
      vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    seurat_object_merge <- NormalizeData(seurat_object_merge, verbose = FALSE)
    seurat_object_merge <- FindVariableFeatures(seurat_object_merge, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    seurat_object_merge <- ScaleData(seurat_object_merge,
      vars.to.regress = vars_to_regress)
  }
  print("PCA")
  seurat_object_merge <- RunPCA(seurat_object_merge, verbose = FALSE)
  print("harmony")
  seurat_object_merge <- RunHarmony(seurat_object_merge, group.by.vars = "orig.ident",
    assay.use = "SCT")
  print("UMAP")
  seurat_object_merge <- RunUMAP(seurat_object_merge, reduction = "harmony", dims = 1:ndims,
    umap.method = "umap-learn", metric = "correlation")
  print("Neighbors")
  seurat_object_merge <- FindNeighbors(seurat_object_merge, reduction = "harmony",
    dims = 1:ndims)
  print("clusters")
  seurat_object_merge <- FindClusters(seurat_object_merge, resolution = resolution)
  col_by_list <- c("cluster", "orig.ident")
  if (plot_ADT){
    ADT_genes <- rownames(GetAssayData(object = seurat_object_merge, slot = "data",
      assay = "ADT"))
    #ADT_genes <- paste0("adt_", ADT_genes)
    col_by_list <- c(col_by_list, ADT_genes)
  }
  plotDimRed(sample_object = seurat_object_merge, save_plot = save_plot,
    col_by = col_by_list, return_plot = FALSE, ...)

  return(seurat_object_merge)
}

merge_objects <- function(seurat_list, sample_name, sctransform = TRUE,
                          nfeatures = 2000, vars_to_regress = NULL,
                          selection_method = "vst"){
  print("merging")
  print(seurat_list)
  seurat_object_merge <- merge(x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)])
  print(seurat_object_merge)
  if(sctransform){
    print("transform")
    seurat_object_merge <- SCTransform(seurat_object_merge,
      vars.to.regress = vars_to_regress,
      verbose = FALSE)
  } else {
    seurat_object_merge <- NormalizeData(seurat_object_merge, verbose = FALSE)
    seurat_object_merge <- FindVariableFeatures(seurat_object_merge, verbose = FALSE,
      selection.method = selection_method, nfeatures = nfeatures)
    seurat_object_merge <- ScaleData(seurat_object_merge,
      vars.to.regress = vars_to_regress)
  }
  return(seurat_object_merge)
}

run_harmony <- function(seurat_object_merge, save_dir, sample_name,
                        resolution = 0.8, plot_ADT = TRUE,
                        ndims = 30, seed = 42,
                        variable_features = VariableFeatures(seurat_object_merge)){
  save_plot <- paste0(save_dir, "images/UMAP_", sample_name, ".pdf")
  pdf(paste0(save_dir, "images/PCA", sample_name, ".pdf"))
  print("PCA")
  seurat_object_merge <- RunPCA(seurat_object_merge, verbose = FALSE,
    features = variable_features)
  print("new")
  plot(VizDimLoadings(seurat_object_merge, dims = 1:2, reduction = "pca"))
  plot(DimPlot(seurat_object_merge, reduction = "pca"))
  dev.off()
  print("harmony")
  seurat_object_merge <- RunHarmony(seurat_object_merge, group.by.vars = "orig.ident",
    assay.use = "SCT")
  print("UMAP")
  seurat_object_merge <- RunUMAP(seurat_object_merge, reduction = "harmony", dims = 1:ndims,
    umap.method = "umap-learn", metric = "correlation", seed.use = seed)
  print("Neighbors")
  seurat_object_merge <- FindNeighbors(seurat_object_merge, reduction = "harmony",
    dims = 1:ndims)
  print("clusters")
  seurat_object_merge <- FindClusters(seurat_object_merge, resolution = resolution)
  col_by_list <- c("cluster", "orig.ident")
  if (plot_ADT){
    ADT_genes <- rownames(GetAssayData(object = seurat_object_merge, slot = "data",
      assay = "ADT"))
    #ADT_genes <- paste0("adt_", ADT_genes)
    col_by_list <- c(col_by_list, ADT_genes)
  }
  plotDimRed(sample_object = seurat_object_merge, save_plot = save_plot,
    col_by = col_by_list, return_plot = FALSE)
  return(seurat_object_merge)
}

##################################################################################
#                           ATAC FUNCTIONS                                       #
##################################################################################

create_ATAC_object <- function(sample, atac_dir, annotation_file, activity = "seurat"){
  peaks <- Read10X_h5(paste0(atac_dir, sample, "/outs/filtered_peak_bc_matrix.h5"))
  meta <- read.table(paste0(atac_dir, sample, "/outs/singlecell.csv"), sep = ",",
    header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  seurat_object <- CreateSeuratObject(counts = peaks,
                                    assay = "ATAC",
                                    project = sample,
                                    meta.data = meta,
                                    min.cells = 1)
  fragment_path <- paste0(atac_dir, sample, "/outs/fragments.tsv.gz")
  
  fragment_file_filtered <- paste0(atac_dir, sample,
    "/outs/filtered_fragments.tsv")

  if(!file.exists(paste0(fragment_file_filtered, ".bgz"))){

    FilterFragments(
        fragment.path = fragment_path,
        cells = colnames(seurat_object),
        output.path = fragment_file_filtered
    )

  }

  seurat_object <- SetFragments(
    object = seurat_object,
    file = paste0(fragment_file_filtered, ".bgz")
    )

  seurat_object$tech <- "atac"

  if(activity == "seurat"){

    activity_matrix <- CreateGeneActivityMatrix(peak.matrix = peaks,
      annotation.file = annotation_file, seq.levels = c(1:22, "X", "Y"),
      upstream = 2000, verbose = TRUE)
  } else if (activity == "signac"){
    # Extract gene coordinates, this comes from Seurat::CreateGeneActivityMatrix
    gtf <- rtracklayer::import(con = annotation.file)
    gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = c(1:22, "X", "Y"),
      pruning.mode = 'coarse')
    if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == 
      GenomeInfoDb::seqlevelsStyle(x = peaks))) {
        GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks)
    }
    gtf.genes <- gtf[gtf$type == 'gene']
    genebodyandpromoter.coords <- Extend(x = gene.coords,
      upstream = 2000, downstream = 0)

    # create a gene by cell matrix
    activity_matrix <- FeatureMatrix(
      fragments = paste0(fragment_file_filtered, ".bgz"),
      features = genebodyandpromoter.coords,
      cells = colnames(pbmc),
      chunk = 20)

    # convert rownames from chromsomal coordinates into gene names
    gene_key <- genebodyandpromoter.coords$gene_name
    names(gene_key) <- GRangesToString(grange = genebodyandpromoter.coords)
    rownames(activity_matrix) <- gene_key[rownames(activity_matrix)]
  }

  seurat_object[['ACTIVITY']] <- CreateAssayObject(counts = activity_matrix)
  seurat_object <- NormalizeData(
    object = seurat_object,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat_object$nCount_ACTIVITY)
  )

  return(seurat_object)
}

ATAC_qc <- function(seurat_object, sample, save_dir){
  seurat_object <- NucleosomeSignal(object = seurat_object)
  seurat_object$pct_reads_in_peaks <- seurat_object$peak_region_fragments / 
    seurat_object$passed_filters * 100
  seurat_object$blacklist_ratio <- seurat_object$blacklist_region_fragments / 
    seurat_object$peak_region_fragments

  pdf(paste0(save_dir, "images/qualityPlots_", sample, ".pdf"))
  plot(VlnPlot(
    object = seurat_object,
    features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1))

  plot(VlnPlot(
    object = seurat_object,
    features = 'peak_region_fragments',
    pt.size = 0.1, log = TRUE))

  plot(FeatureScatter(seurat_object,"peak_region_fragments",'nucleosome_signal', pt.size = 0.1))
  plot(FeatureScatter(seurat_object,"peak_region_fragments",'blacklist_ratio', pt.size = 0.1))
  
  seurat_object$nucleosome_group <- ifelse(seurat_object$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
  plot(PeriodPlot(object = seurat_object, group.by = 'nucleosome_group'))
  dev.off()
  return(seurat_object)

}

ATAC_initial_processing <- function(seurat_object, peak_region_fragments_low = 1000,
  peak_region_fragments_high = 20000, pct_peaks = 15,
  blacklist_rat = 0.05, nuc_signal = 10, features = 'q0'){

  # Subset based on qc metrics
  seurat_object <- subset(x = seurat_object, 
    subset = peak_region_fragments > peak_region_fragments_low &
      peak_region_fragments < peak_region_fragments_high &
      pct_reads_in_peaks > pct_peaks &
      blacklist_ratio < blacklist_rat &
      nucleosome_signal < nuc_signal)

  # Normilization using term frequency-inverse document frequency normalization
  seurat_object <- RunTFIDF(seurat_object)

  # Feature selection, using top n% of peaks. 'q75' does top 25% of peaks
  seurat_object <- FindTopFeatures(seurat_object, min.cutoff = features)
  seurat_object <- RunSVD(
    object = seurat_object,
    assay = 'ATAC',
    reduction.key = 'LSI_',
    reduction.name = 'lsi'
  )

  return(seurat_object)
}

ATAC_group_cells <- function(seurat_object, sample_name, save_dir){
  save_plot <- paste0(save_dir, "images/UMAP_", sample_name, ".pdf")
  seurat_object <- RunUMAP(object = seurat_object, reduction = 'lsi', dims = 1:30,
    umap.method = "umap-learn", metric = "correlation")
  seurat_object <- FindNeighbors(object = seurat_object, reduction = 'lsi', dims = 1:30)
  seurat_object <- FindClusters(object = seurat_object, verbose = FALSE)
  plotDimRed(sample_object = seurat_object, col_by = c("cluster", "orig.ident"),
    save_plot = save_plot, return_plot = FALSE)
  return(seurat_object)
}

label_atac_cells <- function(rna_object, atac_object, sample, save_dir,
  filter = TRUE, rna_meta = "cell_labels", ...){
  atac_object$tech <- "ATAC"
  rna_object$tech <- "RNA"

  transfer_anchors <- FindTransferAnchors(reference = rna_object,
    query = atac_object, features = VariableFeatures(object = rna_object), 
      reference.assay = "RNA", query.assay = "MERGEDACTIVITY",
      reduction = "cca")


  data <- rna_object[[rna_meta]]
  refdata <- data[[1]]
  names(refdata) <- rownames(data)
  celltype_predictions <- TransferData(anchorset = transfer_anchors,
    refdata = refdata, 
      weight.reduction = atac_object[["lsi"]])

  # if(sum(celltype_predictions$prediction.score.max < prediction_score_max) > 0){
  #   celltype_predictions[celltype_predictions$prediction.score.max < 
  #     prediction_score_max, ]$predicted.id <- "low_confidence"
  # }
  # atac_object[[paste0("max_predicted_", rna_meta)]] <- 
  #   celltype_predictions$prediction.score.max

  atac_object <- AddMetaData(atac_object,
    metadata = celltype_predictions$prediction.score.max,
    col.name = paste0("max_predicted_", rna_meta))
  #atac_object <- AddMetaData(atac_object, metadata = celltype_predictions)

  celltype_predictions$predicted.id <- factor(celltype_predictions$predicted.id)

  atac_object <- AddMetaData(atac_object,
    metadata = celltype_predictions$predicted.id,
    col.name = paste0("predicted_", rna_meta))

  # atac_object[[paste0("predicted_", rna_meta)]] <- celltype_predictions$predicted.id

  # table(atac_object$prediction.score.max > 0.5)

  # if (filter) {
  #   atac_object <- subset(atac_object, subset = prediction.score.max > 0.5)
  # }

  save_plot <- paste0(save_dir, "images/", sample, "_RNA_", rna_meta, ".pdf")

  plotDimRed(atac_object, col_by = paste0("predicted_", rna_meta),
    save_plot = save_plot, return_plot = FALSE, ...)

  return(atac_object)
}

integrate_atac_cells <- function(rna_object, atac_object, save_dir, ndims = 30,
  sctransform = TRUE){
  plan("multiprocess", workers = 4)
  options(future.globals.maxSize = 10000 * 1024^2)
  rna_object$protocol <- "RNA"
  atac_object$protocol <- "ATAC"
  if (sctransform){
    object_list <- list(rna_object, atac_object)
    object_features <- SelectIntegrationFeatures(object.list = object_list,
      nfeatures = 3000)
    object_list <- PrepSCTIntegration(object.list = object_list,
      anchor.features = object_features, verbose = TRUE)
    normalization_method <- "SCT"
  } else {
    DefaultAssay(atac_object) <- "MERGEDACTIVITY"
    normalization_method <- "LogNormalize"
    atac_object <- NormalizeData(atac_object, verbose = TRUE)
    atac_object <- FindVariableFeatures(atac_object, selection.method = "vst", 
        nfeatures = 2000, verbose = TRUE)
    object_list <- list(rna_object, atac_object)
  }
  integration_anchors <- FindIntegrationAnchors(object.list = object_list,
    dims = 1:ndims, normalization.method = normalization_method, verbose = TRUE)
  seurat_integrated <- IntegrateData(anchorset = integration_anchors,
    dims = 1:ndims, normalization.method = normalization_method, verbose = TRUE)
  DefaultAssay(seurat_integrated) <- "integrated"
  if (!sctransform){
    seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
  }
  seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
  seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:ndims, 
    umap.method = "umap-learn", metric = "correlation", )
  save_plot <- paste0(save_dir, "images/", sample, "_RNA_and_ATAC_UMAP.pdf")
  plotDimRed(seurat_integrated, col_by = c("protocol", "cell_labels",
    "predicted_cell_labels", save_plot = save_plot, return_plot = FALSE))

  return(seurat_integrated)
}


merge_seurat_ATAC <- function(seurat_1, seurat_2, sample_1 = "s1",
  sample_2 = "s2", sep_1 = c(":", "-"), sep_2 = c(":", "-"), rename_cells = FALSE,
  merge_regions = FALSE, ...){
  
  # find peaks that intersect in both datasets
  intersecting_regions <- GetIntersectingFeatures(
    object.1 = seurat_1,
    object.2 = seurat_2,
    sep.1 = sep_1,
    sep.2 = sep_2
  )

  # # choose a subset of intersecting peaks
  # peaks_use <- sample(intersecting_regions[[1]], size = 10000, replace = FALSE)

  # count fragments per cell overlapping the set of peaks in the 10x data
  overlap_peaks_1 <- FeatureMatrix(
    fragments = GetFragments(object = seurat_1, assay = 'ATAC'),
    features = StringToGRanges(intersecting_regions[[1]], sep = sep_1),
    cells = colnames(seurat_1)
  )

  overlap_peaks_2 <- FeatureMatrix(
    fragments = GetFragments(object = seurat_2, assay = 'ATAC'),
    features = StringToGRanges(intersecting_regions[[1]], sep = sep_2),
    cells = colnames(seurat_2)
  )

  # create a new assay and add it to the 10x dataset
  seurat_1[['overlapPeaks']] <- CreateAssayObject(counts = overlap_peaks_1,
    min.cells = 1)
  seurat_1 <- RunTFIDF(object = seurat_1, assay = 'overlapPeaks')

  seurat_2[['overlapPeaks']] <- CreateAssayObject(counts = overlap_peaks_2,
    min.cells = 1)
  seurat_2 <- RunTFIDF(object = seurat_2, assay = 'overlapPeaks')

  # Currently merging only changes the names in the data slot meaning some meta
  # data cells disappear
  if (rename_cells){
    seurat_1 <- RenameCells(object = seurat_1, add.cell.id = sample_1)
    seurat_2 <- RenameCells(object = seurat_2, add.cell.id = sample_2)
  }
  if (merge_regions){
    merged_seurat <- merge_regions(seurat_1 = seurat_1, seurat_2 = seurat_2)
  } else {
    merged_seurat <- merge_normal(seurat_1 = seurat_1, seurat_2 = seurat_2)
  }
  return(merged_seurat)

}


atac_activity <- function(seurat_object, annotation_file, assay_name = "ACTIVITY"){
  activity_matrix <- CreateGeneActivityMatrix(
    peak.matrix = GetAssayData(object = seurat_object),
    annotation.file = annotation_file,
    seq.levels = c(1:22, "X", "Y"),
    upstream = 2000,
    verbose = TRUE)
  seurat_object[[assay_name]] <- CreateAssayObject(counts = activity_matrix)
  all_counts <- seurat_object[[paste0("nCount_", assay_name)]]
  scale_factor <- median(all_counts[,1])
  seurat_object <- NormalizeData(
      object = seurat_object,
      assay = assay_name,
      normalization.method = 'LogNormalize',
      scale.factor = scale_factor
  )
  return(seurat_object)
}

get_overlap_peaks <- function(seurat_object, intersecting_regions, sample,
    sep = c(":", "-")){
  overlap_peaks <- FeatureMatrix(
    fragments = GetFragments(object = seurat_object, assay = 'ATAC'),
    features = StringToGRanges(intersecting_regions[[1]], sep = sep),
    cells = colnames(seurat_object)
  )

  seurat_object[['overlapPeaks']] <- CreateAssayObject(counts = overlap_peaks,
  min.cells = 1)

  seurat_object <- RunTFIDF(object = seurat_object, assay = 'overlapPeaks')

  seurat_object <- RenameCells(object = seurat_object, add.cell.id = sample)
  return(seurat_object)

}

merge_regions <- function(seurat_1, seurat_2, save_dir = NULL, combined_name = "merged",
  assay_1 = "overlapPeaks", assay_2 = "overlapPeaks"){

  merged_seurat <- MergeWithRegions(
    object.1 = seurat_1,
    object.2 = seurat_2,
    sep.1 = c("-", "-"),
    sep.2 = c("-", "-"),
    assay.1 = assay_1,
    assay.2 = assay_2
  )

  merged_seurat <- RunTFIDF(merged_seurat)
  merged_seurat <- FindTopFeatures(merged_seurat, min.cutoff = 50)
  merged_seurat <- RunSVD(merged_seurat, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
  merged_seurat <- RunUMAP(merged_seurat, reduction = 'lsi', dims = 1:30,
    umap.method = "umap-learn", metric = "correlation")
  merged_seurat <- FindNeighbors(object = merged_seurat, reduction = 'lsi', dims = 1:30)
  merged_seurat <- FindClusters(object = merged_seurat, verbose = FALSE)
  if (!is.null(save_dir)){
    save_plot <- paste0(save_dir, "images/UMAP_", combined_name, ".pdf")
    col_by_list <- c("orig.ident", "cluster")
    plotDimRed(sample_object = merged_seurat, save_plot = save_plot,
        col_by = col_by_list, return_plot = FALSE)
  }
  return(merged_seurat)
}

merge_normal <- function(seurat_1, seurat_2, save_dir = NULL, combined_name = "merged"){
  DefaultAssay(object = seurat_1) <- "overlapPeaks"
  DefaultAssay(object = seurat_2) <- "overlapPeaks"

  merged_seurat <- merge(
    x = seurat_1,
    y = seurat_2,
    merge.data = FALSE
  )

  merged_seurat <- RunTFIDF(merged_seurat)
  merged_seurat <- FindTopFeatures(merged_seurat, min.cutoff = 50)
  merged_seurat <- RunSVD(merged_seurat, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
  merged_seurat <- RunUMAP(merged_seurat, reduction = 'lsi', dims = 1:30,
    umap.method = "umap-learn", metric = "correlation")
  merged_seurat <- FindNeighbors(object = merged_seurat, reduction = 'lsi', dims = 1:30)
  merged_seurat <- FindClusters(object = merged_seurat, verbose = FALSE)
  if (!is.null(save_dir)){
    save_plot <- paste0(save_dir, "images/UMAP_", combined_name, ".pdf")
    col_by_list <- c("orig.ident", "cluster")
    plotDimRed(sample_object = merged_seurat, save_plot = save_plot,
        col_by = col_by_list, return_plot = FALSE)
  }
  return(merged_seurat)
}

integrate_data <- function(seurat_1, seurat_2, intersecting_regions,
  save_dir = NULL, combined_name = "integrated"){
  # find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = list(seurat_1, seurat_2),
    anchor.features = intersecting_regions[[1]],
    assay = c("overlapPeaks", "overlapPeaks"),
    k.filter = NA
  )
  integrated_seurat <- IntegrateData(
    anchorset = anchors,
    weight.reduction = seurat_1[['lsi']],
    preserve.order = TRUE
  )
  integrated_seurat <- RunSVD(object = integrated_seurat,
    n = 30, reduction.name = 'integratedLSI')
  integrated_seurat <- RunUMAP(object = integrated_seurat,
    dims = 1:30, reduction = "integratedLSI",
    umap.method = "umap-learn", metric = "correlation")
  integrated_seurat <- FindNeighbors(object = integrated_seurat,
    reduction = 'integratedLSI', dims = 1:30)
    integrated_seurat <- FindClusters(object = integrated_seurat, verbose = FALSE)
  if (!is.null(save_dir)){
      save_plot <- paste0(save_dir, "images/UMAP_", combined_name, ".pdf")
      col_by_list <- c("orig.ident", "cluster")
      plotDimRed(sample_object = merged_seurat, save_plot = save_plot,
          col_by = col_by_list, return_plot = FALSE)
  }
}

harmonize_data <- function(seurat_1, seurat_2, save_dir = NULL, combined_name = "merged",
  assay_1 = "overlapPeaks", assay_2 = "overlapPeaks", merge_with_regions = FALSE){

  if (merge_with_regions){
    merged_seurat <- MergeWithRegions(
      object.1 = seurat_1,
      object.2 = seurat_2,
      sep.1 = c("-", "-"),
      sep.2 = c("-", "-"),
      assay.1 = assay_1,
      assay.2 = assay_2
      )
      assay_use <- "peaks"
    } else {
      DefaultAssay(object = seurat_1) <- "overlapPeaks"
      DefaultAssay(object = seurat_2) <- "overlapPeaks"

      merged_seurat <- merge(
        x = seurat_1,
        y = seurat_2,
        merge.data = FALSE
     )
      assay_use <- "overlapPeaks"
    }

  merged_seurat <- RunTFIDF(merged_seurat)
  merged_seurat <- FindTopFeatures(merged_seurat, min.cutoff = 50)
  merged_seurat <- RunSVD(merged_seurat, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
  merged_seurat <- RunUMAP(merged_seurat, reduction = 'lsi', dims = 1:30,
    umap.method = "umap-learn", metric = "correlation")
  merged_seurat <- RunHarmony(
    object = merged_seurat,
    group.by.vars = "orig.ident",
    reduction = "lsi",
    assay.use = assay_use,
    project.dim = FALSE)


  merged_seurat <- RunUMAP(merged_seurat, reduction = 'harmony', dims = 1:30,
    umap.method = "umap-learn", metric = "correlation")
  merged_seurat <- FindNeighbors(object = merged_seurat, reduction = 'harmony', dims = 1:30)
  merged_seurat <- FindClusters(object = merged_seurat, verbose = FALSE)
  if (!is.null(save_dir)){
    save_plot <- paste0(save_dir, "images/UMAP_", combined_name, ".pdf")
    col_by_list <- c("orig.ident", "cluster")
    plotDimRed(sample_object = merged_seurat, save_plot = save_plot,
        col_by = col_by_list, return_plot = FALSE)
  }
  return(merged_seurat)
}



################################## Pseudotime functions ######################

plot_sling_pseudotime <- function(seurat_object, y_val, col_by,
                                  pseudotime_curve, color = NULL,
                                  save_plot = NULL, range = NULL,
                                  plot_type = "dot_plot",
                                  height = 3, width = 7,
                                  ggrastr = FALSE,
                                  size = 0.25, raster_width = NULL) {
  print(save_plot)
  # pseudotime <- data.frame(slingshot::slingPseudotime(sling_object))
  # pseudotime_df <- data.frame(pseudotime = pseudotime[[pseudotime_curve]],
  #   row.names = rownames(pseudotime))
  # pseudotime_df <- pseudotime_df[!is.na(pseudotime_df$pseudotime), , drop = FALSE]
  plot_data <- make_plot_df(seurat_object = seurat_object, y_val = y_val,
                            x_val = pseudotime_curve, col_by = col_by)
  # plot_data <- merge(pseudotime_df, plot_data, by = "row.names", all = FALSE)
  if (is.null(color)){
    nColors <- length(levels(factor(plot_data$col_by)))
    color <- RColorBrewer::brewer.pal(nColors, "Set1")
  }
  print(head(plot_data))
  base_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_value,
                                                       y = y_value))
  if (plot_type == "dot_plot"){
    if(ggrastr){
      base_plot <- base_plot + ggrastr::geom_point_rast(ggplot2::aes(color = col_by),
        raster.width = raster_width, size = size) +
        ggplot2::geom_smooth(se = FALSE, color = "black")
        
    } else {
      base_plot <- base_plot + ggplot2::geom_point(ggplot2::aes(color = col_by)) +
        ggplot2::geom_smooth(se = FALSE, color = "black")
    }
  } else if (plot_type == "density"){
    base_plot <- base_plot + ggridges::geom_density_ridges(ggplot2::aes(fill = col_by)) +
      ggplot2::scale_fill_manual(values = color, name = col_by)
  } else {
    stop("plot_type must be dot_plot or density")
  }
  if (!(is.null(range))) {
    base_plot <- base_plot + ggplot2::xlim(range)
  }

  if (color == "viridis") {
    base_plot <- base_plot + viridis::scale_color_viridis(option = "plasma",
      name = "pseudotime")
  } else {
    base_plot <- base_plot + 
      ggplot2::scale_color_manual(values = color, name = col_by)
  }

  base_plot <- base_plot + ggplot2::ylab(y_val) + ggplot2::xlab("pseudotime")

  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot, height = height, width = width)
  }

  return(base_plot)
}


make_plot_df <- function(seurat_object, y_val, x_val,
                          col_by = NULL) {

  # Add y_value to a data frame used for plotting. This value can be a gene, ADT,
  # or a value from meta data like nGene
  if (y_val %in% rownames(seurat_object) |
    y_val %in% colnames(seurat_object[[]])){
    plot_data <- FetchData(object = seurat_object, vars = y_val)
    print("fetching data")
  }else if (y_val %in% rownames(seurat_object[["ADT"]])){
    plot_data <- FetchData(object = seurat_object, vars = paste0("adt_", y_val))
    print("fetching data from ADT slot")
  }else {
    stop("y_val must be a gene, metric from meta data or 'cluster'")
  }



  # Name the column
  names(plot_data) <- "y_value"
  
  # Add a column contining the x_value. This should be something discrete
  # Like timepoint or cluster
  if (x_val == "cluster" | x_val == "Cluster"){
    plot_data$x_value <- as.data.frame(Idents(object = seurat_object))
  }else if (x_val %in% rownames(seurat_object) |
    x_val %in% colnames(seurat_object[[]])){
    plot_data$x_value <- FetchData(object = seurat_object, vars = x_val)[, 1]
  }else if (x_val %in% rownames(seurat_object[["ADT"]])){
    plot_data$x_value <- FetchData(object = seurat_object,
      vars = paste0("adt_", x_val))[, 1]
  }else {
    stop("x_val must be a gene, metric from meta data or 'cluster'")
  }


  # Determine how to color the plot. Default is the x_value but can be any
  # discrete value.
  if (is.null(col_by)) {
    plot_data$col_by <- plot_data$x_value
  } else if (col_by %in% colnames(seurat_object[[]])) {
    plot_data$col_by <- FetchData(object = seurat_object, vars = col_by)[,1]
  } else if (col_by == "cluster") {
    plot_data$col_by <- Idents(object = seurat_object)
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  colnames(plot_data) <- c("y_value", "x_value", "col_by")
  return(plot_data)
}

ATAC_plot_pseudotime <- function(ArchRProj, archr_matrix, gene, trajectory,
  col_by, color = NULL, impute = FALSE, imputeWeights = NULL,
  save_plot = NULL, range = NULL, plot_type = "dot_plot",
  height = 3, width = 7, ggrastr = FALSE, raster_width = NULL,
  useSeqnames = "z", size = 0.25){
  dfT <- getCellColData(ArchRProj, select = c(trajectory, col_by))
  names(dfT) <- c("pseudotime", "col_by")
  dfT <- dfT[!is.na(dfT[,1]), , drop = FALSE]
  if(archr_matrix == "MotifMatrix"){
    ATAC_df <- getMatrixFromProject(ArchRProj, useMatrix = archr_matrix,
      useSeqnames = useSeqnames)
    counts <- assays(ATAC_df)[[useSeqnames]]
    rownames(counts) <- rowData(ATAC_df)$name
    rownames(counts) <- paste0(useSeqnames, ":", rownames(counts))
  } else {
    ATAC_df <- getMatrixFromProject(ArchRProj, useMatrix = archr_matrix)
    counts <- assays(ATAC_df)[[archr_matrix]]
    rownames(counts) <- rowData(ATAC_df)$name
  }
  gene_counts <- counts[gene, ]
  if(impute){
    gene_counts <- matrix(gene_counts, nrow = 1)
    colnames(gene_counts) <- colnames(counts)
    gene_counts <- imputeMatrix(mat = gene_counts,
      imputeWeights = imputeWeights)
    gene_counts <- data.frame(gene_counts)
  } else {
      gene_counts <- t(gene_counts)
  }
  plot_data <- merge(dfT, gene_counts, by = "row.names", all = FALSE)
  plot_data <- as.data.frame(plot_data)
  print(head(plot_data))
  if (is.null(color)){
    nColors <- length(levels(factor(plot_data$col_by)))
    color <- RColorBrewer::brewer.pal(nColors, "Set1")
  }
  base_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pseudotime,
                                                       y = gene_counts))
  if (plot_type == "dot_plot"){
    if(ggrastr){
      base_plot <- base_plot + ggrastr::geom_point_rast(ggplot2::aes(color = col_by),
        raster.width = raster_width, size = size) +
        ggplot2::geom_smooth(se = FALSE, color = "black") 
    } else {
      base_plot <- base_plot + ggplot2::geom_point(ggplot2::aes(color = col_by)) +
        ggplot2::geom_smooth(se = FALSE, color = "black") 
    }
  } else if (plot_type == "density"){
    base_plot <- base_plot + ggridges::geom_density_ridges(ggplot2::aes(fill = col_by)) +
      ggplot2::scale_fill_manual(values = color, name = col_by)
  } else {
    stop("plot_type must be dot_plot or density")
  }
  if (!(is.null(range))) {
    base_plot <- base_plot + ggplot2::xlim(range)
  }
  if (color == "viridis") {
    base_plot <- base_plot + viridis::scale_color_viridis(option = "plasma",
      name = "pseudotime")
  } else {
    base_plot <- base_plot + 
      ggplot2::scale_color_manual(values = color, name = col_by)
  }
  base_plot <- base_plot + ggplot2::ylab(gene)
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot, height = height, width = width)
  }

  return(base_plot)
}



getTrajectoryRNA <- function(
  seurat_object = NULL,
  name = "Trajectory",
  groupEvery = 1,
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = 11,
  groupList = NULL
  ){


  # pull out trajectory values
  trajectory <- seurat_object[[name]]
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]

  breaks <- seq(0, 100, groupEvery)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }


  if(is.null(groupList)){
    groupList <- lapply(seq_along(breaks), function(x){
        if(x == 1){
            NULL
        }else{
            cells <- rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] &
              trajectory[,1] <= breaks[x])]
            # if(breaks[x-1] == 0){
            #     zero_cells <- rownames(trajectory)[which(trajectory[,1] == 0)]
            #     cells <- c(cells, zero_cells)
            # }
            return(cells)
        }
    })[-1]

    names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  }

  assay_data <- GetAssayData(object = seurat_pseudotime, slot = "data")
  message("Creating Trajectory Group Matrix..")

  mat_list <- lapply(names(groupList), function(x){
    cells <- groupList[[x]]
    subset_assay <- assay_data[ , cells, drop = FALSE]
    if(ncol(subset_assay) >= 2){
      new_assay <- data.frame(rowSums(subset_assay))
      rownames(new_assay) <- rownames(subset_assay)
    } else if(ncol(subset_assay) == 1){
      new_assay <- subset_assay
    } else {
      new_assay <- data.frame(rep(0, ncol(subset_assay)))
      rownames(new_assay) <- rownames(subset_assay)
    }
    colnames(new_assay) <- x
    return(new_assay)
  })

  groupMat <- do.call(cbind, mat_list)

  groupMat <- groupMat[ , colSums(groupMat) > 0]


  #Scale
  if(!is.null(scaleTo)){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }else{
      # If there are any points in pseudotime that don't have cells, 
      # this adds in some NA values that then throw off the next section
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }

  if(log2Norm){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }else{
      groupMat <- log2(groupMat + 1)
    }
  }

  featureDF <- data.frame(name = rownames(groupMat))

  if(!is.null(smoothWindow)){
    
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) ArchR:::.centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          smoothMat = as.matrix(smoothGroupMat), 
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )

  }else{

    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }

  useMatrix <- "RNA"
  matrixClass <- "RNA_matrix"
  metadata(seTrajectory)$Params <- list(
    useMatrix = useMatrix, 
    matrixClass = matrixClass,
    scaleTo = scaleTo, 
    log2Norm = log2Norm, 
    smoothWindow = smoothWindow, 
    date = Sys.Date()
  )

  return(seTrajectory)

}


getTrajectory_Kristen <- function(
  ArchRProj = NULL,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  groupEvery = 1,
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = 11,
  threads = getArchRThreads(),
  groupList = NULL
  ){

  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = groupEvery, name = "groupEvery", valid = c("numeric"))
  ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  ArchR:::.validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))

  trajectory <- getCellColData(ArchRProj, name)
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(0, 100, groupEvery)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }

  if(is.null(groupList)){
    groupList <- lapply(seq_along(breaks), function(x){
        if(x == 1){
            NULL
        }else{
            cells <- rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] &
              trajectory[,1] <= breaks[x])]
            # if(breaks[x-1] == 0){
            #     zero_cells <- rownames(trajectory)[which(trajectory[,1] == 0)]
            #     cells <- c(cells, zero_cells)
            # }
            return(cells)
        }
    })[-1]
    names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  }

  featureDF <- ArchR:::.getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))

  message("Creating Trajectory Group Matrix..")
  groupMat <- ArchR:::.getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = featureDF,
      groupList = groupList, 
      threads = threads, 
      verbose = FALSE, 
      useMatrix = useMatrix
  )

  # # Added by Kristen - Remove columns that have sums of 0 no cells
  #groupMat <- groupMat[ , colSums(groupMat) > 0]
  groupMat <- groupMat[ , !(lengths(groupList) == 0)]
  #Scale
  if(!is.null(scaleTo)){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }else{
      # If there are any points in pseudotime that don't have cells, 
      # this adds in some NA values that then throw off the next section
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo

      # Added by Kristen - replace any NAs with 0
      groupMat[is.na(groupMat)] <- 0
    }
  }

  if(log2Norm){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }else{
      groupMat <- log2(groupMat + 1)
    }
  }

  if(!is.null(smoothWindow)){
    
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) ArchR:::.centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          smoothMat = as.matrix(smoothGroupMat), 
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }else{

    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
        assays = SimpleList(
          mat = as.matrix(groupMat)
        ), 
        rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }

  }

  metadata(seTrajectory)$Params <- list(
    useMatrix = useMatrix, 
    matrixClass = matrixClass,
    scaleTo = scaleTo, 
    log2Norm = log2Norm, 
    smoothWindow = smoothWindow, 
    date = Sys.Date()
  )

  return(seTrajectory)

}

# Remake plots with new names
draw_heatmap <- function(max_ident, mat, pal, seTrajectory){
  colData <- data.frame(max_ident)
  names(colData) <- "cell_type"
  ht <- ArchR:::.ArchRHeatmap(mat = mat, scale = FALSE, limits = c(min(mat),
    max(mat)), color = pal, clusterCols = FALSE, clusterRows = FALSE,
    labelRows = FALSE, labelCols = FALSE, showColDendrogram = TRUE,
    name = metadata(seTrajectory)$Params$useMatrix, draw = FALSE,
    colorMap = colorMap, colData = colData)
  return(ht)
}

plotTrajectoryHeatmap_kristen <- function (seTrajectory = NULL,
  varCutOff = 0.9, maxFeatures = 25000, scaleRows = TRUE,
  limits = c(-1.5, 1.5), grepExclude = NULL, pal = NULL, 
  labelMarkers = NULL, labelTop = 50, labelRows = FALSE,
  rowOrder = NULL, useSeqnames = NULL, returnMat = FALSE,
  force = FALSE, max_ident = NULL, colorMap = colorMap,
  colData = NULL) 
{
    ArchR:::.validInput(input = seTrajectory, name = "seTrajectory", 
        valid = c("SummarizedExperiment"))
    ArchR:::.validInput(input = varCutOff, name = "varCutOff", valid = c("numeric", 
        "null"))
    ArchR:::.validInput(input = maxFeatures, name = "maxFeatures", valid = c("integer", 
        "null"))
    ArchR:::.validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
    ArchR:::.validInput(input = limits, name = "limits", valid = c("numeric"))
    ArchR:::.validInput(input = grepExclude, name = "grepExclude", valid = c("character", 
        "null"))
    ArchR:::.validInput(input = pal, name = "pal", valid = c("palette", 
        "null"))
    ArchR:::.validInput(input = labelMarkers, name = "labelMarkers", 
        valid = c("character", "null"))
    ArchR:::.validInput(input = labelTop, name = "labelTop", valid = c("integer"))
    ArchR:::.validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
    ArchR:::.validInput(input = rowOrder, name = "rowOrder", valid = c("vector", 
        "null"))
    ArchR:::.validInput(input = useSeqnames, name = "useSeqnames", valid = c("character", 
        "null"))
    ArchR:::.validInput(input = returnMat, name = "returnMat", valid = c("boolean"))
    ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
    if (metadata(seTrajectory)$Params$matrixClass == "Sparse.Assays.Matrix") {
        if (is.null(useSeqnames) || length(useSeqnames) > 1) {
            message("useSeqnames is NULL or greater than 1 with a Sparse.Assays.Matrix trajectory input.")
            if (force) {
                message("force=TRUE thus continuing")
            } else {
                useSeqnames <- rev(unique(rowData(seTrajectory)$seqnames))[1]
                message(paste0("force=FALSE thus continuing with subsetting useSeqnames = ", 
                  useSeqnames))
            }
        }
    }
    if(is.null(colData)){
      colData <- data.frame(max_ident)
      names(colData) <- "cell_type"
    }
    if (!is.null(useSeqnames)) {
        seTrajectory <- seTrajectory[paste0(rowData(seTrajectory)$seqnames) %in% 
            paste0(useSeqnames), ]
    }
    if (nrow(seTrajectory) == 0) {
        stop("No features left in seTrajectory, please check input!")
    }
    mat <- assay(seTrajectory)
    rSNA <- rowSums(is.na(mat))
    if (sum(rSNA > 0) > 0) {
        message("Removing rows with NA values...")
        mat <- mat[rSNA == 0, ]
    }
    varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))
    orderedVar <- FALSE
    if (is.null(rowOrder)) {
        mat <- mat[order(varQ, decreasing = TRUE), ]
        orderedVar <- TRUE
        if (is.null(varCutOff) & is.null(maxFeatures)) {
            n <- nrow(mat)
        } else if (is.null(varCutOff)) {
            n <- maxFeatures
        } else if (is.null(maxFeatures)) {
            n <- (1 - varCutOff) * nrow(mat)
        } else {
            n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
        }
        n <- min(n, nrow(mat))
        mat <- mat[head(seq_len(nrow(mat)), n), ]
    }
    if (!is.null(labelTop)) {
        if (orderedVar) {
            idxLabel <- rownames(mat)[seq_len(labelTop)]
        } else {
            idxLabel <- rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
        }
    } else {
        idxLabel <- NULL
    }
    if (!is.null(labelMarkers)) {
        idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)), 
            nomatch = 0)
        idxLabel2 <- idxLabel2[idxLabel2 > 0]
    } else {
        idxLabel2 <- NULL
    }
    idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])
    if (scaleRows) {
        mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
            `/`)
        mat[mat > max(limits)] <- max(limits)
        mat[mat < min(limits)] <- min(limits)
    }
    if (nrow(mat) == 0) {
        stop("No Features Remaining!")
    }
    if (is.null(pal)) {
        if (is.null(metadata(seTrajectory)$Params$useMatrix)) {
            pal <- paletteContinuous(set = "solarExtra", n = 100)
        } else if (tolower(metadata(seTrajectory)$Params$useMatrix) == 
            "genescorematrix") {
            pal <- paletteContinuous(set = "blueYellow", n = 100)
        } else {
            pal <- paletteContinuous(set = "solarExtra", n = 100)
        }
    }
    if (!is.null(rowOrder)) {
        idx <- rowOrder
    } else {
        idx <- order(apply(mat, 1, which.max))
    }
    ht <- tryCatch({
        ArchR:::.ArchRHeatmap(mat = mat[idx, ], scale = FALSE, limits = c(min(mat), 
            max(mat)), color = pal, clusterCols = FALSE, clusterRows = FALSE, 
            labelRows = labelRows, labelCols = FALSE, customRowLabel = match(idxLabel, 
                rownames(mat[idx, ])), showColDendrogram = TRUE, 
            name = metadata(seTrajectory)$Params$useMatrix, draw = FALSE,
            colorMap = colorMap, colData = colData)
    }, error = function(e) {
        errorList = list(mat = mat[idx, ], scale = FALSE, limits = c(min(mat), 
            max(mat)), color = pal, clusterCols = FALSE, clusterRows = FALSE, 
            labelRows = labelRows, labelCols = FALSE, customRowLabel = match(idxLabel, 
                rownames(mat[idx, ])), showColDendrogram = TRUE, 
            name = metadata(seTrajectory)$Params$useMatrix, draw = FALSE)
        stop("Heatmap error")
    })
    if (returnMat) {
        return(mat[idx, ])
    } else {
        return(ht)
    }
}


findCorrelations <- function(ArchRProj, atac_cluster, cell_type,
  seZ, useMatrix = "GeneScoreMatrix", positive_TF_cor = 0.5,
  group_name = NULL){
  if(is.null(group_name)){
    group_name <- cell_type
  }
  # Pull out cells just in the one cluster
  col_data <- getCellColData(ArchRProj = ArchRProj)
  cell_idx <- col_data[[atac_cluster]] %in% cell_type

  cluster_cells <- ArchRProj$cellNames[cell_idx]

  if(length(cluster_cells) >= 150){



      projTonsils_cluster <- subsetCells(ArchRProj = ArchRProj,
        cellNames = cluster_cells)




      # Identify correlated motifs
      corGM_MM <- correlateMatrices(
          ArchRProj = projTonsils_cluster,
          useMatrix1 = useMatrix,
          useMatrix2 = "MotifMatrix",
          reducedDims = "Harmony"
      )

      # Add Maximum Delta Deviation to the Correlation Data Frame
      corGM_MM$maxDelta <- rowData(seZ)[match(corGM_MM$MotifMatrix_name, 
        rowData(seZ)$name), "maxDelta"]

      corGM_MM$cell_type <- group_name

      #######
      # Identify positive TF regulators
      corGM_MM <- corGM_MM[order(abs(corGM_MM$cor), decreasing = TRUE), ]
      corGM_MM <- corGM_MM[which(!duplicated(gsub("\\-.*","",
        corGM_MM[,"MotifMatrix_name"]))), ]
      corGM_MM$TFRegulator <- "NO"
      corGM_MM$TFRegulator[which(corGM_MM$cor > positive_TF_cor & 
        corGM_MM$padj < 0.01 & corGM_MM$maxDelta > 
        quantile(corGM_MM$maxDelta, 0.75))] <- "YES"
      corGM_MM$TFpos_neg <- corGM_MM$TFRegulator
      # Anti regulators
      corGM_MM$TFpos_neg[which(corGM_MM$cor < -positive_TF_cor &
        corGM_MM$padj < 0.01 & corGM_MM$maxDelta >
        quantile(corGM_MM$maxDelta, 0.75))] <- "ANTI"

      return(corGM_MM)
  } else {
    return(NULL)
  }
}

plot_correlations <- function(ArchRProj, atac_cluster, useMatrix,
  cell_type_levels, seZ, seurat_object, rna_cluster){
  cell_types <- getCellColData(ArchRProj, select = atac_cluster)
  cell_type_list <- lapply(unique(cell_types[ , 1]), function(x){
    print(x)
    cor_df <- findCorrelations(ArchRProj = ArchRProj, atac_cluster = atac_cluster,
      cell_type = x, useMatrix = useMatrix, seZ = seZ)
    return(cor_df)
    })

  # cell_type_list <- cell_type_list[-which(sapply(cell_type_list, is.null))]

  cell_type_df <- do.call(rbind, cell_type_list)
  print(head(cell_type_df))
  TF_regs <- unique(cell_type_df[
    cell_type_df$TFRegulator == "YES", ]$MotifMatrix_matchName)

  tf_reg_df <- cell_type_df[cell_type_df$MotifMatrix_matchName %in% TF_regs, ]

  tf_reg_df <- data.frame(tf_reg_df)
  print(head(tf_reg_df))

  txn_atac_ordering <- lapply(TF_regs, function(x){
    plot_data_genes <- tf_reg_df[tf_reg_df$MotifMatrix_matchName == x, ]
    return(plot_data_genes[plot_data_genes$padj == 
      min(plot_data_genes$padj), ]$cell_type)
    })

  txn_atac_ordering <- unlist(txn_atac_ordering)

  txn_atac_ordering <- as.character(txn_atac_ordering)

  names(txn_atac_ordering) <- TF_regs

  # Make the order of genes the same as the order of cell type in the plot
  txn_atac_ordering <- factor(txn_atac_ordering,
    levels = cell_type_levels)
  txn_atac_ordering <- txn_atac_ordering[order(txn_atac_ordering, decreasing = TRUE)]

  tf_reg_df$cell_type <- factor(tf_reg_df$cell_type,
    levels = cell_type_levels)
  tf_reg_df$MotifMatrix_matchName <- factor(tf_reg_df$MotifMatrix_matchName,
    levels = names(txn_atac_ordering))
  return(tf_reg_df)

}

plot_RNA_expression <- function(seurat_object, gene_list, cluster_name){
  gene_list <- unique(gene_list[gene_list %in% rownames(seurat_object)])

  txn_plot <- DotPlot(object = seurat_object, features = gene_list,
      group.by = cluster_name,
      cols = c("lightgrey", "blue"))

  txn_plot_data <- txn_plot$data
  txn_plot_data <- txn_plot_data[!is.na(txn_plot_data$features.plot),]


  # only plot genes expressed in at least 25% of one cell type
  min_RNA_exp <- 25

  txn_genes_plot <- lapply(gene_list, function(x){
    plot_data_genes <- txn_plot_data[txn_plot_data$features.plot == x,]
    if(max(plot_data_genes$pct.exp) > min_RNA_exp){
      return(TRUE)
    } else {
      return(FALSE)
    }
    })

  txn_genes_plot <- unlist(txn_genes_plot)

  txn_regulators_filtered <- gene_list[txn_genes_plot]

  # Order genes based on in what cell type they are most highly expressed
  txn_ordering <- lapply(txn_regulators_filtered, function(x){
    plot_data_genes <- txn_plot_data[txn_plot_data$features.plot == x, ]
    return(plot_data_genes[plot_data_genes$pct.exp == 
      max(plot_data_genes$pct.exp), ]$id)
    })

  txn_ordering <- unlist(txn_ordering)

  txn_ordering <- as.character(txn_ordering)

  names(txn_ordering) <- txn_regulators_filtered

  # Make the order of genes the same as the order of cell type in the plot
  txn_ordering <- factor(txn_ordering,
    levels = levels(seurat_object[[]][[cluster_name]]))
  txn_ordering <- txn_ordering[order(txn_ordering)]

  txn_plot_final <- DotPlot(object = seurat_object, features = names(txn_ordering),
      group.by = cluster_name,
      cols = c("lightgrey", "blue"))

  return(txn_plot_final)


}

plot_correlations_geomtile <- function(ArchRProj, atac_cluster, useMatrix,
  cell_type_levels, seZ, remove_cell_type = NULL, positive_TF_cor = 0.5){
  seGroupMotif <- getGroupSE(ArchRProj = ArchRProj,
  useMatrix = "MotifMatrix", groupBy = atac_cluster)

  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs

  cell_types <- getCellColData(ArchRProj, select = atac_cluster)
  cell_type_list <- lapply(unique(cell_types[ , 1]), function(x){
    print(x)
    cor_df <- findCorrelations(ArchRProj = ArchRProj, atac_cluster = atac_cluster,
      cell_type = x, useMatrix = useMatrix, positive_TF_cor = positive_TF_cor,
      seZ = seZ)
    return(cor_df)
    })

  # cell_type_list <- cell_type_list[-which(sapply(cell_type_list, is.null))]

  cell_type_df <- do.call(rbind, cell_type_list)
  if(!is.null(remove_cell_type)){
      cell_type_df <- cell_type_df[!cell_type_df$cell_type %in% remove_cell_type, ]
  }
  print(head(cell_type_df))
  TF_regs <- unique(cell_type_df[
    cell_type_df$TFRegulator == "YES", ]$MotifMatrix_matchName)

  tf_reg_df <- cell_type_df[cell_type_df$MotifMatrix_matchName %in% TF_regs, ]

  tf_reg_df <- data.frame(tf_reg_df)
  print(head(tf_reg_df))
  tf_reg_df$new_cor <- "not_significant"
  tf_reg_df$new_cor[tf_reg_df$cor >= positive_TF_cor] <- "positive_cor"
  tf_reg_df$new_cor[tf_reg_df$cor <= -positive_TF_cor] <- "negative_cor"
  tf_reg_df[tf_reg_df$padj >= 0.05, ]$new_cor <- "not_significant"
  colors <- c("positive_cor" = "blue", "negative_cor" = "red",
    "not_significant" = "black")
  tf_reg_df$new_cor <- factor(tf_reg_df$new_cor, levels = names(colors))


  txn_atac_ordering <- lapply(TF_regs, function(x){
    plot_data_genes <- tf_reg_df[tf_reg_df$MotifMatrix_matchName == x, ]
    return(plot_data_genes[plot_data_genes$cor == 
      max(plot_data_genes$cor), ]$cell_type)
    })

  txn_atac_ordering <- unlist(txn_atac_ordering)

  txn_atac_ordering <- as.character(txn_atac_ordering)

  names(txn_atac_ordering) <- TF_regs

  # Make the order of genes the same as the order of cell type in the plot
  txn_atac_ordering <- factor(txn_atac_ordering,
    levels = cell_type_levels)
  txn_atac_ordering <- txn_atac_ordering[order(txn_atac_ordering, decreasing = TRUE)]

  tf_reg_df$cell_type <- factor(tf_reg_df$cell_type,
    levels = cell_type_levels)
  tf_reg_df$MotifMatrix_matchName <- factor(tf_reg_df$MotifMatrix_matchName,
    levels = names(txn_atac_ordering))

  plot <- ggplot(data = tf_reg_df, aes(x = cell_type, y = MotifMatrix_matchName,
    fill = cor))+
    geom_tile() +
    # scale_fill_manual(values = colors)
    scale_fill_gradientn(colours = c("blue", "white", "white", "red"),
      values = c(0, 0.2, 0.5, 0.8, 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  return(plot)

}



# From ArchR
