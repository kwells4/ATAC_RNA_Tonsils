library(ArchR)
set.seed(1)

inputFiles <- c(Tonsil1a = "fragment_files/scATAC_tonsil1a_fragments.tsv.gz",
	            Tonsil1b = "fragment_files/scATAC_tonsil1b_fragments.tsv.gz",
	            Tonsil2a = "fragment_files/scATAC_tonsil2a_fragments.tsv.gz",
	            Tonsil2b = "fragment_files/scATAC_tonsil2b_fragments.tsv.gz",
	            Tonsil3a = "fragment_files/scATAC_tonsil3a_fragments.tsv.gz",
	            Tonsil3b = "fragment_files/scATAC_tonsil3b_fragments.tsv.gz",
                BCP003 = "tonsil_fragments/BCP003_fragments.tsv.gz",
                BCP004 = "tonsil_fragments/BCP004_fragments.tsv.gz",
                BCP005 = "tonsil_fragments/BCP005_fragments.tsv.gz",
                BCP006 = "tonsil_fragments/BCP006_fragments.tsv.gz")

addArchRGenome("hg38")

output_dir <- "All_Tonsils_unanalyzed"

# Make arrow files and perform quality analysis
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# # Compute doublet scores
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)


ArrowFiles_all <- c("Tonsil1a.arrow", "Tonsil1b.arrow",
                    "Tonsil2a.arrow", "Tonsil2b.arrow",
                    "Tonsil3a.arrow", "Tonsil3b.arrow",
                    "BCP003.arrow", "BCP004.arrow",
                    "BCP005.arrow", "BCP006.arrow")

# Create ArchR project
projTonsils <- ArchRProject(
  ArrowFiles = ArrowFiles_all, 
  outputDirectory = output_dir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# Ridgeplot of TSS enrichment
p1 <- plotGroups(
    ArchRProj = projTonsils, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )


# Violin plot of TSS enrichment
p2 <- plotGroups(
    ArchRProj = projTonsils, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )


# Ridge plot of unique nuclear frags 
p3 <- plotGroups(
    ArchRProj = projTonsils, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

# Violin plot of unique nuclear frags
p4 <- plotGroups(
    ArchRProj = projTonsils, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projTonsils,
  addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = projTonsils)

p2 <- plotTSSEnrichment(ArchRProj = projTonsils)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf",
  ArchRProj = projTonsils, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projTonsils,
  outputDirectory = "All_Tonsils_unanalyzed", load = FALSE)

# Filter doublets 
projTonsils2 <- filterDoublets(projTonsils, filterRatio = 1.4)

# Run itterative LSI to do dimensionality reduction
projTonsils2 <- addIterativeLSI(
    ArchRProj = projTonsils2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:28
)



# Batch correction
projTonsils2 <- addHarmony(
    ArchRProj = projTonsils2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Find clusters 
projTonsils2 <- addClusters(
    input = projTonsils2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

projTonsils2 <- addClusters(
    input = projTonsils2,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters_Harmony",
    resolution = 0.8
)

table(projTonsils2$Clusters)


cM <- confusionMatrix(paste0(projTonsils2$Clusters), paste0(projTonsils2$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)


cMH <- confusionMatrix(paste0(projTonsils2$Clusters_Harmony),
  paste0(projTonsils2$Sample))

cMH <- cMH / Matrix::rowSums(cMH)
pH <- pheatmap::pheatmap(
    mat = as.matrix(cMH), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

pdf(paste0(output_dir, "/Plots/cluster_confusion_matrix1.pdf"))
p
dev.off()


pdf(paste0(output_dir, "/Plots/cluster_confusion_matrixH.pdf"))
pH
dev.off()


# UMAP on batch correction
projTonsils2 <- addUMAP(
    ArchRProj = projTonsils2, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.2, 
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Sample", embedding = "UMAP_Harmony")

p4 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Clusters_Harmony", embedding = "UMAP_Harmony")

plotPDF(p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf",
  ArchRProj = projTonsils2, addDOC = FALSE, width = 5, height = 5)


doublet_plot <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "DoubletScore", embedding = "UMAP_Harmony")

pdf(paste0(output_dir, "/Plots/Plot-Doublet-score-UMAP.pdf"))
doublet_plot

dev.off()

doublet_enrich_plot <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "DoubletEnrichment", embedding = "UMAP_Harmony")

pdf(paste0(output_dir, "/Plots/Plot-Doublet-enrichment-UMAP.pdf"))
doublet_enrich_plot

dev.off()


