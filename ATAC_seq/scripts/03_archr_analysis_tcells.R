library(ArchR)
set.seed(1)

addArchRGenome("hg38")

output_dir <- "T_cells"
output_dir2 <- "T_cells"

print(output_dir)

seurat_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/output/allSamples_nomt_TCells_snakemake/files/"



markerGenes <- c(
    "PRDM1", "RORA", "RGS1", "CYTOR", "FOXP3", "CXCL13", "CXCR5", "FKBP5",
    "PDCD1", "BCL6", "PASK", "IL7R", "NABP1", "ANK3", "ZFP36L2",
    "CCR7", "LEF1", "ITGA6", "KLF2", "SELL", "GZMA", "CCL5", "GZMK", "CST7",
    "PRF1", "ANXA1", "MKI67", "PCNA", "TOP2A", "GNLY", "XCL1",
    "CTSW", "ID2", "CD4", "CD8A"
    )




projTonsils <- loadArchRProject("All_Tonsils_unanalyzed")

# Only keep T cells
idx_t_cells <- projTonsils$Clusters_Harmony %in% c("C14", "C15", "C16", "C17")

cells_t_cells <- projTonsils$cellNames[idx_t_cells]

projTonsils2 <- subsetArchRProject(ArchRProj = projTonsils,
  cells = cells_t_cells, outputDirectory = output_dir)

getOutputDirectory(projTonsils2)

# Repeat analysis without doublet cluster
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
    dimsToUse = 1:25,
    force = TRUE
)


# Batch correction by sample
projTonsils2 <- addHarmony(
    ArchRProj = projTonsils2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

# Find clusters 
projTonsils2 <- addClusters(
    input = projTonsils2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
)


projTonsils2 <- addClusters(
    input = projTonsils2,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters_Harmony",
    resolution = 0.8,
    force = TRUE
)


table(projTonsils2$Clusters_Harmony)

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


# Run UMAP
projTonsils2 <- addUMAP(
    ArchRProj = projTonsils2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
    metric = "cosine",
    force = TRUE
)

p1 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
    name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projTonsils2,
    addDOC = FALSE, width = 5, height = 5)


# UMAP on batch correction
projTonsils2 <- addUMAP(
    ArchRProj = projTonsils2, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.1, 
    metric = "cosine",
    force = TRUE
)

p3 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Sample", embedding = "UMAP_Harmony")


p5 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Clusters_Harmony", embedding = "UMAP_Harmony")

plotPDF(p3,p5, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf",
  ArchRProj = projTonsils2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(projTonsils2, outputDirectory = output_dir)

# Find markers with Harmony
markersGS_harmony <- getMarkerFeatures(
    ArchRProj = projTonsils2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters_Harmony",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


markerListHarmony <- getMarkers(markersGS_harmony,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

lapply(names(markerListHarmony), function(x){
  marker_df <- markerListHarmony[[x]]
  write.table(marker_df, file = paste0(output_dir, "/files/",
    x, "_markers_harmony.txt"),
    sep = ",")
  })

heatmapGS_harmony <- markerHeatmap(
  seMarker = markersGS_harmony, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)


pdf(paste0(output_dir, "/Plots/marker_heatmap_genescore_harmony.pdf"),
  width = 8, height = 6)

heatmapGS_harmony

dev.off()

p2 <- plotEmbedding(
    ArchRProj = projTonsils2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)


# Impuation with magic

projTonsils2 <- addImputeWeights(projTonsils2)

p2 <- plotEmbedding(
    ArchRProj = projTonsils2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils2)
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)


set.seed(1)
#projTonsils2 <- loadArchRProject("Tonsil_Harmony2")
# Load in seurat_object and name cells
library(Seurat)
seurat_object <- readRDS(paste0(seurat_dir, "TCells_dubRm_namedClust_seurat.rda"))
DefaultAssay(seurat_object) <- "RNA"
seurat_object[["SCT"]] <- NULL
seurat_object[["ADT"]] <- NULL




seurat_object$tcell_cluster_char <- "00"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Tfh_CXCL13"] <- "01"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Tfh"] <- "02"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Naive_TCells"] <- "03"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Treg"] <- "04"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Tcm_CD4"] <- "05"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Circulating"] <- "06"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "Tcm_CD8"] <- "07"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "CTL"] <- "08"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "CTL_PRF1"] <- "09"
seurat_object$tcell_cluster_char[seurat_object$cell_type_tCell == "NK_cells"] <- "10"

seurat_object$num_cell_type_tCell <- paste0(seurat_object$tcell_cluster_char,
    "_", seurat_object$cell_type_tCell)

# ArchR accepts "unmodified Seurat objects as input to the integration workflow"
# I don't need to normalize because jeff just wants counts, he performs the
# normalization on his own

# This is to fix one bug in the ArchR code
outDir <- getOutputDirectory(projTonsils2)

# Unsupervised linkage
projTonsils2 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = FALSE,
    groupRNA = "num_cell_type_tCell",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# Further refine using constrained analysis 

# Identify which cell types from scRNA are most abundant in each scATAC clusters
cM <- as.matrix(confusionMatrix(projTonsils2$Clusters_Harmony,
  projTonsils2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

unique(unique(projTonsils2$predictedGroup_Un))

##################### Try 1
cCD8 <- paste0(paste0("0", c(7, 8, 9)), collapse = "|")

#cNK <- "10"

cGroup1 <- paste0(paste0("0", c(1:6)), collapse = "|")

cGroup1 <- paste0(c(cGroup1, "10"), collapse = "|")

# cGroup2 <- "06"


clustCD8 <- rownames(cM)[grep(cCD8, preClust)]

#clustNK <- rownames(cM)[grep(cNK, preClust)]

clustGroup1 <- rownames(cM)[grep(cGroup1, preClust)]

# clustGroup2 <- rownames(cM)[grep(cGroup2, preClust)]

rnaCD8 <- colnames(seurat_object)[grep(cCD8, seurat_object$num_cell_type_tCell)]

#rnaNK <- colnames(seurat_object)[grep(cNK, seurat_object$num_cell_type_tCell)]

rnaGroup1 <- colnames(seurat_object)[grep(cGroup1, seurat_object$num_cell_type_tCell)]

# rnaGroup2 <- colnames(seurat_object)[grep(cGroup2, seurat_object$num_cell_type_tCell)]


groupList <- SimpleList(
    CD8 = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustCD8],
        RNA = rnaCD8
    ),
    # NK = SimpleList(
    #     ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustNK],
    #     RNA = rnaNK
    # ),
    Group1 = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGroup1],
        RNA = rnaGroup1
    )
    # Group2 = SimpleList(
    #     ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGroup2],
    #     RNA = rnaGroup2)    
)


projTonsils2 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "num_cell_type_tCell",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)


pal <- paletteDiscrete(values = seurat_object$num_cell_type_tCell)

p1 <- plotEmbedding(
    projTonsils2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal,
    embedding = "UMAP_Harmony"
)

p2 <- plotEmbedding(
    projTonsils2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal,
    embedding = "UMAP_Harmony"
)

plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projTonsils2,
  addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)

#projTonsils3 <- loadArchRProject(output_dir)
set.seed(1)
projTonsils2 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "num_cell_type_tCell",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)


getAvailableMatrices(projTonsils2)

# Add pseudo-scRNA-seq profiles
projTonsils2 <- addImputeWeights(projTonsils2)

p1 <- plotEmbedding(
    ArchRProj = projTonsils2, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils2)
)


plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)


# Label scATAC-seq with scRNA-seq info
cM <- confusionMatrix(projTonsils2$Clusters_Harmony, projTonsils2$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

remapClust <- c("01_Tfh_CXCL13" = "Tfh",
                "03_Naive_TCells" = "Naive_TCells",
                "04_Treg" = "Treg",
                "05_Tcm_CD4" = "Tcm_CD4",
                "06_Circulating" = "Circulating",
                "07_Tcm_CD8" = "Tcm_CD8",
                "08_CTL" = "CTL")

labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust),
  newLabels = remapClust)
labelNew2

projTonsils2$Clusters2 <- mapLabels(projTonsils2$Clusters_Harmony,
  newLabels = labelNew2, oldLabels = labelOld)


cell_types <- levels(seurat_object$cell_type_tCell)

colors_cell_type_new <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type_new) <- cell_types

colors_cell_type_new <- colors_cell_type_new[order(names(colors_cell_type_new))]

colors_cell_type_new <- colors_cell_type_new[names(colors_cell_type_new) %in% labelNew2]

p1 <- plotEmbedding(projTonsils2, colorBy = "cellColData", name = "Clusters2",
    pal = colors_cell_type_new, embedding = "UMAP_Harmony")

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projTonsils2,
  addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)

pathToMacs2 <- findMacs2()

projTonsils2 <- addReproduciblePeakSet(
    ArchRProj = projTonsils2, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(projTonsils2)

peak_set <- getPeakSet(projTonsils2)

saveRDS(peak_set, paste0(outDir, "/files/macs2_peaks.rds"))

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)

# Add peak matrix
projTonsils2 <- addPeakMatrix(projTonsils2)

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)


getAvailableMatrices(projTonsils2)

table(projTonsils2$Clusters2)

# Find marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = projTonsils2, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

lapply(names(markerList), function(x){
  marker_df <- markerList[[x]]
  write.table(marker_df, file = paste0(outDir, "/files/", x, "_peak_markers.txt"),
    sep = ",")
  })

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6,
  ArchRProj = projTonsils2, addDOC = FALSE)


#############################
#        Chromvar           #
#############################

projTonsils2 <- addMotifAnnotations(ArchRProj = projTonsils2,
  motifSet = "cisbp", name = "Motif")

# Motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projTonsils2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Change n based on what is good
# Can also update this to be the values called by the correlation analysis
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 30, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 12,
  height = 6, ArchRProj = projTonsils2, addDOC = FALSE)

# Motif deviations

# Check if motif annotations added
if("Motif" %ni% names(projTonsils2@peakAnnotation)){
    projTonsils2 <- addMotifAnnotations(ArchRProj = projTonsils2,
      motifSet = "cisbp", name = "Motif")
}

projTonsils2 <- addBgdPeaks(projTonsils2)

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)



projTonsils2 <- addDeviationsMatrix(
  ArchRProj = projTonsils2, 
  peakAnnotation = "Motif"
)

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)


#####################
# I saved here and quit. Saving somehow removes the peak set
getPeakSet(projTonsils2)

plotVarDev <- getVarDeviations(projTonsils2, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5,
  height = 5, ArchRProj = projTonsils2, addDOC = FALSE)

motifs <- markerGenes
markerMotifs <- getFeatures(projTonsils2, select = paste(motifs, collapse="|"),
  useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

projTonsils2 <- addImputeWeights(projTonsils2)

p <- plotGroups(ArchRProj = projTonsils2, 
  groupBy = "Clusters2", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projTonsils2)
)

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5,
  height = 5, ArchRProj = projTonsils2, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projTonsils2, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils2)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-Deviations-W-Imputation.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)

# Motif footprinting
motifPositions <- getPositions(projTonsils2)
motifPositions

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
  value = TRUE)))
markerMotifs

# Can make this for only a subset of cell types using the useGroups argument
seFoot <- getFootprints(
  ArchRProj = projTonsils2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)

# Subtract Tn5 bias
p1 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils2, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Divide by Tn5 bias
p2 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils2, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# No normalization
p3 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils2, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Co-accessibility
projTonsils2 <- addCoAccessibility(
    ArchRProj = projTonsils2,
    reducedDims = "Harmony"
)

cA <- getCoAccessibility(
    ArchRProj = projTonsils2,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

p <- plotBrowserTrack(
    ArchRProj = projTonsils2, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projTonsils2)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)

# Peak to gene links
projTonsils2 <- addPeak2GeneLinks(
    ArchRProj = projTonsils2,
    reducedDims = "Harmony"
)


p <- plotBrowserTrack(
    ArchRProj = projTonsils2, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projTonsils2)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = projTonsils2, 
    addDOC = FALSE, width = 5, height = 5)

set.seed(0)
p <- plotPeak2GeneHeatmap(ArchRProj = projTonsils2, groupBy = "Clusters2", k = 20,
    palGroup = colors_cell_type_new)
output_dir <- getOutputDirectory(projTonsils2)
pdf(paste0(output_dir, "/Plots/Plot-Peak2Gene_Heatmap.pdf"), height = 10, width = 15)
p
dev.off()

# Find positive TF regulators
seGroupMotif <- getGroupSE(ArchRProj = projTonsils2,
  useMatrix = "MotifMatrix", groupBy = "Clusters2")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = projTonsils2,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGSM_MM

corGIM_MM <- correlateMatrices(
    ArchRProj = projTonsils2,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGIM_MM

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name,
  rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name,
  rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
tfReg_GSM_MM <- sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
write.table(tfReg_GSM_MM, paste0(output_dir, "/files/tf_regulator_GSM_MM.txt"))


corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
tfReg_GIM_MM <- sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
write.table(tfReg_GIM_MM, paste0(output_dir, "/files/tf_regulator_GIM_MM.txt"))



saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)

