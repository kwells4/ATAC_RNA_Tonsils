library(ArchR)
set.seed(1)

addArchRGenome("hg38")

output_dir2 <- "All_Tonsils_pre_peaks"
output_dir <- "All_Tonsils_new_chromvar"
#output_dir <- "All_Tonsils_new_chromvar_confident_peaks"
print(output_dir)
if(output_dir == "All_Tonsils_new_chromvar") {
    peak_reproducability <- "2"
} else if (output_dir == "All_Tonsils_new_chromvar_confident_peaks"){
    peak_reproducability <- "(n+1)/2"
}

seurat_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/output/allSamples_nomt_snakemake/files/"

markerGenes <- c(
    "TOP2A", "PCNA", "MKI67", "ICAM1", "SEMA7A", "BACH2", "AICDA", #Germinal center dz
    "MS4A1", "STX7", "HMMR", "PLK1", "POU2F2", "CXCR4", # germinal center lz
    "XBP1", "PRDM1", # plasmablasts
    "CD27", "AIM2", "SCIMP", # MBC
    "CXCR6", "FCRL4", "NEAT1", # FCRL4 MBC
    "CD72", "TCL1A", # Naive B cells
    "NFKB1", "CCND2", "IRF4", "CD83", "EGR3", "TCF4", # preGC/Acrivated B
    "CCR7", "CXCL13", "FKBP5", "LEF1", "KLF2", "KLF6", "RORA", "IL7R", "NKG7",
    "ANXA1", "TXNIP", "GZMA", "CLDND1", "BANK1", "ANXA2", "XCL1",
    "CD81", "CD8A", "CD4", "CD3D", "CD3G", # T cell
    "CLEC4C", "CD14", "FCGR2A", "IFNGR1" # DC
    )

# Remove the doublet cluster
projTonsils2 <- loadArchRProject("All_Tonsils_unanalyzed")
idx_doublet_cluster <- projTonsils2$Clusters_Harmony != "C7"

cells_doublet_cluster <- projTonsils2$cellNames[idx_doublet_cluster]

projTonsils3 <- subsetArchRProject(ArchRProj = projTonsils2,
  cells = cells_doublet_cluster, outputDirectory = output_dir)

output_dir <- getOutputDirectory(projTonsils3)

# Repeat analysis without doublet cluster
projTonsils3 <- addIterativeLSI(
    ArchRProj = projTonsils3,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)


# Batch correction
projTonsils3 <- addHarmony(
    ArchRProj = projTonsils3,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

# Find clusters 
projTonsils3 <- addClusters(
    input = projTonsils3,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
)


projTonsils3 <- addClusters(
    input = projTonsils3,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters_Harmony",
    resolution = 0.8,
    force = TRUE
)

table(projTonsils3$Clusters)

cM <- confusionMatrix(paste0(projTonsils3$Clusters), paste0(projTonsils3$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

cMH <- confusionMatrix(paste0(projTonsils3$Clusters_Harmony),
  paste0(projTonsils3$Sample))

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
projTonsils3 <- addUMAP(
    ArchRProj = projTonsils3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.3, 
    metric = "cosine",
    force = TRUE
)

p1 <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projTonsils3, addDOC = FALSE, width = 5, height = 5)


# UMAP on batch correction
projTonsils3 <- addUMAP(
    ArchRProj = projTonsils3, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

p3 <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData",
  name = "Sample", embedding = "UMAP_Harmony")

p4 <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData",
  name = "Clusters_Harmony", embedding = "UMAP_Harmony")

plotPDF(p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf",
  ArchRProj = projTonsils3, addDOC = FALSE, width = 5, height = 5)


doublet_plot <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData",
  name = "DoubletScore", embedding = "UMAP_Harmony")

pdf(paste0(output_dir, "/Plots/Plot-Doublet-score-UMAP.pdf"))
doublet_plot

dev.off()

doublet_enrich_plot <- plotEmbedding(ArchRProj = projTonsils3, colorBy = "cellColData",
  name = "DoubletEnrichment", embedding = "UMAP_Harmony")

pdf(paste0(output_dir, "/Plots/Plot-Doublet-enrichment-UMAP.pdf"))
doublet_enrich_plot

dev.off()


# Find markers with Harmony
markersGS_harmony <- getMarkerFeatures(
    ArchRProj = projTonsils3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters_Harmony",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


markerListHarmony <- getMarkers(markersGS_harmony,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

dir.create(paste0(output_dir, "/files"))

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
    ArchRProj = projTonsils3, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils3, 
    addDOC = FALSE, width = 5, height = 5)


# Impuation with magic

projTonsils3 <- addImputeWeights(projTonsils3)

p2 <- plotEmbedding(
    ArchRProj = projTonsils3, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils3)
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils3, 
    addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projTonsils3, outputDirectory = output_dir,
  load = FALSE)


set.seed(1)
#projTonsils2 <- loadArchRProject("Tonsil_Harmony2")
# Load in seurat_object and name cells
library(Seurat)
seurat_object <- readRDS(paste0(seurat_dir,
"allSamples_nomt_dubRm_namedClust_seurat.rda"))
DefaultAssay(seurat_object) <- "RNA"
seurat_object[["SCT"]] <- NULL
seurat_object[["ADT"]] <- NULL

# ArchR accepts "unmodified Seurat objects as input to the integration workflow"
# I don't need to normalize because jeff just wants counts, he performs the
# normalization on his own

# This is to fix one bug in the ArchR code
outDir <- getOutputDirectory(projTonsils3)

# Unsupervised linkage
projTonsils3 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = FALSE,
    groupRNA = "num_cell_type",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# Further refine using constrained analysis 

# Identify which cell types from scRNA are most abundant in each scATAC clusters
cM <- as.matrix(confusionMatrix(projTonsils3$Clusters_Harmony,
  projTonsils3$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

unique(unique(projTonsils3$predictedGroup_Un))

# cBcells <- paste0(c(paste0("0", 0:5), 12), collapse="|")

# cActiveB <- paste0(c("06", "07"), collapse = "|")

cBcells <- paste0(c(paste0("0", 0:6), 13), collapse="|")

cNonB <- paste0(paste0(14:22), collapse="|")

# cGCB <- paste0(c(paste0("0", 8:9), 10:11), collapse = "|")

cLight_GCB <- paste0(c(paste0("0", 7:8)), collapse = "|")

cDark_GCB <- paste0(c(paste0("0", 8:9), 10:12), collapse = "|")

clustBcells <- rownames(cM)[grep(cBcells, preClust)]

# clustActiveB <- rownames(cM)[grep(cActiveB, preClust)]

clustNonB <- rownames(cM)[grep(cNonB, preClust)]

# clustGCB <- rownames(cM)[grep(cGCB, preClust)]

clust_Dark_GCB <- rownames(cM)[grep(cDark_GCB, preClust)]

clust_Light_GCB <- rownames(cM)[grep(cLight_GCB, preClust)]

rnaBcells <- colnames(seurat_object)[grep(cBcells, seurat_object$num_cell_type)]

# rnaActiveB <- colnames(seurat_object)[grep(cActiveB, seurat_object$num_cell_type)]

rnaNonB <- colnames(seurat_object)[grep(cNonB, seurat_object$num_cell_type)]

# rnaGCB <- colnames(seurat_object)[grep(cGCB, seurat_object$num_cell_type)]

rna_Dark_GCB <- colnames(seurat_object)[grep(cDark_GCB, seurat_object$num_cell_type)]

rna_Light_GCB <- colnames(seurat_object)[grep(cLight_GCB, seurat_object$num_cell_type)]



groupList <- SimpleList(
    Bcells = SimpleList(
        ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clustBcells],
        RNA = rnaBcells
    ),
    nonB = SimpleList(
        ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clustNonB],
        RNA = rnaNonB
    ),
    # GCB = SimpleList(
    #     ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clustGCB],
    #     RNA = rnaGCB
    # ),
    # ActiveB = SimpleList(
    #     ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clustActiveB],
    #     RNA = rnaActiveB)
    LightGCB = SimpleList(
        ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clust_Light_GCB],
        RNA = rna_Light_GCB
    ),
    DarkGCB = SimpleList(
        ATAC = projTonsils3$cellNames[projTonsils3$Clusters_Harmony %in% clust_Dark_GCB],
        RNA = rna_Dark_GCB)    
)

projTonsils3 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "num_cell_type",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)


pal <- paletteDiscrete(values = seurat_object$num_cell_type)

p1 <- plotEmbedding(
    projTonsils3, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal,
    embedding = "UMAP_Harmony"
)

p2 <- plotEmbedding(
    projTonsils3, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal,
    embedding = "UMAP_Harmony"
)

plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projTonsils3,
  addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = projTonsils3, outputDirectory = output_dir,
  load = FALSE)

# Start here June 15, 2020
projTonsils3 <- loadArchRProject(output_dir)
set.seed(1)
projTonsils3 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "num_cell_type",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)


getAvailableMatrices(projTonsils3)

# Add pseudo-scRNA-seq profiles
projTonsils3 <- addImputeWeights(projTonsils3)

p1 <- plotEmbedding(
    ArchRProj = projTonsils3, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils3)
)


plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
    ArchRProj = projTonsils3, 
    addDOC = FALSE, width = 5, height = 5)


# Label scATAC-seq with scRNA-seq info
cM <- confusionMatrix(projTonsils3$Clusters_Harmony, projTonsils3$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

remapClust <- c("00_Naive_B" = "Naive_B",
                "01_Naive_B" = "Naive_B",
                "02_Naive_B" = "Naive_B",
                "03_Naive_B" = "Naive_B",
                "04_Naive_B" = "Naive_B",
                "05_Activated_B" = "Activated_B",
                "06_Activated_B" = "Activated_B",
                "07_Light_GC_B" = "Light_GC_B",
                "08_Light_GC_B" = "Light_GC_B",
                "09_Dark_GC_B" = "Dark_GC_B",
                "10_Dark_GC_B" = "Dark_GC_B",
                "11_Dark_GC_B" = "Dark_GC_B",
                "12_Dark_GC_B" = "Dark_GC_B",
                "13_Plasma_B" = "Plasma_B",
                "14_Plasma_B" = "Plasma_B",
                "15_CD4_T" = "CD4_T",
                "16_CD4_T" = "CD4_T",
                "17_CD4_T" = "CD4_T",
                "18_CD4_T" = "CD4_T",
                "19_CD4_T" = "CD4_T",
                "20_CD8_T" = "CD8_T",
                "21_Dendritic" = "Dendritic",
                "22_Monocyte_DC" = "Monocyte_DC")

labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust),
  newLabels = remapClust)
labelNew2

projTonsils3$Clusters2 <- mapLabels(projTonsils3$Clusters_Harmony,
  newLabels = labelNew2, oldLabels = labelOld)


cell_types <- c("Naive_B", "Activated_B", "Light_GC_B", "Dark_GC_B",
                "Plasma_B", "CD4_T", "CD8_T", "Monocyte_DC", "Dendritic")

colors_cell_type_new <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type_new) <- cell_types

colors_cell_type_new <- colors_cell_type_new[order(names(colors_cell_type_new))]

p1 <- plotEmbedding(projTonsils3, colorBy = "cellColData", name = "Clusters2",
    pal = colors_cell_type_new, embedding = "UMAP_Harmony")

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projTonsils3,
  addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = projTonsils3, outputDirectory = output_dir,
  load = FALSE)

saveArchRProject(ArchRProj = projTonsils3, outputDirectory = output_dir2,
  load = FALSE)


# Make pseudobulk replicates
projTonsils4 <- addGroupCoverages(ArchRProj = projTonsils3, groupBy = "Clusters2")


# Call peaks with Macs2
pathToMacs2 <- findMacs2()

projTonsils4 <- addReproduciblePeakSet(
    ArchRProj = projTonsils4, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2,
    reproducability = peak_reproducability
)
getPeakSet(projTonsils4)

peak_set <- getPeakSet(projTonsils4)

saveRDS(peak_set, paste0(outDir, "/files/macs2_peaks.rds"))

saveArchRProject(ArchRProj = projTonsils4, outputDirectory = output_dir,
  load = FALSE)

# Add peak matrix
projTonsils5 <- addPeakMatrix(projTonsils4)

saveArchRProject(ArchRProj = projTonsils4, outputDirectory = output_dir,
  load = FALSE)

projTonsils4 <- loadArchRProject(output_dir)
outDir <- getOutputDirectory(projTonsils4)

# Add peaks back because they are deleted when I save and quit
peakSet <- readRDS(paste0(outDir, "/files/macs2_peaks.rds"))

projTonsils4 <- addPeakSet(ArchRProj = projTonsils4, peakSet = peakSet,
  force = TRUE)

# Add peak matrix
projTonsils5 <- addPeakMatrix(projTonsils4)


getAvailableMatrices(projTonsils5)

table(projTonsils5$Clusters2)

# Find marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = projTonsils5, 
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
  ArchRProj = projTonsils5, addDOC = FALSE)


#############################
#        Chromvar           #
#############################

projTonsils5 <- addMotifAnnotations(ArchRProj = projTonsils5,
  motifSet = "cisbp", name = "Motif")

# Motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projTonsils5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Change n based on what is good
# Can also update this to be the values called by the correlation analysis
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 30, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 12,
  height = 6, ArchRProj = projTonsils5, addDOC = FALSE)

# Motif deviations

# Check if motif annotations added
if("Motif" %ni% names(projTonsils5@peakAnnotation)){
    projTonsils5 <- addMotifAnnotations(ArchRProj = projTonsils5,
      motifSet = "cisbp", name = "Motif")
}

projTonsils5 <- addBgdPeaks(projTonsils5)

# outDir <- getOutputDirectory(projTonsils5)
# projTonsils5 <- loadArchRProject("All_Tonsils")
# peakSet <- readRDS(paste0(outDir, "/files/macs2_peaks.rds"))

# projTonsils5 <- addPeakSet(ArchRProj = projTonsils5, peakSet = peakSet,
#   force = TRUE)

anno <- getPeakAnnotation(projTonsils5, "Motif")
matches <- readRDS(anno$Matches)

# This gene was previously returning "NA" values, so I'm removing it
matches <- matches[,colnames(matches) != "ENSG00000250542_156"]

projTonsils5 <- addDeviationsMatrix(
  ArchRProj = projTonsils5, 
  peakAnnotation = "Motif",
  matches = matches
)

saveArchRProject(ArchRProj = projTonsils5, outputDirectory = output_dir,
  load = FALSE)

projTonsils5 <- loadArchRProject(output_dir)

#####################
# I saved here and quit. Saving somehow removes the peak set
outDir <- getOutputDirectory(projTonsils5)

peakSet <- readRDS(paste0(outDir, "/files/macs2_peaks.rds"))

projTonsils <- addPeakSet(ArchRProj = projTonsils, peakSet = peakSet,
  force = TRUE)

getPeakSet(projTonsils5)

plotVarDev <- getVarDeviations(projTonsils5, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5,
  height = 5, ArchRProj = projTonsils5, addDOC = FALSE)

motifs <- markerGenes
markerMotifs <- getFeatures(projTonsils5, select = paste(motifs, collapse="|"),
  useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

projTonsils5 <- addImputeWeights(projTonsils5)

p <- plotGroups(ArchRProj = projTonsils5, 
  groupBy = "Clusters2", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projTonsils5)
)

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5,
  height = 5, ArchRProj = projTonsils5, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projTonsils5, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils5)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-Deviations-W-Imputation.pdf", 
    ArchRProj = projTonsils5, 
    addDOC = FALSE, width = 5, height = 5)

# Motif footprinting
motifPositions <- getPositions(projTonsils5)
motifPositions

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
  value = TRUE)))
markerMotifs

# Can make this for only a subset of cell types using the useGroups argument
seFoot <- getFootprints(
  ArchRProj = projTonsils5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)

# Subtract Tn5 bias
p1 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Divide by Tn5 bias
p2 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# No normalization
p3 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils5, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Co-accessibility
projTonsils5 <- addCoAccessibility(
    ArchRProj = projTonsils5,
    reducedDims = "Harmony"
)

cA <- getCoAccessibility(
    ArchRProj = projTonsils5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

p <- plotBrowserTrack(
    ArchRProj = projTonsils5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projTonsils5)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projTonsils5, 
    addDOC = FALSE, width = 5, height = 5)

# Peak to gene links
projTonsils5 <- addPeak2GeneLinks(
    ArchRProj = projTonsils5,
    reducedDims = "Harmony"
)


p <- plotBrowserTrack(
    ArchRProj = projTonsils5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projTonsils5)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = projTonsils5, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotPeak2GeneHeatmap(ArchRProj = projTonsils5, groupBy = "Clusters2", k = 20,
    palGroup = colors_cell_type_new)
output_dir <- getOutputDirectory(projTonsils5)
pdf(paste0(output_dir, "/Plots/Plot-Peak2Gene_Heatmap.pdf"), height = 10, width = 15)
p
dev.off()

# Find positive TF regulators
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

corGSM_MM

corGIM_MM <- correlateMatrices(
    ArchRProj = projTonsils5,
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
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])



saveArchRProject(ArchRProj = projTonsils5, outputDirectory = output_dir,
  load = FALSE)

