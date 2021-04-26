library(ArchR)
set.seed(1)

addArchRGenome("hg38")

output_dir <- "B_cells"
output_dir2 <- "B_cells"

print(output_dir)

seurat_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/output/allSamples_nomt_BCells_snakemake/files/"



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




projTonsils <- loadArchRProject("All_Tonsils_unanalyzed")

# Only keep B cells
idx_b_cells <- projTonsils$Clusters_Harmony %in% c("C11", "C12", "C13",
                                                   "C2", "C3", "C4", "C5",
                                                   "C6", "C8", "C9", "C10")

cells_b_cells <- projTonsils$cellNames[idx_b_cells]

projTonsils2 <- subsetArchRProject(ArchRProj = projTonsils,
  cells = cells_b_cells, outputDirectory = output_dir)

getOutputDirectory(projTonsils2)

projTonsils2$Batch <- projTonsils2$Sample

projTonsils2$Batch[grepl("BCP", projTonsils2$Batch)] <- "Batch2"

projTonsils2$Batch[grepl("Tonsil", projTonsils2$Batch)] <- "Batch1"


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
    dimsToUse = 1:30,
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

# Batch correction by batch
projTonsils2 <- addHarmony(
    ArchRProj = projTonsils2,
    reducedDims = "IterativeLSI",
    name = "Harmony_batch",
    groupBy = "Batch",
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

projTonsils2 <- addClusters(
    input = projTonsils2,
    reducedDims = "Harmony_batch",
    method = "Seurat",
    name = "Clusters_Harmony_batch",
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
    minDist = 0.3, 
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
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

p3 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Sample", embedding = "UMAP_Harmony")

p4 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Batch", embedding = "UMAP_Harmony")

p5 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Clusters_Harmony", embedding = "UMAP_Harmony")

plotPDF(p3,p4,p5, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf",
  ArchRProj = projTonsils2, addDOC = FALSE, width = 5, height = 5)


# UMAP on batch correction
projTonsils2 <- addUMAP(
    ArchRProj = projTonsils2, 
    reducedDims = "Harmony_batch", 
    name = "UMAP_Harmony_batch", 
    nNeighbors = 30, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

p3 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Sample", embedding = "UMAP_Harmony_batch")

p4 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Batch", embedding = "UMAP_Harmony_batch")

p5 <- plotEmbedding(ArchRProj = projTonsils2, colorBy = "cellColData",
  name = "Clusters_Harmony", embedding = "UMAP_Harmony_batch")

plotPDF(p3,p4,p5, name = "Plot-UMAP2Harmony-Sample-Clusters-Batch.pdf",
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
seurat_object <- readRDS(paste0(seurat_dir, "BCells_dubRm_namedClust_seurat.rda"))
DefaultAssay(seurat_object) <- "RNA"
seurat_object[["SCT"]] <- NULL
seurat_object[["ADT"]] <- NULL




seurat_object$bcell_cluster_char <- "00"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "FCRL4+MCB"] <- "01"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "MCB"] <- "02"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Naive_B"] <- "03"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Activated_B"] <- "04"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Interferon_Active_B"] <- "05"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Circulating_Dark_GC_B"] <- "06"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Dark_GC_B"] <- "07"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Light_GC_B"] <- "08"
seurat_object$bcell_cluster_char[seurat_object$cell_type_bCell == "Plasmablasts"] <- "09"

seurat_object$num_cell_type_bCell <- paste0(seurat_object$bcell_cluster_char,
    "_", seurat_object$cell_type_bCell)

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
    groupRNA = "num_cell_type_bCell",
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
cActiveB <- paste0(paste0("0", c(1:3)), collapse = "|")

cNaiveB <- paste0(paste0("0", c(4:5)), collapse = "|")

cGCB <- paste0(paste0("0", c(6:8)), collapse = "|")

cPlasmablasts <- "09"

clustActiveB <- rownames(cM)[grep(cActiveB, preClust)]

clustNaiveB <- rownames(cM)[grep(cNaiveB, preClust)]

clustGCB <- rownames(cM)[grep(cGCB, preClust)]

clustPlasmablasts <- rownames(cM)[grep(cPlasmablasts, preClust)]

rnaActiveB <- colnames(seurat_object)[grep(cActiveB, seurat_object$num_cell_type_bCell)]

rnaNaiveB <- colnames(seurat_object)[grep(cNaiveB, seurat_object$num_cell_type_bCell)]

rnaGCB <- colnames(seurat_object)[grep(cGCB, seurat_object$num_cell_type_bCell)]

rnaPlasmablasts <- colnames(seurat_object)[grep(cPlasmablasts, seurat_object$num_cell_type_bCell)]


groupList <- SimpleList(
    ActiveB = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustActiveB],
        RNA = rnaActiveB
    ),
    NaiveB = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustNaiveB],
        RNA = rnaNaiveB
    ),
    GCB = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGCB],
        RNA = rnaGCB
    ),
    Plasmablasts = SimpleList(
        ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustPlasmablasts],
        RNA = rnaPlasmablasts)    
)

# Try 2
# cInterferon <- paste0(paste0("0", c(1, 5)), collapse = "|")

# cActiveB <- "04"

# cNaiveB <- paste0(paste0("0", c(2, 3)), collapse = "|")

# cGCB <- paste0(paste0("0", c(7:8)), collapse = "|")

# cCycling <- "06"

# cPlasmablasts <- "09"

# clustInterferon <- rownames(cM)[grep(cInterferon, preClust)]

# clustActiveB <- rownames(cM)[grep(cActiveB, preClust)]

# clustNaiveB <- rownames(cM)[grep(cNaiveB, preClust)]

# clustGCB <- rownames(cM)[grep(cGCB, preClust)]

# clustCycling <- rownames(cM)[grep(cCycling, preClust)]

# clustPlasmablasts <- rownames(cM)[grep(cPlasmablasts, preClust)]

# rnaInterferon <- colnames(seurat_object)[grep(cInterferon, seurat_object$num_cell_type_bCell)]

# rnaActiveB <- colnames(seurat_object)[grep(cActiveB, seurat_object$num_cell_type_bCell)]

# rnaNaiveB <- colnames(seurat_object)[grep(cNaiveB, seurat_object$num_cell_type_bCell)]

# rnaGCB <- colnames(seurat_object)[grep(cGCB, seurat_object$num_cell_type_bCell)]

# rnaCycling <- colnames(seurat_object)[grep(cCycling, seurat_object$num_cell_type_bCell)]

# rnaPlasmablasts <- colnames(seurat_object)[grep(cPlasmablasts, seurat_object$num_cell_type_bCell)]


# groupList <- SimpleList(
#     Interferon = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustInterferon],
#         RNA = rnaInterferon
#     ),
#     ActiveB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustActiveB],
#         RNA = rnaActiveB
#     ),
#     NaiveB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustNaiveB],
#         RNA = rnaNaiveB
#     ),
#     GCB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGCB],
#         RNA = rnaGCB
#     ),
#     Cycling = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustCycling],
#         RNA = rnaCycling
#     ),
#     Plasmablasts = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustPlasmablasts],
#         RNA = rnaPlasmablasts)    
# )

# Try 3
# cActiveB <- paste0(paste0("0", c(1:3)), collapse = "|")

# cNaiveB <- paste0(paste0("0", c(4:5)), collapse = "|")

# cGCB <- paste0(paste0("0", c(7:8)), collapse = "|")

# cCycling <- "06"

# cPlasmablasts <- "09"

# clustActiveB <- rownames(cM)[grep(cActiveB, preClust)]

# clustNaiveB <- rownames(cM)[grep(cNaiveB, preClust)]

# clustGCB <- rownames(cM)[grep(cGCB, preClust)]

# clustCycling <- rownames(cM)[grep(cCycling, preClust)]

# clustPlasmablasts <- rownames(cM)[grep(cPlasmablasts, preClust)]

# rnaActiveB <- colnames(seurat_object)[grep(cActiveB, seurat_object$num_cell_type_bCell)]

# rnaNaiveB <- colnames(seurat_object)[grep(cNaiveB, seurat_object$num_cell_type_bCell)]

# rnaGCB <- colnames(seurat_object)[grep(cGCB, seurat_object$num_cell_type_bCell)]

# rnaCycling <- colnames(seurat_object)[grep(cCycling, seurat_object$num_cell_type_bCell)]

# rnaPlasmablasts <- colnames(seurat_object)[grep(cPlasmablasts, seurat_object$num_cell_type_bCell)]


# groupList <- SimpleList(
#     ActiveB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustActiveB],
#         RNA = rnaActiveB
#     ),
#     NaiveB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustNaiveB],
#         RNA = rnaNaiveB
#     ),
#     GCB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGCB],
#         RNA = rnaGCB
#     ),
#     Cycling = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustCycling],
#         RNA = rnaCycling
#     ),
#     Plasmablasts = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustPlasmablasts],
#         RNA = rnaPlasmablasts)    
# )

# Try 4

# cNonGCB <- paste0(paste0("0", c(1:5)), collapse = "|")

# cGCB <- paste0(paste0("0", c(6:8)), collapse = "|")

# cPlasmablasts <- "09"

# clustNonGCB <- rownames(cM)[grep(cNonGCB, preClust)]

# clustGCB <- rownames(cM)[grep(cGCB, preClust)]

# clustPlasmablasts <- rownames(cM)[grep(cPlasmablasts, preClust)]

# rnaNonGCB <- colnames(seurat_object)[grep(cNonGCB, seurat_object$num_cell_type_bCell)]

# rnaGCB <- colnames(seurat_object)[grep(cGCB, seurat_object$num_cell_type_bCell)]

# rnaPlasmablasts <- colnames(seurat_object)[grep(cPlasmablasts, seurat_object$num_cell_type_bCell)]


# groupList <- SimpleList(
#     NonGCB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustNonGCB],
#         RNA = rnaNonGCB
#     ),
#     GCB = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustGCB],
#         RNA = rnaGCB
#     ),
#     Plasmablasts = SimpleList(
#         ATAC = projTonsils2$cellNames[projTonsils2$Clusters_Harmony %in% clustPlasmablasts],
#         RNA = rnaPlasmablasts)    
# )


projTonsils2 <- addGeneIntegrationMatrix(
    ArchRProj = projTonsils2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seurat_object,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "num_cell_type_bCell",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)


pal <- paletteDiscrete(values = seurat_object$num_cell_type_bCell)

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
    groupRNA = "num_cell_type_bCell",
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

remapClust <- c("01_FCRL4+MCB" = "FCRL4+MCB",
                "02_MCB" = "MCB",
                "03_Naive_B" = "Naive_B",
                "04_Activated_B" = "Activated_B",
                "05_Interferon_Active_B" = "Interferon_Active_B",
                "06_Circulating_Dark_GC_B" = "Circulating_Dark_GC_B",
                "07_Dark_GC_B" = "Dark_GC_B",
                "08_Light_GC_B" = "Light_GC_B",
                "09_Plasmablasts" = "Plasmablasts")

labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust),
  newLabels = remapClust)
labelNew2

projTonsils2$Clusters2 <- mapLabels(projTonsils2$Clusters_Harmony,
  newLabels = labelNew2, oldLabels = labelOld)


cell_types <- levels(seurat_object$cell_type_bCell)

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

# output_dir <- output_dir2

# saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
#   load = TRUE)

# Make pseudobulk replicates
projTonsils2 <- addGroupCoverages(ArchRProj = projTonsils2, groupBy = "Clusters2")

saveArchRProject(ArchRProj = projTonsils2, outputDirectory = output_dir,
  load = FALSE)

# Call peaks with Macs2
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

######## Get peak2gene links from heatmap

# ArchRProj <- projTonsils5
# groupBy <- "Clusters2"
# k <- 20
# corCutOff <- 0.45 
# FDRCutOff <- 0.0001
# nPlot <- 25000
# palGroup <- colors_cell_type_new
# limitsATAC <- c(-2, 2)
# limitsRNA <- c(-2, 2)
# groupBy <- "Clusters"
# palATAC <- paletteContinuous("solarExtra")
# palRNA <- paletteContinuous("blueYellow")


# ccd <- getCellColData(ArchRProj, select = groupBy)
# p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
# p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
# mATAC <- assay(readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ])
# mRNA <- assay(readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ])
# gc()
# KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
# KNNGroups <- lapply(seq_along(KNNList), function(x){
#     KNNx <- KNNList[[x]]
#     names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
# }) %>% unlist
# cD <- DataFrame(row.names=paste0("K", seq_len(ncol(mATAC))), groupBy = KNNGroups)
# pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
# if(!is.null(palGroup)){
#     pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
# }
# colorMap <- list(groupBy = pal)
# attr(colorMap[[1]], "discrete") <- TRUE

# rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
#   z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
#   if(limit){
#     z[z > max] <- max
#     z[z < min] <- min
#   }
#   return(z)
# }

# mATAC <- rowZscores(mATAC)
# mRNA <- rowZscores(mRNA)
# mATAC_names <- mATAC
# mRNA_names <- mRNA
# rownames(mATAC) <- NULL
# rownames(mRNA) <- NULL
# colnames(mATAC) <- paste0("K", seq_len(ncol(mATAC)))
# colnames(mRNA) <- paste0("K", seq_len(ncol(mRNA)))
# rownames(mATAC) <- paste0("P2G", seq_len(nrow(mATAC)))
# rownames(mRNA) <- paste0("P2G", seq_len(nrow(mRNA)))

# k1 <- kmeans(mATAC, k)
# if(nrow(mATAC) > nPlot){
#     nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
#     splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
#     kDF <- lapply(seq_along(splitK), function(x){
#       idx <- sample(splitK[[x]], floor(nPK[x]))
#       k <- rep(x, length(idx))
#       DataFrame(k = k, idx = idx)
#     }) %>% Reduce("rbind", .)
# }else{
#     kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
# }


# bS <- ArchR:::.binarySort(t(ArchR:::.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
# rowOrder <- rownames(bS[[1]])
# colOrder <- colnames(bS[[1]])
# kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))

# gene_set <- metadata(p2g)$geneSet

# peak_set <- metadata(p2g)$peakSet

# peak_set2 <- getPeakSet(projTonsils5)

# mATAC <- mATAC[kDF[,2],colOrder]
# mRNA <- mRNA[kDF[,2],colOrder]

# htATAC <- ArchR:::.ArchRHeatmap(
#     mat = mATAC[kDF[,2],colOrder],
#     scale = FALSE,
#     limits = limitsATAC,
#     color = palATAC, 
#     colData = cD[colOrder,,drop=FALSE],
#     colorMap = colorMap,
#     clusterCols = FALSE,
#     clusterRows = FALSE,
#     split = kDF[,3],
#     labelRows = FALSE,
#     labelCols = FALSE,
#     draw = FALSE,
#     name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
# )

# htRNA <- ArchR:::.ArchRHeatmap(
#     mat = mRNA[kDF[,2],colOrder], 
#     scale = FALSE,
#     limits = limitsRNA,
#     color = palRNA, 
#     colData = cD[colOrder,,drop=FALSE],
#     colorMap = colorMap,
#     clusterCols = FALSE,
#     clusterRows = FALSE,
#     split = kDF[,3],
#     labelRows = FALSE,
#     labelCols = FALSE,
#     draw = FALSE,
#     name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
# )


# output_dir <- getOutputDirectory(projTonsils5)
# pdf(paste0(output_dir, "/Plots/Plot-Peak2Gene_Heatmap-Kristen.pdf"), height = 10, width = 15)
# htATAC + htRNA
# dev.off()


# bS_names <- ArchR:::.binarySort(t(ArchR:::.groupMeans(t(mATAC_names[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
# rowOrder_names <- rownames(bS_names[[1]])
# colOrder_names <- colnames(bS_names[[1]])
# kDF_names <- kDF
# kDF_names[,3] <- as.integer(mapLabels(paste0(kDF_names[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))


# mATAC_names <- mATAC_names[kDF_names[,2],colOrder_names]
# mRNA_names <- mRNA_names[kDF_names[,2],colOrder_names]

# kDF_names$RNA_names <- rownames(mRNA_names)
# kDF_names$ATAC_names <- rownames(mATAC_names)

# kDF_names$RNA_names <- gsub("f", "", kDF_names$RNA_names)
# kDF_names$ATAC_names <- gsub("f", "", kDF_names$ATAC_names)

# indexes <- kDF_names$RNA_names
# gene_list <- lapply(indexes, function(x){
#     x <- as.integer(x)
#     gene_subset <- gene_set[x,]
#     gene_name <- gene_subset$name
#     return(gene_name)
#     })

# ATAC_indexes <- kDF_names$ATAC_names

# peak_list <- lapply(ATAC_indexes, function(x){
#     x <- as.integer(x)
#     peak_subset <- peak_set2[x,]
#     return(peak_subset)
#     })

# peak_list <- lapply(ATAC_indexes, function(x){
#     x <- as.integer(x)
#     peak_subset <- peak_set[x,]
#     return(peak_subset)
#     })


# gene_list <- unlist(gene_list)

# kDF_names$gene <- gene_list

# kDF_markers <- kDF_names[kDF_names$gene %in% markerGenes,]

# peak_df <- do.call(rbind, peak_list)

# gene_set[5,]

# peak_set[60,]