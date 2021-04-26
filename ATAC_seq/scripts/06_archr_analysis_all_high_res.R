library(ArchR)
set.seed(1)

addArchRGenome("hg38")

output_dir <- commandArgs(trailingOnly = TRUE)

#output_dir <- "All_Tonsils_high_res_new_chromvar"
#output_dir <- "All_Tonsils_new_chromvar_confident_peaks"
print(output_dir)
if(output_dir == "All_Tonsils_high_res_new_chromvar") {
    peak_reproducability <- "2"
    remove_gene <- TRUE
} else if (output_dir == "All_Tonsils_high_res_confident_peaks_new_chromvar"){
    peak_reproducability <- "(n+1)/2"
    remove_gene <- TRUE
} else if (output_dir == "All_Tonsils_high_res"){
    peak_reproducability <- "2"
    remove_gene <- FALSE
} else if (output_dir == "All_Tonsils_high_res_confident_peaks"){
    peak_reproducability <- "(n+1)/2"
    remove_gene <- FALSE
}

seurat_dir <- "/oak/stanford/groups/wjg/zshipony/EZH2_scRNA/ADT_190924_KLW_test/output/allSamples_nomt_snakemake/files/"

markerGenes <- c("CLEC4C", "CD14", "FCGR2A", "IFNGR1","GZMB", #Dendritic
    "LYZ", "VIM", "CST3", #monocytes
    "KLRG1", "GZMA", "IL32", "NKG7", "CCL5", #CD8T
    "IL6R", "CD28", "CD3D", #CD4T
    "PRDM1", "XBP1", "JCHAIN", "MZB1", #plasma B
    "LPP", "NEIL1", "FCRL3", "CD40", "MME", "HMMR", "LMO2", #Light GC
    "TUBA1B", "HMGB2", "STMN1", "TOP2A", #Dark GC
    "FCER2", "CD72", "TCL1A", "BANK1", # Naive B
    "FCRL4", "CCR6", "PLAC8", "TNFRSF13B", "CD83", # Activated B
    "CCR6", "FCRL4", "PTPN1", # FCRL4+MCB
    "ZBTB20", "CD27", "AHNAK", #MBC
    "CD72", "FAM129C", "KLF2", #Naive B
    "CD83", "CCND2", "FCER2", # Activated B
    "XAF1", "MX1", "IFI44L", #Interferon active B
    "PCNA", "TUBA1B", "TOP2A", # ciruclating dark GCB
    "AICDA", "MME", "SUGCT", # Dark zone B
    "NEIL1", "FGD6", "LMO2", #Light zone b
    "XBP1", "JCHAIN", "MZB1", #plasmablasts
    "PRDM1", "RORA", "RGS1", "CYTOR", "FOXP3", # T reg
    "CXCL13", # Tfh CXCL13
    "CXCR5", "FKBP5", "PDCD1", "BCL6", "PASK", "ST8SIA1", # Tfh
    "IL7R", "NABP1", "ANK3", "ZFP36L2", # Naive T cells
    "CCR7", "LEF1", "ITGA6", "KLF2", "SELL", #Tcm
    "GZMA", "CCL5", "GZMK", "CST7", #CTL
    "PRF1", "ANXA1", #CTL_PRF1
    "MKI67", "PCNA", "TOP2A", #Circulating
    "GNLY", "XCL1", "CTSW", "ID2" # NK
    )

markerGenes <- unique(markerGenes)

# load the project
projTonsils <- loadArchRProject("All_tonsils_high_res_unanalyzed")

projTonsils <- saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = TRUE)

projTonsils <- addImputeWeights(projTonsils)

p2 <- plotEmbedding(
    ArchRProj = projTonsils, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils)
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-New-Marker-Genes-W-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils, 
    addDOC = FALSE, width = 5, height = 5)


p2 <- plotEmbedding(
    ArchRProj = projTonsils, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils)
)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-New-Marker-Genes-W-Imputation-Harmony.pdf", 
    ArchRProj = projTonsils, 
    addDOC = FALSE, width = 5, height = 5)


cell_types <- unique(projTonsils$RNA_high_res)

cell_types <- cell_types[order(cell_types)]

colors_cell_type_new <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

names(colors_cell_type_new) <- cell_types

p1 <- plotEmbedding(projTonsils, colorBy = "cellColData", name = "RNA_high_res",
    pal = colors_cell_type_new, embedding = "UMAP_Harmony")

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projTonsils,
  addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

projTonsils <- loadArchRProject(output_dir)

saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

# Make pseudobulk replicates
projTonsils <- addGroupCoverages(ArchRProj = projTonsils, groupBy = "RNA_high_res")


# Call peaks with Macs2
pathToMacs2 <- findMacs2()

projTonsils <- addReproduciblePeakSet(
    ArchRProj = projTonsils, 
    groupBy = "RNA_high_res", 
    pathToMacs2 = pathToMacs2,
    reproducability = peak_reproducability
)
getPeakSet(projTonsils)

peak_set <- getPeakSet(projTonsils)
outDir <- getOutputDirectory(projTonsils)

saveRDS(peak_set, paste0(outDir, "/files/macs2_peaks.rds"))

saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

# Add peak matrix
projTonsils <- addPeakMatrix(projTonsils)

saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

outDir <- getOutputDirectory(projTonsils)


getAvailableMatrices(projTonsils)

table(projTonsils$RNA_high_res)

# Find marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = projTonsils, 
    useMatrix = "PeakMatrix", 
    groupBy = "RNA_high_res",
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
  ArchRProj = projTonsils, addDOC = FALSE)


#############################
#        Chromvar           #
#############################

projTonsils <- addMotifAnnotations(ArchRProj = projTonsils,
  motifSet = "cisbp", name = "Motif")

# Motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projTonsils,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Change n based on what is good
# Can also update this to be the values called by the correlation analysis
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 30, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 12,
  height = 6, ArchRProj = projTonsils, addDOC = FALSE)

# Motif deviations

# Check if motif annotations added
if("Motif" %ni% names(projTonsils@peakAnnotation)){
    projTonsils <- addMotifAnnotations(ArchRProj = projTonsils,
      motifSet = "cisbp", name = "Motif")
}

projTonsils <- addBgdPeaks(projTonsils)


anno <- getPeakAnnotation(projTonsils, "Motif")
matches <- readRDS(anno$Matches)

if(remove_gene){
    # This gene was previously returning "NA" values, so I'm removing it
    matches <- matches[,colnames(matches) != "ENSG00000250542_156"]
}

projTonsils <- addDeviationsMatrix(
  ArchRProj = projTonsils, 
  peakAnnotation = "Motif",
  matches = matches
)

saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

projTonsils <- loadArchRProject(output_dir)

#####################
# I saved here and quit. Saving somehow removes the peak set
outDir <- getOutputDirectory(projTonsils)
peakSet <- readRDS(paste0(outDir, "/files/macs2_peaks.rds"))

projTonsils <- addPeakSet(ArchRProj = projTonsils, peakSet = peakSet,
  force = TRUE)

getPeakSet(projTonsils)

plotVarDev <- getVarDeviations(projTonsils, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5,
  height = 5, ArchRProj = projTonsils, addDOC = FALSE)

motifs <- markerGenes
markerMotifs <- getFeatures(projTonsils, select = paste(motifs, collapse="|"),
  useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

projTonsils <- addImputeWeights(projTonsils)

p <- plotGroups(ArchRProj = projTonsils, 
  groupBy = "RNA_high_res", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projTonsils)
)

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5,
  height = 5, ArchRProj = projTonsils, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projTonsils, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projTonsils)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-Deviations-W-Imputation.pdf", 
    ArchRProj = projTonsils, 
    addDOC = FALSE, width = 5, height = 5)

# Motif footprinting
motifPositions <- getPositions(projTonsils)
motifPositions

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
  value = TRUE)))
markerMotifs

# Can make this for only a subset of cell types using the useGroups argument
seFoot <- getFootprints(
  ArchRProj = projTonsils, 
  positions = motifPositions[markerMotifs], 
  groupBy = "RNA_high_res"
)

# Subtract Tn5 bias
p1 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Divide by Tn5 bias
p2 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# No normalization
p3 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = projTonsils, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5,
  pal = colors_cell_type_new
)

# Co-accessibility
projTonsils <- addCoAccessibility(
    ArchRProj = projTonsils,
    reducedDims = "Harmony"
)

cA <- getCoAccessibility(
    ArchRProj = projTonsils,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

p <- plotBrowserTrack(
    ArchRProj = projTonsils, 
    groupBy = "RNA_high_res", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projTonsils)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projTonsils, 
    addDOC = FALSE, width = 5, height = 5)

# Peak to gene links
projTonsils <- addPeak2GeneLinks(
    ArchRProj = projTonsils,
    reducedDims = "Harmony"
)


p <- plotBrowserTrack(
    ArchRProj = projTonsils, 
    groupBy = "RNA_high_res", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projTonsils)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = projTonsils, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotPeak2GeneHeatmap(ArchRProj = projTonsils, groupBy = "RNA_high_res", k = 20,
    palGroup = colors_cell_type_new)
output_dir <- getOutputDirectory(projTonsils)
pdf(paste0(output_dir, "/Plots/Plot-Peak2Gene_Heatmap.pdf"), height = 10, width = 15)
p
dev.off()

# Find positive TF regulators
seGroupMotif <- getGroupSE(ArchRProj = projTonsils,
  useMatrix = "MotifMatrix", groupBy = "RNA_high_res")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = projTonsils,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)

corGSM_MM

corGIM_MM <- correlateMatrices(
    ArchRProj = projTonsils,
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



saveArchRProject(ArchRProj = projTonsils, outputDirectory = output_dir,
  load = FALSE)

