library(ArchR)
set.seed(1)


addArchRGenome("hg38")

function_path <- "/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/scripts/functions.R"
gwas_snps <- read.csv("/oak/stanford/groups/wjg/zshipony/EZH2_scATAC/GWAS_SNPs_hg38.csv")

ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))


#ATAC_dir <- commandArgs(trailingOnly = TRUE)
ATAC_dir <- "All_Tonsils_high_res_new_chromvar"
#cluster <- "Clusters2"
cluster <- "RNA_high_res"

if(ATAC_dir == "All_Tonsils_new_chromvar"){
	seed <- 4
} else {
	seed <- 1
}


markerGenes <- c(
    "TOP2A", "PCNA", "MKI67", "ICAM1", "SEMA7A", "BACH2", #Germinal center dz
    "MS4A1", "STX7", "HMMR", "PLK1", "POU2F2", "CXCR4", # germinal center lz
    "XBP1", "PRDM1", # plasmablasts
    "CD27", "AIM2", "SCIMP", # MBC
    "CXCR6", "FCRL4", "NEAT1", # FCRL4 MBC
    "CD72", "TCL1A", # Naive B cells
    "NFKB1", "CCND2", "IRF4", "CD83", "EGR3", "TCF4", # preGC/Acrivated B
    "CCR7", "CXCL13", "FKBP5", "LEF1", "KLF2", "KLF6", "RORA", "IL7R", "NKG7",
    "ANXA1", "TXNIP", "GZMA", "CLDND1", "BANK1", "ANXA2", "XCL1",
    "CD81", "CD8A", "CD4", "CD3D", "CD3G", # T cell
    "CLEC4C", "CD14", "FCGR2A", "IFNGR1", # DC
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
    "GZMA", "CCL5", "GZMK", "CST7", "CCL4", #CTL
    "PRF1", "ANXA1", #CTL_PRF1
    "MKI67", "PCNA", "TOP2A", #Circulating
    "GNLY", "XCL1", "CTSW", "ID2" # NK
)

markerGenes <- unique(markerGenes)
if(cluster == "RNA_high_res"){
    cell_types <- c("MCB", "FCRL4+MCB", "Naive_B", "Activated_B",
    "Circulating_Dark_GC_B", "Light_GC_B", "Plasmablasts",
    "Circulating", "Treg", "Naive_TCells", "Tfh", "Tcm_CD4",
    "Tcm_CD8",  "CTL", "Monocyte_DC", "Dendritic")

    colors_cell_type_new <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

    names(colors_cell_type_new) <- cell_types

    colors_cell_type_new <- sort(colors_cell_type_new)
    
} else if (cluster == "Clusters2"){

    cell_types <- c("Naive_B", "Activated_B", "Light_GC_B", "Dark_GC_B",
                    "Plasma_B", "CD4_T", "CD8_T", "Monocyte_DC", "Dendritic")

    colors_cell_type_new <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(length(cell_types))

    names(colors_cell_type_new) <- cell_types

    colors_cell_type_new <- colors_cell_type_new[1:8]

    colors_cell_type_new <- sort(colors_cell_type_new)
}

source(function_path)

projTonsils5 <- loadArchRProject(ATAC_dir)

outDir <- getOutputDirectory(projTonsils5)
output_dir <- getOutputDirectory(projTonsils5)


peak_set <- readRDS(paste0(output_dir, "/files/macs2_peaks.rds"))

projTonsils5 <- addPeakSet(projTonsils5, peakSet = peak_set, force = TRUE)

######## Get peak2gene links from heatmap
projTonsils5 <- addPeak2GeneLinks(
  ArchRProj = projTonsils5,
  reducedDims = "Harmony"
)


p2gLinks <- getPeak2GeneLinks(projTonsils5, returnLoops = FALSE)

# Check if indexes of peak2gene links align with the peak set
identical(metadata(p2gLinks)[[1]]$ranges, peak_set$ranges)


# get matrix of all peak2gene links
ArchRProj <- projTonsils5
k <- 20
corCutOff <- 0.45 
FDRCutOff <- 0.0001
#nPlot <- 25000
palGroup <- colors_cell_type_new
limitsATAC <- c(-2, 2)
limitsRNA <- c(-2, 2)
groupBy <- cluster
palATAC <- paletteContinuous("solarExtra")
palRNA <- paletteContinuous("blueYellow")


ccd <- getCellColData(ArchRProj, select = groupBy)
p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
mATAC <- assay(readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ])
mRNA <- assay(readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ])

# Update nPlot to be all peak 2 gene links
nPlot <- nrow(mATAC)

gc()
KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
KNNGroups <- lapply(seq_along(KNNList), function(x){
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
}) %>% unlist
cD <- DataFrame(row.names=paste0("K", seq_len(ncol(mATAC))), groupBy = KNNGroups)
pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
if(!is.null(palGroup)){
    pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
}
colorMap <- list(groupBy = pal)
attr(colorMap[[1]], "discrete") <- TRUE

mATAC <- ArchR:::.rowZscores(mATAC)
mRNA <- ArchR:::.rowZscores(mRNA)
mATAC_names <- mATAC
mRNA_names <- mRNA
rownames(mATAC) <- NULL
rownames(mRNA) <- NULL
colnames(mATAC) <- paste0("K", seq_len(ncol(mATAC)))
colnames(mRNA) <- paste0("K", seq_len(ncol(mRNA)))
rownames(mATAC) <- paste0("P2G", seq_len(nrow(mATAC)))
rownames(mRNA) <- paste0("P2G", seq_len(nrow(mRNA)))

mATAC_names_df <- data.frame(ATAC_names = rownames(mATAC_names),
    ATAC_new_names = rownames(mATAC))

mATAC_names_df$ATAC_new_names <- as.character(mATAC_names_df$ATAC_new_names)

mRNA_names_df <- data.frame(RNA_names = rownames(mRNA_names),
    RNA_new_names = rownames(mRNA))

mRNA_names_df$RNA_new_names <- as.character(mRNA_names_df$RNA_new_names)

set.seed(seed)
k1 <- kmeans(mATAC, k)
if(nrow(mATAC) > nPlot){
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
}else{
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
}

set.seed(seed)
subset_plot <- 25000
nPK <- subset_plot * table(k1$cluster) / length(k1$cluster) 
splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
kDF_subset <- lapply(seq_along(splitK), function(x){
  idx <- sample(splitK[[x]], floor(nPK[x]))
  k <- rep(x, length(idx))
  DataFrame(k = k, idx = idx)
}) %>% Reduce("rbind", .)

bS <- ArchR:::.binarySort(t(ArchR:::.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
rowOrder <- rownames(bS[[1]])
colOrder <- colnames(bS[[1]])
kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))

if(!identical(mATAC_names_df$ATAC_new_names, rownames(kDF))){
    mATAC_names_df <- mATAC_names_df[match(rownames(kDF),
        mATAC_names_df$ATAC_new_names)]
}
kDF$ATAC_idx <- mATAC_names_df$ATAC_names

if(!identical(mRNA_names_df$RNA_new_names, rownames(kDF))){
    mRNA_names_df <- mRNA_names_df[match(rownames(kDF),
        mRNA_names_df$RNA_new_names)]
}
kDF$RNA_idx <- mRNA_names_df$RNA_names


gene_set <- metadata(p2g)$geneSet

peak_set <- metadata(p2g)$peakSet

saveRDS(peak_set, paste0(outDir, "/files/p2g_peak_set.rda"))

kDF$RNA_idx <- gsub("f", "", kDF$RNA_idx)
kDF$ATAC_idx <- gsub("f", "", kDF$ATAC_idx)

indexes <- kDF$RNA_idx
gene_list <- lapply(indexes, function(x){
    x <- as.integer(x)
    gene_subset <- gene_set[x,]
    gene_name <- gene_subset$name
    return(gene_name)
    })


gene_list <- unlist(gene_list)

kDF$gene <- gene_list

write.csv(kDF, paste0(outDir, "/files/peak2gene_links_", cluster, "_merged.csv"))
saveRDS(colOrder, paste0(outDir, "/files/peak2gene_links_", cluster , 
    "col_order.rda"))
write.csv(kDF_subset, paste0(outDir, "/files/peak2gene_subset_peaks.csv"))
saveRDS(mATAC, paste0(outDir, "/files/ATAC_single_cell_mat_for_gwas.rda"))
saveRDS(mRNA, paste0(outDir, "/files/RNA_single_cell_mat_for_gwas.rda"))

peak_2_gene <- read.csv(paste0(outDir, "/files/peak2gene_links_", cluster, "_merged.csv"))
colOrder <- readRDS(paste0(outDir, "/files/peak2gene_links_", cluster , 
    "col_order.rda"))
kDF_subset <- read.csv(paste0(outDir, "/files/peak2gene_subset_peaks.csv"))


##################
# Gwas
# Get gwas info
# gwas_snps$start <- gwas_snps$pos
# gwas_snps$end <- gwas_snps$pos

gwas_snps$start <- gwas_snps$hg38_pos
gwas_snps$end <- gwas_snps$hg38_pos

granges_gwas <- makeGRangesFromDataFrame(gwas_snps, keep.extra.columns = TRUE)

# Find overlaps between the snps and peaks
overlap <- GenomicRanges::findOverlaps(query = granges_gwas,
    subject = peak_set)


overlap <- as.data.frame(overlap)

overlap$ATAC_idx <- overlap$subjectHits

# merge with DF of peak to gene links
# Merging is better because there may be several snps in the same peak
# so we want to duplicate that
overlap_peaks <- merge(overlap, peak_2_gene, by = "ATAC_idx",
    no.dups = TRUE)

gwas_snps$gwas_idx <- 1:nrow(gwas_snps)

# subset the gwas DF
# same here, we want to allow for duplicates
indexes <- overlap_peaks$queryHits
gene_list <- lapply(indexes, function(x){
    x <- as.integer(x)
    gwas_subset <- gwas_snps[x,]
    return(gwas_subset)
    })

gwas_df <- do.call(rbind, gene_list)

gwas_df <- gwas_df[,c(1:11, 47, 50)]


# Combine with the peak to gene links
overlap_peaks <- cbind(overlap_peaks, gwas_df)

# Order by disease
disease <- c(1:length(unique(overlap_peaks$Disease)))
names(disease) <- as.character(unique(overlap_peaks$Disease))

kdf_disease <- as.character(overlap_peaks$Disease)

disease_num <- disease[kdf_disease]

overlap_peaks$disease_num <- disease_num

print(disease)

write.csv(disease, paste0(outDir, "/files/gwas_peak2gene_table.csv"))

write.csv(overlap_peaks, paste0(outDir, "/files/gwas_overlap_peaks.csv"))

library(tidyverse)

# Remove any rows where the RNA index, ATAC index, and disease are the same
overlap_peaks_unique <- overlap_peaks[!duplicated(overlap_peaks[, c("ATAC_idx",
    "RNA_idx", "Disease")]), ]


write.csv(overlap_peaks_unique,
    paste0(outDir, "/files/peak2gene_links_gwas_hits.csv"))


overlap_peaks_unique_2 <- overlap_peaks_unique[!duplicated(
    overlap_peaks_unique[, c("RNA_idx", "Disease")]), ]

# Determine if there is an enrichment from the total gwas for immune related
total_disease_snps <- data.frame(table(gwas_snps$Disease))

overlap_disease_snps <- data.frame(table(overlap_peaks$Disease))

total_disease_snps$type <- "all_snps"

overlap_disease_snps$type <- "overlapping_snps"

total_disease_snps$percent <- total_disease_snps$Freq/sum(total_disease_snps$Freq)

overlap_disease_snps$percent <- overlap_disease_snps$Freq/sum(overlap_disease_snps$Freq)

disease_snps <- rbind(total_disease_snps, overlap_disease_snps)

colnames(disease_snps) <- c("Disease", "Frequency", "Percent", "Type")

autoimmune_disease <- c("Asthma", "Atopic_dermatitis", "Allergy", 
    "Primary_sclerosing_cholangitis", "Alopecia_areata", 
    "Juvenile_idiopathic_arthritis", "Systemic_sclerosis",
    "Vitiligo", "Type_1_diabetes", "Autoimmune_thyroiditis",
    "Rheumatoid_arthritis", "Multiple_sclerosis", "Celiac_disease",
    "Primary_biliary_cirrhosis", "Systemic_lupus_erythematosus",
    "Kawasaki_disease", "Behcets_disease", "Psoriasis",
    "Ankylosing_spondylitis", "Crohns_disease", "Ulcerative_colitis")

autoimmune_bool <- levels(total_disease_snps$Var1) %in% autoimmune_disease

text_color <- rep("black", length(autoimmune_bool))

text_color[autoimmune_bool] <- "red"

autoimmune_snps <- disease_snps[disease_snps$Disease %in% autoimmune_disease, ]

# Make a bar plot of this
pdf(paste0(output_dir, "/Plots/GWAS_disease_percent.pdf"),
    height = 8, width = 11.5)
ggplot2::ggplot(disease_snps, aes(y = Percent, x = Disease, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
        color = text_color))

ggplot2::ggplot(autoimmune_snps, aes(y = Percent, x = Disease, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
dev.off()


# Only plot autoimmune diseases

# Remove any rows where the RNA index, ATAC index, and disease are the same
overlap_peaks_autoimmune <- overlap_peaks_unique[overlap_peaks_unique$Disease %in%
    autoimmune_disease, ]


write.csv(overlap_peaks_autoimmune,
    paste0(outDir, "/files/peak2gene_links_gwas_hits_autoimmune.csv"))


