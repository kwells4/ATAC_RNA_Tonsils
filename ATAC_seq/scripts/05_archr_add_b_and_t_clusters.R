library(ArchR)
set.seed(1)

addArchRGenome("hg38")


# all_data <- "All_Tonsils_pre_peaks"

# output_directory <- "All_tonsils_high_res_unanalyzed"

all_data <- "All_Tonsils_new_chromvar"

output_directory <- "All_Tonsils_new_chromvar"

b_cell <- "B_cells"

t_cell <- "T_cells"

# Remove the doublet cluster
proj_all <- loadArchRProject(all_data)

proj_bcells <- loadArchRProject(b_cell)

proj_tcells <- loadArchRProject(t_cell)

# Transfer labels from b cells to all cells
all_meta_data <- getCellColData(ArchRProj = proj_all)

Bcell_meta_data <- getCellColData(ArchRProj = proj_bcells, select = "Clusters2")

Tcell_meta_data <- getCellColData(ArchRProj = proj_tcells, select = "Clusters2")

all_meta_data$in_others <- FALSE
all_meta_data[intersect(rownames(all_meta_data),
	rownames(Bcell_meta_data)), ]$in_others <- TRUE

all_meta_data[intersect(rownames(all_meta_data),
	rownames(Tcell_meta_data)), ]$in_others <- TRUE

all_meta_data[all_meta_data$Clusters2 %in% 
	c("Dendritic", "Monocyte_DC"), ]$in_others <- TRUE

colnames(Bcell_meta_data) <- "Clusters2_Bcell"

colnames(Tcell_meta_data) <- "Clusters2_Tcell"

all_meta_data <- merge(all_meta_data, Bcell_meta_data, by = "row.names",
	all.x = TRUE)

rownames(all_meta_data) <- all_meta_data$Row.names
all_meta_data$Row.names <- NULL

all_meta_data$RNA_high_res <- "other"

all_meta_data$RNA_high_res[!is.na(all_meta_data$Clusters2_Bcell)] <-
	all_meta_data$Clusters2_Bcell[!is.na(all_meta_data$Clusters2_Bcell)]

# Add T cells
all_meta_data <- merge(all_meta_data, Tcell_meta_data, by = "row.names",
	all.x = TRUE)

# Same problem here, it doesn't perfectly map
all_meta_data$RNA_high_res[!is.na(all_meta_data$Clusters2_Tcell)] <-
	all_meta_data$Clusters2_Tcell[!is.na(all_meta_data$Clusters2_Tcell)]

rownames(all_meta_data) <- all_meta_data$Row.names
all_meta_data$Row.names <- NULL

all_meta_data$RNA_high_res[all_meta_data$Clusters2 %in% 
	c("Dendritic", "Monocyte_DC")] <- all_meta_data$Clusters2[all_meta_data$Clusters2 %in% 
	c("Dendritic", "Monocyte_DC")]

cM <- as.matrix(confusionMatrix(all_meta_data$Clusters_Harmony,
  all_meta_data$RNA_high_res))

preClust <- colnames(cM)[apply(cM, 1 , which.max)]
preClust <- data.frame(cbind(preClust, rownames(cM))) #Assignments

for (i in grep("other", all_meta_data$RNA_high_res)){
	cluster <- all_meta_data[i,]$Clusters_Harmony
	new_cluster <- as.character(preClust[preClust$V2 == cluster, ]$preClust)
	all_meta_data$RNA_high_res[i] <- new_cluster
}

all_meta_data <- all_meta_data[match(rownames(proj_all), rownames(all_meta_data)), ]

proj_all$RNA_high_res <- all_meta_data$RNA_high_res

saveArchRProject(ArchRProj = proj_all, outputDirectory = output_directory)