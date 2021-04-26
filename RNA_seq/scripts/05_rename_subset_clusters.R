library(Seurat)

seurat_path <- snakemake@input[[1]]
function_path <- snakemake@params[[1]]
save_dir <- snakemake@params[[2]]
subset <- snakemake@params[[3]]
save_seurat <- snakemake@output[[1]]

source(function_path)


print(save_dir)
seurat_object <- readRDS(seurat_path)
ggplot2::theme_set(ggplot2::theme_classic())

if(subset == "GCBcells"){

	cell_labels <- c(2, 1, 13, 3, 9, 6, 8, 4, 7, 5, 12, 10, 0, 2)

	names(cell_labels) <- 0:13

	idents <- Idents(seurat_object)

	all_labels <- cell_labels[idents]

	all_labels <- factor(all_labels)

	# add to metadata 
	seurat_object$new_cluster <- all_labels

	print(head(seurat_object[[]]))

	saveRDS(seurat_object, file = paste0(save_seurat))


	plotDimRed(seurat_object, col_by = c("cell_type", "new_cluster"), 
	    save_plot = paste0(save_dir, "images/new_clusters.pdf"), size = 0.5)

} else if (subset == "BCells"){
	cell_labels <- c(2, 1, 0, 3, 7, 8, 5, 5, 5, 6, 5, 2, 8, 2, 7, 4, 8)

	names(cell_labels) <- 0:16

	idents <- Idents(seurat_object)

	all_labels <- cell_labels[idents]

	all_labels <- factor(all_labels)

	# add to metadata 
	seurat_object$new_cluster <- all_labels

	cell_type_labels <- c("Naive_B", #2
		                  "MCB", #1
		                  "FCRL4+MCB", #0
		                  "Activated_B", # 3
		                  "Light_GC_B", #7
		                  "Plasmablasts", #8
		                  "Circulating_Dark_GC_B", #5
		                  "Circulating_Dark_GC_B", #5
		                  "Circulating_Dark_GC_B", #5
		                  "Dark_GC_B", #6
		                  "Circulating_Dark_GC_B", #5
		                  "Naive_B", #2
		                  "Plasmablasts", #8
		                  "Naive_B", #2
		                  "Light_GC_B", #7
		                  "Interferon_Active_B", #4
		                  "Plasmablasts" #8
		                  )

	names(cell_type_labels) <- 0:16

	all_labels_cell_type <- cell_type_labels[idents]

	all_labels_cell_type <- factor(all_labels_cell_type,
		levels = c("FCRL4+MCB", "MCB", "Naive_B", "Activated_B",
			       "Interferon_Active_B", "Circulating_Dark_GC_B",
			       "Dark_GC_B", "Light_GC_B", "Plasmablasts"))

	# add to metadata 
	seurat_object$cell_type_bCell <- all_labels_cell_type


	print(head(seurat_object[[]]))

	saveRDS(seurat_object, file = paste0(save_seurat))


	plotDimRed(seurat_object, col_by = c("cell_type", "new_cluster", "cell_type_bCell"), 
	    save_plot = paste0(save_dir, "images/new_clusters.pdf"), size = 0.5)

} else if (subset == "TCells"){
	cell_labels <- c(1, 3, 4, 0, 6, 2, 6, 5, 8, 9, 7)

	names(cell_labels) <- 0:10

	idents <- Idents(seurat_object)

	all_labels <- cell_labels[idents]

	all_labels <- factor(all_labels)

	# add to metadata 
	seurat_object$new_cluster <- all_labels

	cell_type_labels <- c("Tfh_CXCL13", #2
		                  "Naive_TCells", #3
		                  "Tcm_CD4", #4
		                  "Treg", # 0
		                  "CTL", #6
		                  "Tfh", #2
		                  "CTL", #6
		                  "Tcm_CD8", #5
		                  "Circulating", #8
		                  "NK_cells", #9
		                  "CTL_PRF1" #7
		                  )

	names(cell_type_labels) <- 0:10

	all_labels_cell_type <- cell_type_labels[idents]

	all_labels_cell_type <- factor(all_labels_cell_type,
		levels = c("Treg", "Tfh", "Tfh_CXCL13", "Naive_TCells",
			       "Tcm_CD4", "Tcm_CD8",
			       "CTL", "CTL_PRF1", "Circulating", "NK_cells"))

	# add to metadata 
	seurat_object$cell_type_tCell <- all_labels_cell_type


	print(head(seurat_object[[]]))

	saveRDS(seurat_object, file = paste0(save_seurat))


	plotDimRed(seurat_object, col_by = c("cell_type", "new_cluster", "cell_type_tCell"), 
	    save_plot = paste0(save_dir, "images/new_clusters.pdf"), size = 0.5)

} else if (grepl("_no_cycling", subset)){
	saveRDS(seurat_object, file = paste0(save_seurat))
	
} else {
	print(paste0("NO RENAMING PROVIDED.
		PLEASE UPDATE 'scripts/final_scripts/05_rename_subset_clusters.R'
		to include ", subset, " subset"))
}