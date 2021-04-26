library(Seurat)

seurat_path <- snakemake@input[[1]]
function_path <- snakemake@params[[1]]
save_dir <- snakemake@params[[2]]
save_seurat <- snakemake@output[[1]]

source(function_path)


print(save_dir)
seurat_object <- readRDS(seurat_path)

cell_names <- c("Naive_B", #0
                 "Activated_B", #1
                 "Light_GC_B", #2
                 "Activated_B", #3
                 "Naive_B", #4
                 "CD4_T", #5
                 "CD4_T", #6
                 "Plasma_B",#7
                 "CD4_T", #8
                 "CD8_T", #9
                 "Dark_GC_B", #10
                 "Naive_B", #11
                 "Dark_GC_B", #12
                 "CD4_T", #13
                 "Naive_B", #14
                 "Dark_GC_B", #15
                 "Naive_B", #16
                 "Monocyte_DC", #17
                 "CD4_T", #18
                 "Dendritic", #19
                 "Naive_B"  #20
                 )

# Will likely need to change this line to make a character
names(cell_names) <- 0:20

idents <- seurat_object$seurat_clusters

cell_labels_char <- c("01", "07", "08", "06", "04", "16", "13", "12",
    "14", "18", "11", "03", "10", "15", "02", "09", "00", "20", "17",
    "19", "05")

names(cell_labels_char) <- 0:20

all_labels_char <- cell_labels_char[idents]

all_labels_char <- factor(all_labels_char)

seurat_object$new_cluster_char <- all_labels_char

cell_labels <- c(1, 7, 8, 6, 4, 16, 13, 12, 14, 18, 11, 3, 10, 15, 2, 9, 0,
  20, 17, 19, 5)

names(cell_labels) <- 0:20

all_labels <- cell_labels[idents]

all_labels <- factor(all_labels)

# add to metadata 
seurat_object$new_cluster <- all_labels

idents <- Idents(seurat_object)

all_labels <- cell_names[idents]

all_labels <- factor(all_labels)

seurat_object$cell_type <- all_labels

seurat_object$num_cell_type <- paste0(seurat_object$new_cluster_char, "_",
                                      seurat_object$cell_type)

print(head(seurat_object[[]]))

saveRDS(seurat_object, file = paste0(save_seurat))


plotDimRed(seurat_object, col_by = c("cell_type", "num_cell_type", "cluster"), 
    save_plot = paste0(save_dir, "images/named_clusters.pdf"))

