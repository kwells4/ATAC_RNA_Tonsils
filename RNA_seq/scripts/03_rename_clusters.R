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
                 "Activated_B", #2
                 "Naive_B", #3
                 "Light_GC_B", #4
                 "CD4_T", #5
                 "Dark_GC_B", #6
                 "CD4_T",#7
                 "Plasma_B", #8
                 "CD8_T", #9
                 "Dark_GC_B", #10
                 "CD4_T", #11
                 "CD4_T", #12
                 "Naive_B", #13
                 "Dark_GC_B", #14
                 "Naive_B", #15
                 "Monocyte_DC", #16
                 "Dark_GC_B", #17
                 "Light_GC_B", #18
                 "CD4_T", #19
                 "Dendritic",  #20
                 "Naive_B", #21
                 "Plasma_B" #22
                 )

# Will likely need to change this line to make a character
names(cell_names) <- 0:22

idents <- seurat_object$seurat_clusters

cell_labels_char <- c("01", "05", "06", "04", "07", "17", "10", "19",
    "13", "20", "12", "18", "16", "00", "09", "03", "22", "11", "08",
    "15", "21", "02", "14")

names(cell_labels_char) <- 0:22

all_labels_char <- cell_labels_char[idents]

all_labels_char <- factor(all_labels_char)

seurat_object$new_cluster_char <- all_labels_char

cell_labels <- c(1, 5, 6, 4, 7, 17, 10, 19, 13, 20, 12, 18, 16, 0, 9, 3,
                 22, 11, 8, 15, 21, 2, 14)
names(cell_labels) <- 0:22

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

