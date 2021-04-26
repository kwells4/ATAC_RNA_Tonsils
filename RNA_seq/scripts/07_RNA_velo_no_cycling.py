from PIL import Image
import scvelo as scv
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt

scv.settings.set_figure_params('scvelo')


# # Germinal center cells
# input_loom = "output/allSamples/files/allSamples_RNA_velo.loom"
# cluster_png = "output/allSamples_GCBcells_snakemake/images/rna_velocity_GCBcells_cluster_harmony.png"
# cell_type_png = "output/allSamples_GCBcells_snakemake/images/rna_velocity_GCBcells_cell_type_harmony.png"
# seurat_cells = "output/allSamples_GCBcells_snakemake/files/GCBcells_cell_info.txt"
# seurat_clusters = "output/allSamples_GCBcells_snakemake/files/GCBcells_cluster_info.txt"
# seurat_umap = "output/allSamples_GCBcells_snakemake/files/GCBcells_UMAP_info.txt"
# seurat_harmony = "output/allSamples_GCBcells_snakemake/files/GCBcells_harmony_info.txt"
# seurat_pca = "output/allSamples_GCBcells_snakemake/files/GCBcells_pca_info.txt"
# subgroup = "GCBcells"

# input_loom = "output/allSamples/files/allSamples_RNA_velo.loom"
# cluster_png = "output/allSamples_nomt_GCBcells_snakemake/images/rna_velocity_GCBcells_cluster_harmony.png"
# cell_type_png = "output/allSamples_nomt_GCBcells_snakemake/images/rna_velocity_GCBcells_cell_type_harmony.png"
# seurat_cells = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_cell_info.txt"
# seurat_clusters = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_cluster_info.txt"
# seurat_umap = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_UMAP_info.txt"
# seurat_harmony = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_harmony_info.txt"
# seurat_pca = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_pca_info.txt"
# subgroup = "GCBcells"
# save_dir = "output/allSamples_nomt_GCBcells_snakemake/"
# pseudotime_output = "output/allSamples_nomt_GCBcells_snakemake/files/GCBcells_pseudotime.txt"

# All B Cells
# input_loom = "output/allSamples/files/allSamples_RNA_velo.loom"
# cluster_png = "output/allSamples_nomt_BCells_snakemake/images/rna_velocity_BCells_cluster_harmony.png"
# cell_type_png = "output/allSamples_nomt_BCells_snakemake/images/rna_velocity_BCells_cell_type_harmony.png"
# seurat_cells = "output/allSamples_nomt_BCells_snakemake/files/BCells_cell_info.txt"
# seurat_clusters = "output/allSamples_nomt_BCells_snakemake/files/BCells_cluster_info.txt"
# seurat_umap = "output/allSamples_nomt_BCells_snakemake/files/BCells_UMAP_info.txt"
# seurat_harmony = "output/allSamples_nomt_BCells_snakemake/files/BCells_harmony_info.txt"
# seurat_pca = "output/allSamples_nomt_BCells_snakemake/files/BCells_pca_info.txt"
# subgroup = "BCells"
# save_dir = "output/allSamples_nomt_BCells_snakemake/"
# pseudotime_output = "output/allSamples_nomt_BCells_snakemake/files/BCells_pseudotime.csv"


input_loom = snakemake.input[0]
#input_loom = "output/allSamples/files/allSamples_RNA_velo.loom"
pseudotime_output = snakemake.output[2]

cluster_png = snakemake.output[0]
cell_type_png = snakemake.output[1]
subgroup=snakemake.params[0]
save_dir = snakemake.params[1]
#output_png = "output/allSamples/images/rna_velocity_allsamples_harmony.png"

seurat_clusters = snakemake.input["cluster"]
#seurat_clusters = "output/allSamples/files/cluster_info.txt"
seurat_clusters = pd.read_csv(seurat_clusters, index_col = 0)

seurat_cells = snakemake.input["cell_type"]
#seurat_cells = "output/allSamples/files/cluster_info.txt"
seurat_cells = pd.read_csv(seurat_cells, index_col = 0)

seurat_pca = snakemake.input["pca"]
#seurat_pca = "output/allSamples/files/pca_info.txt"
seurat_pca = pd.read_csv(seurat_pca, index_col = 0)

seurat_umap = snakemake.input["umap"]
#seurat_umap = "output/allSamples/files/UMAP_info.txt"
seurat_umap = pd.read_csv(seurat_umap, index_col = 0)

seurat_harmony = snakemake.input["harmony"]
#seurat_harmony = "output/allSamples/files/harmony_info.txt"
seurat_harmony = pd.read_csv(seurat_harmony, index_col = 0)

adata = scv.read(input_loom, sparse = True, cache = True)
adata.var_names_make_unique()

adata_names = adata.obs.index.tolist()

adata_newnames = [re.sub(r"count:", "", name) for name in adata_names]
adata_newnames = [re.sub(r"x", "", name) for name in adata_newnames]

adata.obs["new_names"] = adata_newnames

adata.obs = adata.obs.set_index("new_names")

scv.utils.show_proportions(adata)
# number of unspliced seems really high. splice: 36%, unspliced:59%

seurat_cell_list = list(seurat_cells.index)

# Subset to only cells we have (THIS SHOULD PROBABLY MOVE TO AFTER THE
	# NORMALIZATION)
adata = adata[adata.obs.index.isin(seurat_cell_list)]
scv.utils.show_proportions(adata)
# Unspliced still high: 40%, unspliced: 54%, this is very consistent with what
# We saw in the 10x output.
# This paper does show high intron retention in B cells
# Dynamic changes in intron retention are tightly associated with regulation of
# splicing factors and proliferative activity during B-cell development

print(adata)

#scv.utils.cleanup(adata, clean='all')


n_top_genes = 2000
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_top_genes)
#scv.pp.filter_and_normalize(adata, min_counts=20, min_counts_u=10)
# n_top_genes = 3000
# log = True
# scv.pp.filter_genes(adata, min_counts = 20, min_counts_u = 10)
# scv.pp.normalize_per_cell(adata)
# if n_top_genes is not None:
#     scv.pp.filter_genes_dispersion(adata, n_top_genes = 1000)
# if log:
#     scv.pp.log1p(adata)

# Compute using harmony instead of pcs
seurat_harmony_ord = seurat_harmony.reindex(adata.obs.index)

harmony_coord = np.asarray(seurat_harmony_ord)

adata.obsm["X_pca"] = harmony_coord
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# if subgroup == "GCBcells":
# 	scv.tl.velocity(adata)
# elif subgroup == "BCells":
# 	print("deterministic")
# 	# scv.tl.recover_dynamics(adata)
# 	scv.tl.velocity(adata, use_raw = True)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

# Reorder cluster and cell type dfs to align with adata object
seurat_cell_ord = seurat_cells.reindex(adata.obs.index)
seurat_clust_ord = seurat_clusters.reindex(adata.obs.index)

# Add to adata object
adata.obs["cell_type"] = seurat_cell_ord["cell_type"]
print(seurat_clust_ord)
adata.obs["seurat_cluster"] = seurat_clust_ord["x"]

# Reorder UMAP df to align with adata object
seurat_umap_ord = seurat_umap.reindex(adata.obs.index)

umap_coord = np.asarray(seurat_umap_ord)

adata.obsm["X_umap"] = umap_coord

scv.tl.velocity_embedding(adata, basis='umap')

trans_matrix = scv.tl.transition_matrix(adata)

scv.tl.velocity_confidence(adata)

scv.pl.scatter(adata, color='velocity_confidence', perc=[2,98])



# adata.uns["cell_type_colors"] = ("#E67A7B", "#780607", "#4DAF4A", "#FF7F00", "#A65628", "#999999")

adata.uns["cell_type_colors"] = ("#984EA3", "#4DAF4A")
adata.uns["seurat_cluster_colors"] = ("#E41A1C", "#79577C", "#3C899E", "#49A75A",
"#6F8273", "#9F5196", "#DF6F32", "#FFA60F", "#FFF52F", "#CFA42D", "#B25C3F",
"#E4779C", "#D28AB0", "#999999")


scv.pl.velocity_embedding_stream(adata, basis = "umap", color = "cell_type",
	legend_loc = "on data")

plt.savefig(cell_type_png)

scv.pl.velocity_embedding_stream(adata, basis = "umap", color = "seurat_cluster",
	legend_loc = "on data")

plt.savefig(cluster_png)

scv.tl.rank_velocity_genes(adata, match_with='clusters', resolution=.4)

pd.DataFrame(adata.uns['rank_velocity_genes']['names']).head()

scv.tl.velocity_pseudotime(adata)


scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot')

plt.savefig(save_dir + "images/velocity_pseudotime.pdf")

adata.obs.to_csv(pseudotime_output)

# with PdfPages("output/allSamples_GCBcells_snakemake/images/GCBcells_rna_velocity_images.pdf") as pdf:
# 	plt.figure(None,(8,8))
# 	scv.pl.velocity_embedding(adata, basis='umap', arrow_length=1.2, arrow_size=1.2, dpi=150, color = "cell_type")
# 	pdf.savefig()
# 	plt.close()
# 	plt.figure(None,(8,8))
# 	scv.pl.velocity_embedding(adata, basis='umap', arrow_length=1.2, arrow_size=1.2, dpi=150)
# 	pdf.savefig()
# 	plt.close()
# 	plt.figure(None,(8,8))
# 	scv.pl.scatter(adata, color='velocity_confidence', perc=[2,98])
# 	pdf.savefig()
# 	plt.close()
# 	var_names = ['TUBA1C', 'CD83', 'MAP3K5', 'NFKBIA']
# 	for i in var_names:
# 		plt.figure(None,(8,8))
# 		scv.pl.velocity(adata, var_names=i, colorbar=True)
# 		pdf.savefig()
# 		plt.close()
# 	plt.figure(None,(8,8))
# 	scv.pl.velocity_graph(adata)
# 	pdf.savefig()
# 	plt.close()