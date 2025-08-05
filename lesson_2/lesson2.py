# Lesson 2 of the missing semester
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# load spatial data
adata = sc.datasets.visium_sge()
adata.var_names_make_unique() # optional but good practice

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat',
                            n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]

# reduce dims
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.spatial(adata, color='leiden', spot_size=100)

# Identify marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)

# spatial neighborhood enrichment
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode='moran')

# convert result to df
moran_df = pd.DataFrame(adata.uns['moranI'])

# set gene names as col (instead of index)
moran_df = moran_df.reset_index().rename(columns={"index":"gene"})

# sort by Moran's I value
moran_df_sorted = moran_df.sort_values("I", ascending=False)

# plot top 20 genes
top_n = 20
plt.figure(figsize=(10,6))
plt.barh(
    moran_df_sorted['gene'][:top_n][::-1], # reverse for descending y-axis
    moran_df_sorted['I'][:top_n][::-1]
)
plt.xlabel("Moran's I")
plt.title(f"Top {top_n} Genes by Spatial Autocorrelation")
plt.tight_layout()
plt.show()
