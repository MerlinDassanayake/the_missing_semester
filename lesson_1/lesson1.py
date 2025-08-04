import scanpy as sc
import matplotlib.pyplot as plt 
import leidenalg
import igraph
import numpy as np


adata = sc.datasets.pbmc3k()
adata.raw = adata # always keep the raw data

# annotate MT genes
adata.var["mt"] = adata.var_names.str.startswith("MT-") # works for human data

# calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# plot number of genes/cell
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts",
                       "pct_counts_mt"],
               jitter=0.4, multi_panel=True)

# Filtering cells that express more than 200 genes but less than 2000
# Filtering percent mitochondrial count of less than 5
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 2000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Unsupervised cell clustering with PCA

# normalize total counts to 10,000
sc.pp.normalize_total(adata, target_sum=1e4)

# log transform
sc.pp.log1p(adata)

# identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125,
                            max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# scale each gene to unit variance and zero mean, clip extreme values to +- 10
sc.pp.scale(adata, max_value=10)

# perform PCA
sc.tl.pca(adata, svd_solver="arpack")

# set nearest nieghbors and compute UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# perform leiden clustering
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'], title='Leiden Clustering')

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# manually annotate clusters (using canonical markers)
cluster_annotations = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
    '4': 'NK cells',
    '5': 'FCGR3A+ Monocytes',
    '6': 'Dendritic cells',
    '7': 'Megakaryocytes', # trailing comma good practice
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
sc.pl.umap(adata, color = ['cell_type'], title = 'Annotated Cell Types')
