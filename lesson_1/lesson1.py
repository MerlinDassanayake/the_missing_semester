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

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 2000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

