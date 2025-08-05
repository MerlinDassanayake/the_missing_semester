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
