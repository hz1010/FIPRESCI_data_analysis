import scanpy as sc
import numpy as np
adata=sc.read_h5ad("brain_ann.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)##min_disp=0.5
adata.raw = adata


adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.umap(adata)


#sc.tl.paga(adata, groups="predicted.id",copy=False)
#sc.pl.paga(adata,threshold=0.03,show=False)
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=4, use_rep='X_diffmap')#n_ne=10
sc.tl.draw_graph(adata)
sc.tl.paga(adata, groups="predicted.id",copy=False)
sc.pl.paga(adata,threshold=0.03,show=False)


sc.tl.draw_graph(adata, init_pos="paga")
sc.tl.umap(adata, init_pos='paga')
adata.write(filename="brain_pre_nb4.h5ad")
