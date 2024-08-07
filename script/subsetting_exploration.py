import pandas as pd
import numpy as np 
import scanpy as sc
import os 

import matplotlib.pyplot as plt
plt.style.use('default')

# extract human cells and mouse cells 
folder_list = os.listdir("../data/aligned_samples")
folder_list = [x for x in folder_list if "DS_Store" not in x]

# subset out the human cells from the mouse cells -- raw 
for sample_id in folder_list: 
    print(sample_id)
    out_path = "../output/seperated_sample_matrix_raw/" + sample_id + "/"
    if os.path.isdir(out_path) == False:
        os.makedirs(out_path)

    adata = sc.read_10x_h5("../data/aligned_samples/" + sample_id + "/filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()

    meta_data = pd.read_csv("../data/aligned_samples/" + sample_id + "/analysis/gem_classification.csv", index_col=0)

    mouse_meta_data = meta_data.loc[meta_data['call'] == 'mm10', :]
    human_meta_data = meta_data.loc[meta_data['call'] == 'GRCh38', :]
    
    mouse_adata = adata[mouse_meta_data.index, :]
    mouse_adata = mouse_adata[:, mouse_adata.var['genome'] == 'mm10']
    new_genes = list(mouse_adata.var_names)
    new_genes = [x.removeprefix("mm10___") for x in new_genes]
    mouse_adata.var.index = new_genes
    mouse_adata.write_h5ad(out_path + "mouse_adata.h5ad")

    human_adata = adata[:, adata.var['genome'] == 'GRCh38']
    human_adata = human_adata[human_meta_data.index, :]
    new_genes = list(human_adata.var_names)
    new_genes = [x.removeprefix("GRCh38_") for x in new_genes]
    human_adata.var.index = new_genes
    human_adata.write_h5ad(out_path + "human_adata.h5ad")

# check some of the mouse genes across both human and mouse cells
for sample_id in folder_list: 
    out_path = "../output/mouse_genes_check/" + sample_id + "/"
    if os.path.isdir(out_path) == False:
        os.makedirs(out_path)
    adata = sc.read_10x_h5("../data/aligned_samples/" + sample_id + "/filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    adata = adata[:, adata.var['genome'] == 'mm10']
    meta_data = pd.read_csv("../data/aligned_samples/" + sample_id + "/analysis/gem_classification.csv", index_col=0)
    meta_data = meta_data.loc[adata.obs.index, :]
    adata.obs['call'] = meta_data['call']
    
    new_genes = list(adata.var_names)
    new_genes = [x.removeprefix("mm10___") for x in new_genes]
    adata.var.index = new_genes

    sc.pl.highest_expr_genes(adata, n_top=20)
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='Cd14')
    sc.pl.pca_variance_ratio(adata, log=True)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15)

    sc.tl.umap(adata)

    sc.settings.figdir = out_path
    sc.pl.umap(adata, color=['Cd14'], save = 'Cd14.png')
    sc.pl.umap(adata, color=['call'], save = 'call.png')
    sc.pl.umap(adata, color=['Hbb-bs'], save = 'Hbb.png')

