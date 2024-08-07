import pandas as pd
import numpy as np 
import scanpy as sc
import os 

import matplotlib.pyplot as plt
plt.style.use('default')

# extract human cells and mouse cells 
folder_list = os.listdir("../output/seperated_sample_matrix_raw")
folder_list = [x for x in folder_list if "DS_Store" not in x]

# filter out data for human cells 
for sample_id in folder_list: 
    print(sample_id)
    out_path = "../output/filtered_sample_matrix/" + sample_id + "/"
    if os.path.isdir(out_path) == False:
        os.makedirs(out_path)

    human_adata = sc.read_h5ad("../output/seperated_sample_matrix_raw/" + sample_id + "/human_adata.h5ad")
    human_adata.var['mt'] = human_adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(human_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.settings.figdir = out_path
    sc.pl.violin(human_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save = 'human_basic_statistics.png')
    
    sc.external.pp.scrublet(human_adata)
    
    human_adata = human_adata[human_adata.obs['predicted_doublet'] == False, :]
    human_adata = human_adata[human_adata.obs.n_genes_by_counts > 2000, :] # remove cells with low gene counts
    
    # this study suggest that a threshold of 10% for mitochondrial genes might be the best for human cells 
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8599307/#!po=46.4286
    human_adata = human_adata[human_adata.obs.pct_counts_mt < 10, :]

    print(human_adata.shape)
    human_adata.write_h5ad(out_path + "filtered_human_adata.h5ad")

# filter out data for mouse cells 
for sample_id in folder_list: 
    print(sample_id)
    out_path = "../output/filtered_sample_matrix/" + sample_id + "/"
    if os.path.isdir(out_path) == False:
        os.makedirs(out_path)

    mouse_adata = sc.read_h5ad("../output/seperated_sample_matrix_raw/" + sample_id + "/mouse_adata.h5ad")
    mouse_adata.var['mt'] = mouse_adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(mouse_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.settings.figdir = out_path
    sc.pl.violin(mouse_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save = 'mouse_basic_statistics.png')
    
    sc.external.pp.scrublet(mouse_adata)
    
    mouse_adata = mouse_adata[mouse_adata.obs['predicted_doublet'] == False, :]
    mouse_adata = mouse_adata[mouse_adata.obs.n_genes_by_counts > 500, :] # remove cells with low gene counts

    # for mouse cells, I am going to stay with 5% based on the publication stated above    
    mouse_adata = mouse_adata[mouse_adata.obs.pct_counts_mt < 5, :]

    print(mouse_adata.shape)
    mouse_adata.write_h5ad(out_path + "filtered_mouse_adata.h5ad")