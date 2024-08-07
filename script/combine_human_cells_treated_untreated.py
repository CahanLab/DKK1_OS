import pandas as pd
import numpy as np 
import scanpy as sc
import os 
import matplotlib.pyplot as plt
import anndata
import scanpy.external as sce

import gseapy
import pickle
plt.style.use('default')

folder_list = os.listdir("../output/filtered_sample_matrix")
folder_list = [x for x in folder_list if "DS_Store" not in x]
folder_list = ['F43R', 'F43N', 'OSW3', 'OSW4', 'OSW5', 'OSW6']

big_adata = anndata.AnnData()

out_path = "../output/treated_untreated/human_combined/"
if os.path.isdir(out_path) == False:
    os.makedirs(out_path)

for sample_id in folder_list: 
    temp_adata = sc.read_h5ad("../output/filtered_sample_matrix/" + sample_id + "/filtered_human_adata.h5ad")
    temp_adata.obs['batch'] = sample_id
    temp_adata.obs.index = [x + "_" + sample_id for x in temp_adata.obs.index]
    if big_adata.shape[1] == 0: 
        big_adata = temp_adata
    else:
        big_adata = big_adata.concatenate(temp_adata, index_unique=None, batch_key=None)

big_adata.var = big_adata.var.iloc[:, 0:4]
big_adata.var['mt'] = big_adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(big_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.settings.figdir = out_path # set the figure output path
sc.pl.violin(big_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save = 'human_basic_statistics.png')

# filter out genes with really sparse expression pattern 
# at least have expression in 100 cells. 
sc.pp.filter_genes(big_adata, min_cells=100)

sc.pp.normalize_total(big_adata, target_sum=1e4)
sc.pp.log1p(big_adata)

cell_cycle_genes = open("../data/cell_cycle_genes/cell_cycle_genes.txt", "r")
s_genes = cell_cycle_genes.readline()
g2m_genes = cell_cycle_genes.readline()

s_genes = s_genes.split(", ")
g2m_genes = g2m_genes.split(", ")

sc.tl.score_genes_cell_cycle(big_adata, s_genes, g2m_genes, copy=False)

big_adata.raw = big_adata
sc.pp.highly_variable_genes(big_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

big_adata = big_adata[:, big_adata.var.highly_variable]
sc.pp.scale(big_adata, max_value=10)
sc.tl.pca(big_adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(big_adata, log=True, save='elbow_plot.png')

sc.pp.neighbors(big_adata, n_neighbors=10, n_pcs=10) # used to be 20  n_neighbor = 10, n_pcs = 10 is pretty good. 
sc.tl.umap(big_adata)
sc.pl.umap(big_adata, color=['batch'], save = 'batch_info.png')
sc.pl.umap(big_adata, color=['phase'], save = 'phase_plot.png')
sc.pl.umap(big_adata, color=['n_genes_by_counts'], save = 'n_genes.png')

sc.tl.leiden(big_adata, resolution = 0.09)
sc.pl.umap(big_adata, color = 'leiden', save = 'leidens.png')

big_adata.obs['condition'] = None
big_adata.obs.loc[big_adata.obs['batch'].isin(['F43N', 'F43R']), 'condition'] = 'Untreated'
big_adata.obs.loc[big_adata.obs['batch'].isin(['OSW3', 'OSW4', 'OSW5', 'OSW6']), 'condition'] = 'Treated'
sc.pl.umap(big_adata, color=['condition'], save = 'conditions.png')

big_adata.obs['clusters_id'] = None
big_adata.obs.loc[big_adata.obs['leiden'] == '0', 'clusters_id'] = 'Untreated'
big_adata.obs.loc[big_adata.obs['leiden'] == '1', 'clusters_id'] = 'Treated 2'
big_adata.obs.loc[big_adata.obs['leiden'] == '2', 'clusters_id'] = 'Treated 1'
sc.pl.umap(big_adata, color = 'clusters_id', save = 'clusters_id.png')

sc.tl.rank_genes_groups(big_adata, 'clusters_id', method='wilcoxon', pts = True)

for temp_cat in np.unique(big_adata.obs['clusters_id']):
    print(temp_cat)
    gene_rank = sc.get.rank_genes_groups_df(big_adata, group=temp_cat)
    gene_rank = gene_rank.loc[np.logical_or(gene_rank['pct_nz_group'] > 0.1, gene_rank['pct_nz_reference'] > 0.1), :]

    gene_rank.sort_values(by=['scores'], inplace=True, ascending=False)
    gene_set_names = gseapy.get_library_name(organism='Human')
    gene_rank = gene_rank.loc[:, ['names', 'scores']]
    gene_rank.index = gene_rank['names']
    gene_rank = gene_rank.drop("names", axis = 1)
    res = gseapy.prerank(rnk = gene_rank, 
                        permutation_num=1000,
                        gene_sets='GO_Biological_Process_2021')
    pickle.dump(res, open(os.path.join(out_path + temp_cat + "_gsea_object_biological_process.pickle"), 'wb'))
    terms = res.res2d
    terms.to_csv(os.path.join(out_path, temp_cat + "_gsea_results_go_process.csv"))

    res = gseapy.prerank(rnk=gene_rank, 
                        permutation_num = 1000,
                        gene_sets='MSigDB_Hallmark_2020')
    pickle.dump(res, open(os.path.join(out_path + temp_cat + "_gsea_object_MSigDB_Hallmark_2020.pickle"), 'wb'))
    terms = res.res2d
    terms.to_csv(os.path.join(out_path, temp_cat + "_gsea_results_MSigDB_Hallmark_2020.csv"))

    res = gseapy.prerank(rnk=gene_rank, 
                        permutation_num = 1000,
                        gene_sets='MSigDB_Oncogenic_Signatures')
    pickle.dump(res, open(os.path.join(out_path + temp_cat + "_gsea_object_MSigDB_Oncogenic_Signatures.pickle"), 'wb'))
    terms = res.res2d
    terms.to_csv(os.path.join(out_path, temp_cat + "_gsea_results_MSigDB_Oncogenic_Signatures.csv"))

##### make plots of the bone markers #####
sc.settings.figdir = out_path
sc.pl.umap(big_adata, color = ['SPP1', 'CALCB', 'COL1A1'], save = 'bone_marker_genes.png')
sc.pl.umap(big_adata, color = ['RUNX2', 'TWIST1', 'TGFBR1'], save = 'bone_marker_genes2.png')

##### make the poor-man pca pseudotime #####
pca_coord = big_adata.obsm['X_pca'][:, 0]
pca_coord = (pca_coord - np.min(pca_coord)) / (np.max(pca_coord) - np.min(pca_coord))
pca_coord = 1 - pca_coord
big_adata.obs['pca_pt'] = pca_coord
sc.pl.umap(big_adata, color=['pca_pt'])

for temp_cat in np.unique(big_adata.obs['clusters_id']):
    print(temp_cat)
    print(np.mean(big_adata.obs.loc[big_adata.obs['clusters_id'] == temp_cat, 'pca_pt']))
##### save the scanpy object ###### 
big_adata.write_h5ad(os.path.join(out_path, "scanpy_object.h5ad"))

##### get the dpt #####
sc.tl.diffmap(big_adata)
big_adata.uns['iroot'] = np.flatnonzero(big_adata.obs['clusters_id']  == 'Untreated')[0]

sc.tl.dpt(big_adata, n_branchings=2)



