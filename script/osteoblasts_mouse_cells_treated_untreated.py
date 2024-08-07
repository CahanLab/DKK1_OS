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

out_path = "../output/treated_untreated/mouse_combined/"

big_adata = sc.read_h5ad(os.path.join(out_path, "mouse_scanpy.h5ad"))

##### look at the mouse to human conversion ######
from gseapy import Biomart
bm = Biomart()
# note the dataset and attribute names are different
m2h = bm.query(dataset='mmusculus_gene_ensembl',
               attributes=['ensembl_gene_id','external_gene_name',
                           'hsapiens_homolog_ensembl_gene',
                           'hsapiens_homolog_associated_gene_name'])
m2h = m2h.loc[m2h['hsapiens_homolog_associated_gene_name'].isna() == False, :].copy()
m2h_dict = dict()
for i in m2h.index: 
    m2h_dict[m2h.loc[i, 'external_gene_name']] = m2h.loc[i, 'hsapiens_homolog_associated_gene_name']

##### find the transcriptional differences between endothelial 1 and 2 ##### 
cell_type = 'osteoblasts'
if os.path.isdir(os.path.join(out_path, cell_type + "_subclusters_treated_untreated")) == False:
    os.makedirs(os.path.join(out_path, cell_type + "_subclusters_treated_untreated"))
sub_adata = big_adata[big_adata.obs['cell_types_1'] == 'osteoblasts', :]
sc.tl.rank_genes_groups(sub_adata, 'condition', method='wilcoxon', pts = True)

for sub_type in np.unique(sub_adata.obs['condition']):
    gene_rank = sc.get.rank_genes_groups_df(sub_adata, group=sub_type)
    gene_rank = gene_rank.loc[np.logical_or(gene_rank['pct_nz_group'] > 0.1, gene_rank['pct_nz_reference'] > 0.1), :]
    gene_rank.to_csv(os.path.join(out_path, cell_type + "_subclusters_treated_untreated", sub_type + "_DE_genes.csv"))

    gene_rank.sort_values(by=['scores'], inplace=True, ascending=False)
    gene_set_names = gseapy.get_library_name(organism='Mouse')

    gene_rank = gene_rank.loc[:, ['names', 'scores']]
    for i in gene_rank.index:
        if gene_rank.loc[i, 'names'] in m2h_dict.keys():
            gene_rank.loc[i, 'names'] = m2h_dict[gene_rank.loc[i, 'names']]
    gene_rank.index = gene_rank['names']
    gene_rank = gene_rank.drop("names", axis = 1)
    res = gseapy.prerank(rnk = gene_rank, 
                        permutation_num=1000,
                        gene_sets='GO_Biological_Process_2021')
    terms = res.res2d
    terms.to_csv(os.path.join(out_path, cell_type + "_subclusters_treated_untreated", sub_type + "_enrichment.csv"))

