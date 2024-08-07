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

out_path = "../output/treated_untreated/mouse_combined/"
if os.path.isdir(out_path) == False:
    os.makedirs(out_path)

for sample_id in folder_list: 
    temp_adata = sc.read_h5ad("../output/filtered_sample_matrix/" + sample_id + "/filtered_mouse_adata.h5ad")
    temp_adata.obs['batch'] = sample_id
    temp_adata.obs.index = [x + "_" + sample_id for x in temp_adata.obs.index]
    if big_adata.shape[1] == 0: 
        big_adata = temp_adata
    else:
        big_adata = big_adata.concatenate(temp_adata, index_unique=None, batch_key=None)

big_adata.var = big_adata.var.iloc[:, 0:4]
big_adata.var['mt'] = big_adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(big_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.settings.figdir = out_path # set the figure output path
sc.pl.violin(big_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save = 'human_basic_statistics.png')

# filter out genes with really sparse expression pattern 
# at least have expression in 100 cells. 
sc.pp.filter_genes(big_adata, min_cells=100)

sc.pp.normalize_total(big_adata, target_sum=1e4)
sc.pp.log1p(big_adata)

big_adata.raw = big_adata
sc.pp.highly_variable_genes(big_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

big_adata = big_adata[:, big_adata.var.highly_variable]
sc.pp.scale(big_adata, max_value=10)
sc.tl.pca(big_adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(big_adata, log=True, save='elbow_plot.png')

sc.pp.neighbors(big_adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(big_adata)
sc.pl.umap(big_adata, color=['batch'], save = 'batch_info.png')
sc.pl.umap(big_adata, color=['n_genes_by_counts'], save = 'n_genes.png')

big_adata.obs['condition'] = None
big_adata.obs.loc[big_adata.obs['batch'].isin(['F43N', 'F43R']), 'condition'] = 'Untreated'
big_adata.obs.loc[big_adata.obs['batch'].isin(['OSW3', 'OSW4', 'OSW5', 'OSW6']), 'condition'] = 'Treated'
sc.pl.umap(big_adata, color=['condition'], save = 'conditions.png')

sc.tl.leiden(big_adata, resolution = 0.5)
sc.pl.umap(big_adata, color=['leiden'], save = 'leiden.png')

sc.tl.leiden(big_adata,.15, restrict_to=["leiden",["13"]])
sc.pl.umap(big_adata, color=['leiden_R'], save = 'leiden_R.png')

sc.tl.rank_genes_groups(big_adata, 'leiden_R', method='wilcoxon', pts = True)

if os.path.isdir(os.path.join(out_path, "leiden_marker_genes")) == False: 
    os.makedirs(os.path.join(out_path, "leiden_marker_genes"))

for temp_leiden in np.unique(big_adata.obs['leiden_R']):
    gene_rank = sc.get.rank_genes_groups_df(big_adata, group=temp_leiden)
    gene_rank.to_csv(os.path.join(out_path, "leiden_marker_genes", temp_leiden + "_markergenes.csv"))

# cluster 10, 14, 7, 8 marker genes -- osteoblasts 
sc.pl.umap(big_adata, color = ['Runx2', 'Col1a1', 'Bmp1', 'Spp1']) 

# cluster 4 marker genes -- macrophage M2
# https://www.nature.com/articles/s41467-019-12384-2#:~:text=Macrophage%20M2%20polarization%20involves%20tyrosine,%2C%20Fizz1)%2C%20chitinase%2Dlike
sc.pl.umap(big_adata, color = ['Cd68', 'Arg1', 'Mrc1'])

# cluster 9 marker genes -- macrophage
sc.pl.umap(big_adata, color = ['Cd68', 'Cxcl16', 'Fcgr1', 'Cd38'])

# cluster 11 -- endothelial 1
sc.pl.umap(big_adata, color = ['Pecam1', 'Ccl21a', 'Egfl7', 'Cldn5']) 

# https://link.springer.com/article/10.1007/s10456-011-9251-z
# Cd34 -- tip cells 
# Dll4 - angiogeneiss 

# https://www.sciencedirect.com/science/article/pii/B9780128023853000061
# Kdr -- VEGFR2


# cluster 12 -- tip endothelial
sc.pl.umap(big_adata, color = ['Pecam1', 'Ccl21a', 'Egfl7', 'Cldn5', 'Cd34', 'Dll4', 'Kdr']) 

# https://www.frontiersin.org/articles/10.3389/fimmu.2020.00346/full#:~:text=Similar%20to%20their%20function%20in,and%20Sema3F%20(Table%201).
#  this is probably a anti-angio endothelial cells 
# cluster 6 -- anti-angio endothelial
sc.pl.umap(big_adata, color = ['Pecam1', 'Ccl21a', 'Egfl7', 'Cldn5', 'Sema3a', 'Sema3d']) 

# cluster 16 -- highly proliferative neutrophils 
sc.pl.umap(big_adata, color = ['Top2a', 'Cd177', 'Ly6g', 'S100a8', 'Cd24a'])

# cluster 15 -- unknown
sc.pl.umap(big_adata, color = ['AI662270', 'Ptma', 'Rpsa', 'Dock10', 'Rpl36al', 'Rps6'])

# cluster 13 -- unknown
sc.pl.umap(big_adata, color = ['Rplp1', 'Rpl26', 'Rpl18a', 'Tpt1', 'Rpl10a'])

# cluster 0 marker genes -- neutrophil 1 
# probably cancer associated 
# https://www.biorxiv.org/content/10.1101/2023.07.13.548820v1.full

# https://academic.oup.com/jleukbio/article/94/4/585/6959607
# Ly6g -- migratory maybe? 

# Cd177 -- neutrophil activation?
# https://jlb.onlinelibrary.wiley.com/doi/full/10.1002/JLB.3A0520-081RR
# the paper above suggested that Cd177 allows for transmigration of blood to tissue 

sc.pl.umap(big_adata, color = ['Cd177', 'S100a8', 'Pglyrp1', 'Cd24a', 'Il1b', 'Ly6g'])

# cluster 1 marker genes -- neutrophil 3
# probably normal or initiating neutrophils 
# https://www.biorxiv.org/content/10.1101/2023.07.13.548820v1.full
sc.pl.umap(big_adata, color = ['Cd177', 'S100a8', 'Pglyrp1', 'Cd24a', 'Mmp8', 'S100a6'])

# cluster 5, 3, 2 -- neutrophil 2 
# probably intermeidate 
# https://www.biorxiv.org/content/10.1101/2023.07.13.548820v1.full
# Cxcr2 -- no idea what that does 
sc.pl.umap(big_adata, color = ['Cd177', 'S100a8', 'Pglyrp1', 'Cd24a', 'Mmp8', 'S100a6', 'Cxcr2'])

# cluster 13,1 -- muscle 
sc.pl.umap(big_adata, color = ['Mylpf', 'Tnni2'])
# drug is causing more neutrophil signaling in the tumor cells 
# cancer associated neutrophils are enriched in Il1b -- inhibits tumor growths 


big_adata.obs['cell_types_1'] = None
big_adata.obs.loc[big_adata.obs['leiden'].isin(['10', '14', '7', '8']), 'cell_types_1'] = 'osteoblasts'
big_adata.obs.loc[big_adata.obs['leiden'] == '4', 'cell_types_1'] = 'macrophages M2'
big_adata.obs.loc[big_adata.obs['leiden'] == '9', 'cell_types_1'] = 'macrophages'
big_adata.obs.loc[big_adata.obs['leiden'] == '11', 'cell_types_1'] = 'endothelial'
big_adata.obs.loc[big_adata.obs['leiden'] == '12', 'cell_types_1'] = 'tip endothelial'
big_adata.obs.loc[big_adata.obs['leiden'] == '6', 'cell_types_1'] = 'anti-angio endothelial'
big_adata.obs.loc[big_adata.obs['leiden'] == '16', 'cell_types_1'] = 'proliferative neutrophils'
big_adata.obs.loc[big_adata.obs['leiden'] == '15', 'cell_types_1'] = 'unknown 1'
big_adata.obs.loc[big_adata.obs['leiden_R'] == '13,0', 'cell_types_1'] = 'macrophages M2'
big_adata.obs.loc[big_adata.obs['leiden_R'] == '13,1', 'cell_types_1'] = 'muscles'
big_adata.obs.loc[big_adata.obs['leiden_R'] == '13,2', 'cell_types_1'] = 'unknown 2' 
big_adata.obs.loc[big_adata.obs['leiden'] == '1', 'cell_types_1'] = 'neutrophils 1'
big_adata.obs.loc[big_adata.obs['leiden'] == '5', 'cell_types_1'] = 'neutrophils 2'
big_adata.obs.loc[big_adata.obs['leiden'] == '3', 'cell_types_1'] = 'neutrophils 2'
big_adata.obs.loc[big_adata.obs['leiden'] == '2', 'cell_types_1'] = 'neutrophils 2'
big_adata.obs.loc[big_adata.obs['leiden'] == '0', 'cell_types_1'] = 'neutrophils 3'

sc.pl.umap(big_adata, color = "cell_types_1")


big_adata.write_h5ad(os.path.join(out_path, "mouse_scanpy.h5ad"))

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
for cell_type in ['endothelial', 'neutrophils', 'macrophages']:
    if os.path.isdir(os.path.join(out_path, cell_type + "_subclusters")) == False:
        os.makedirs(os.path.join(out_path, cell_type + "_subclusters"))
    sub_adata = big_adata[[cell_type in x for x in big_adata.obs['cell_types_1']], :]
    sc.tl.rank_genes_groups(sub_adata, 'cell_types_1', method='wilcoxon', pts = True)
    for sub_type in np.unique(sub_adata.obs['cell_types_1']):
        gene_rank = sc.get.rank_genes_groups_df(sub_adata, group=sub_type)
        gene_rank = gene_rank.loc[np.logical_or(gene_rank['pct_nz_group'] > 0.1, gene_rank['pct_nz_reference'] > 0.1), :]
        gene_rank.to_csv(os.path.join(out_path, cell_type + "_subclusters", sub_type + "_DE_genes.csv"))

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
        terms.to_csv(os.path.join(out_path, cell_type + "_subclusters", sub_type + "_enrichment.csv"))

