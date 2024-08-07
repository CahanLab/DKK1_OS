import pandas as pd
import numpy as np 
import scanpy as sc
import os 
import matplotlib.pyplot as plt

import gseapy
import pickle

import seaborn as sns

from palettable.colorbrewer.qualitative import Dark2_7, Set2_4, Set1_6

plt.style.use('default')
sc.set_figure_params(dpi_save = 600, fontsize = 18)
input_path = "../output/treated_untreated/human_combined/"
adata = sc.read_h5ad(os.path.join(input_path, "scanpy_object.h5ad"))

output_path = "../output/treated_untreated/human_combined/figures"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

if os.path.isdir(os.path.join(output_path, "marker_genes")) == False: 
    os.makedirs(os.path.join(output_path, "marker_genes"))

for temp_leiden in np.unique(adata.obs['clusters_id']):
    gene_rank = sc.get.rank_genes_groups_df(adata, group=temp_leiden)
    gene_rank.to_csv(os.path.join(output_path, "marker_genes", temp_leiden + "_markergenes.csv"))

samp_tab = adata.obs
samp_tab.to_csv(os.path.join(output_path, "treated_untreated_sample_tab.csv"))

# plot out the UMAPs 
sc.settings.figdir = output_path # set the figure output path
sc.pl.umap(adata, color=['phase'], save = 'phase_plot.png')
sc.pl.umap(adata, color=['condition'], palette = Set1_6.mpl_colors, frameon = False, save = 'condition_plot.png')
sc.pl.umap(adata, color=['batch'], save = 'batch_plot.png')
sc.pl.umap(adata, color=['clusters_id'], palette = Set2_4.mpl_colors, save = 'clusters_id_plot.png')
sc.pl.pca(adata, color = 'clusters_id', palette = Set2_4.mpl_colors, save = 'pca_clsuters_id_plot.png')
sc.pl.pca(adata, color = 'clusters_id', palette = Set2_4.mpl_colors, save = 'pca_clsuters_id_plot.png')

sc.pl.umap(adata, color = 'pca_pt', save = 'UMAP_PCA_pt.png')
# plot out interesting 

# load in the genests 
GO_genesets = gseapy.parser.get_library("GO_Biological_Process_2021")
for temp_key in GO_genesets.keys():
    if 'osteo' in temp_key: 
        print(temp_key)

osteo_genes = GO_genesets['positive regulation of osteoblast differentiation (GO:0045669)']
gene_rank = sc.get.rank_genes_groups_df(adata, group='Untreated')
gene_rank = gene_rank.loc[np.logical_or(gene_rank['pct_nz_group'] > 0.1, gene_rank['pct_nz_reference'] > 0.1), :]
gene_rank = gene_rank.loc[gene_rank['pvals_adj'] < 0.05, :]
osteo_genes = np.intersect1d(osteo_genes, gene_rank['names'])
sc.pl.dotplot(adata, osteo_genes, groupby = 'clusters_id')
# plot out osteoblast development genes 
selected_osteo_genes = ['YAP1', 'RUNX2', 'COL1A1', 'FERMT2']
sc.pl.umap(adata, color = selected_osteo_genes, save='selected_bone_genes.png')

# plot out wnt signaling genes 
GO_genesets = gseapy.parser.get_library("GO_Biological_Process_2021")
for temp_key in GO_genesets.keys():
    if 'Wnt' in temp_key: 
        print(temp_key)
wnt_genes = GO_genesets['negative regulation of Wnt signaling pathway (GO:0030178)']
gene_rank = sc.get.rank_genes_groups_df(adata, group='Untreated')
gene_rank = gene_rank.loc[np.logical_or(gene_rank['pct_nz_group'] > 0.1, gene_rank['pct_nz_reference'] > 0.1), :]
gene_rank = gene_rank.loc[gene_rank['pvals_adj'] < 0.05, :]
wnt_genes = np.intersect1d(wnt_genes, gene_rank['names'])
sc.pl.dotplot(adata, wnt_genes, groupby = 'clusters_id')

selected_wnt_genes = ['AMER2', 'APOE', 'CAV1', 'DKK1', 'DKK2', 'DKK3', 'IGFBP2']
sc.pl.umap(adata, color = selected_wnt_genes, save='selected_wnt_inhibitor_genes.png')


cell_cycle_genes = open("../data/cell_cycle_genes/cell_cycle_genes.txt", "r")
s_genes = cell_cycle_genes.readline()
g2m_genes = cell_cycle_genes.readline()

s_genes = s_genes.split(", ")
g2m_genes = g2m_genes.split(", ")
g2m_genes = np.intersect1d(g2m_genes, adata.var.index)
sc.pl.dotplot(adata, g2m_genes, groupby = 'clusters_id', save='G2M_genes.png')

s_genes = np.intersect1d(s_genes, adata.var.index)
sc.pl.dotplot(adata, s_genes, groupby = 'clusters_id', save='S_genes.png')
