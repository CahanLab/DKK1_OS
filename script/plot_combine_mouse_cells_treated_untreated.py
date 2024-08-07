import pandas as pd
import numpy as np 
import scanpy as sc
import os 
import matplotlib.pyplot as plt

import gseapy
import pickle

import seaborn as sns

from palettable.colorbrewer.qualitative import Dark2_7

plt.style.use('default')
sc.set_figure_params(dpi_save = 600, fontsize = 18)
input_path = "../output/treated_untreated/mouse_combined/"
adata = sc.read_h5ad(os.path.join(input_path, "mouse_scanpy.h5ad"))

output_path = "../output/treated_untreated/mouse_combined/figures"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

# plot out the UMAPs 
sc.settings.figdir = output_path # set the figure output path
adata.obs['Cell Types'] = adata.obs['cell_types_1']
sc.pl.umap(adata, color=['condition'], palette = Dark2_7.mpl_colors, frameon = False, save = 'condition_plot.png')
sc.pl.umap(adata, color=['batch'], save = 'batch_plot.png')
sc.pl.umap(adata, color=['Cell Types'], save = 'cell_types_plot.png')

##### plot out the marker genes ######
markers_dict = dict()
markers_dict['osteoblasts'] = ['Runx2', 'Col1a1', 'Bmp1']
markers_dict['macrophage'] = ['Cd68', 'Fcgr1']
markers_dict['macrophage M2'] = ['Arg1', 'Mrc1']
markers_dict['endothelial'] = ['Pecam1', 'Egfl7', 'Cldn5']
markers_dict['anti-angio endothelial'] = ['Sema3a', 'Sema3d']
markers_dict['tip endothelial'] = ['Flt1', 'Dll4', 'Kdr']
markers_dict['neutrophils'] = ['Ly6g', 'Cd177', 'S100a8', 'S100a9', 'Cd24a']
markers_dict['prolif. neutrophils'] = ['Top2a', 'Birc5']
markers_dict['muscle'] = ['Mustn1', 'Tpm1', 'Mylpf', 'Tnni2']

sc.set_figure_params(figsize=(8, 6))
sc.pl.dotplot(adata, markers_dict, groupby='cell_types_1', categories_order = ['osteoblasts', 'macrophages', 'macrophages M2', 'endothelial', 'anti-angio endothelial', 'tip endothelial', 'neutrophils 1', 'neutrophils 2', 'neutrophils 3', 'proliferative neutrophils', 'muscles', 'unknown 1', 'unknown 2'], dendrogram=False, save = 'cell_marker_genes.png')

##### plot out the healthy and tumor associated neutrophil genes #####
# https://www.biorxiv.org/content/10.1101/2023.07.13.548820v1.full

sub_adata = adata[['neutrophils' in x for x in adata.obs['cell_types_1']], :]
sub_adata = sub_adata[sub_adata.obs['cell_types_1'] != 'proliferative neutrophils', :]
markers_dict = dict()
markers_dict['healthy_enriched'] = ['Mmp8', 'Ifitm6', 'S100a6', 'Lyz2', 'Chil3', 'G0s2', 'Fpr2']
markers_dict['tumor_enriched'] = ['Cdkn1a', 'Ifitm1', 'Ccl4', 'Cd14', 'Ccl3', 'Il1b']

sc.set_figure_params(figsize=(8, 6))
sc.pl.dotplot(sub_adata, markers_dict, groupby='cell_types_1', dendrogram=False, save = 'healthy_tumor_neutrophils.png')

samp_tab = adata.obs
samp_tab.to_csv(os.path.join(output_path, "sample_tab.csv"))