import pandas as pd
import numpy as np 
import scanpy as sc
import os 
import matplotlib.pyplot as plt
import anndata
import scanpy.external as sce

import gseapy
import pickle

import pySingleCellNet as pySCN

plt.style.use('default')

out_path = '../output/treated_untreated/tualne_pySCN_human/'
if os.path.isdir(out_path) == False: 
    os.makedirs(out_path)

sc.settings.figdir = out_path # set the figure output path

##### load in the raw data ######
folder_list = ['F43R', 'F43N', 'OSW3', 'OSW4', 'OSW5', 'OSW6']
big_adata = anndata.AnnData()

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

# load in the training data 
exp_train = pd.read_csv("../data/Tulane_data/expTab.csv", index_col = 0)
st_train = pd.read_csv("../data/Tulane_data/sampTab.csv", index_col = 0)

train_adata = sc.AnnData(exp_train.T)
st_train = st_train.loc[train_adata.obs.index, :]
train_adata.obs['cell_type'] = st_train['cell_type']
#train_adata = train_adata[train_adata.obs['cell_type'] != 'Undetermined Osteoblasts', :]

cgenes = train_adata.var.index.intersection(big_adata.var.index)
train_adata = train_adata[:, cgenes]

expTrain, expVal = pySCN.splitCommonAnnData(train_adata, ncells=300,dLevel="cell_type")
[cgenesA, xpairs, tspRF] = pySCN.scn_train(expTrain, nTopGenes = 100, nRand = 100, nTrees = 1000 ,nTopGenePairs = 100, dLevel = "cell_type", stratify=True, limitToHVG=True)

pickle.dump([cgenesA, xpairs, tspRF], open(os.path.join(out_path, "classifier_object.pickle"), "wb"))

[cgenesA, xpairs, tspRF] = pickle.load(open(os.path.join(out_path, "classifier_object.pickle"), "rb"))
# plot out the evaluation 
adVal = pySCN.scn_classify(expVal, cgenesA, xpairs, tspRF, nrand = 0)
ax = sc.pl.heatmap(adVal, adVal.var_names.values, groupby='SCN_class', cmap='viridis', dendrogram=False, swap_axes=True, save='validation_heatmap.png')

adVal.obs['cell'] = adVal.obs.index
assessment =  pySCN.assess_comm(expTrain, adVal, resolution = 0.005, nRand = 0, dLevelSID = "cell", classTrain = "cell_type", classQuery = "cell_type")
pickle.dump(assessment, open(out_path + "assessment.pickle", 'wb'))
pySCN.plot_PRs(assessment)

# plot out our data 
ourdata_class = pySCN.scn_classify(big_adata, cgenesA, xpairs, tspRF, nrand=0)
processed_adata = sc.read_h5ad("../output/treated_untreated/human_combined/scanpy_object.h5ad")
processed_adata.obs['pySCN_class'] = ourdata_class.obs['SCN_class']
processed_adata.obs = pd.concat([processed_adata.obs, ourdata_class.to_df()], axis = 1)
ourdata_class.obs['condition'] = processed_adata.obs['condition']
ourdata_class.obs['clusters_id'] = processed_adata.obs['clusters_id']

ourdata_class.write_h5ad(os.path.join(out_path, "pySCN_class_matrix.h5ad"))
processed_adata.write_h5ad(os.path.join(out_path, "pySCN_integrated_scanpy.h5ad"))

sc.set_figure_params(figsize=(6, 8), fontsize = 14)
ax = sc.pl.heatmap(ourdata_class, ourdata_class.var_names.values, groupby='clusters_id', cmap='viridis', dendrogram=False, swap_axes=True, figsize = (8, 1.3), save='condition_heatmap.png')

sc.set_figure_params(figsize=(4, 4), fontsize = 14)
sc.pl.umap(processed_adata, color = 'pySCN_class', palette = 'Dark2', save='pySCN_UMAP.png')

sc.set_figure_params(figsize=(10, 10), fontsize = 40)
sc.pl.violin(processed_adata, groupby='condition', palette = 'Set1', keys = ['Mature Osteoblasts', 'Undetermined Osteoblasts', 'Pre-osteoblasts'], jitter = False, save='pySCN_scoreS_violin.png')