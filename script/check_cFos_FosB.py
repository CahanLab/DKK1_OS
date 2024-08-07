import scanpy as sc
import pandas as pd
import numpy as np 
import os 
import matplotlib.pyplot as plt
from scipy.stats import ranksums

from palettable.colorbrewer.qualitative import Dark2_7, Set2_4, Set1_6

plt.style.use('default')
sc.set_figure_params(dpi_save = 600, fontsize = 18)
input_path = "../output/treated_untreated/human_combined/"
adata = sc.read_h5ad(os.path.join(input_path, "scanpy_object.h5ad"))

output_path = "../output/treated_untreated/human_combined/cFos_FosB"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

cFos = ['MMP3', 'MMP13', 'ADAMTS5', 'IL1B', 'BCL3', 'FOSL1', 'SOX9', 'NFATC1']
FosB = ['FOSB', 'TNC', 'NFKB1', 'CDK5', 'GRIA2', 'CALB1', 'MMP13', 'FOS']
total_genes = cFos + FosB

# get the normalized anndata and sample table 
norm_exp = adata.raw.to_adata().to_df()
samp_tab = adata.obs

# get the average normalized expression of cFos and FosB genes in untreated, treated 1 and treated 2
mean_exp = pd.DataFrame(columns = ['genes', 'category', 'untreated', 'treated 1', 'treated 2'])
for tmp_gene in total_genes: 
    untreated = np.mean(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Untreated', :].index, tmp_gene])
    treated_1 = np.mean(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Treated 1', :].index, tmp_gene])
    treated_2 = np.mean(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Treated 2', :].index, tmp_gene])
    if tmp_gene in cFos: 
        cat = 'cFos'
    else: 
        cat = 'FosB'
    mean_exp.loc[tmp_gene, :] = [tmp_gene, cat, untreated, treated_1, treated_2]
mean_exp.to_csv(os.path.join(output_path, 'mean_exp.csv'), index=False)

# this is where I run logfold change and wilcoxon 
stats_df = pd.DataFrame(columns = ['genes', 
                                   'category', 
                                   'treated_1_untreated_logFC', 
                                   'treated_1_untreated_pval', 
                                   'treated_2_untreated_logFC', 
                                   'treated_2_untreated_pval'])
for tmp_gene in total_genes: 
    if tmp_gene in cFos: 
        cat = 'cFos'
    else: 
        cat = 'FosB'
    untreated = np.array(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Untreated', :].index, tmp_gene])
    treated_1 = np.array(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Treated 1', :].index, tmp_gene])
    treated_2 = np.array(norm_exp.loc[samp_tab.loc[samp_tab['clusters_id'] == 'Treated 2', :].index, tmp_gene])
    
    treated_1_untreated_logFC = np.log2(np.nan_to_num(np.divide(np.mean(treated_1), np.mean(untreated)), nan = 1))
    treated_2_untreated_logFC = np.log2(np.nan_to_num(np.divide(np.mean(treated_2), np.mean(untreated)), nan = 1))

    treated_1_untreated_pval = ranksums(treated_1, untreated)[1]
    treated_2_untreated_pval = ranksums(treated_2, untreated)[1]

    stats_df.loc[tmp_gene, :] = [tmp_gene, cat, 
                                 treated_1_untreated_logFC, treated_1_untreated_pval, 
                                 treated_2_untreated_logFC, treated_2_untreated_pval]

stats_df.to_csv(os.path.join(output_path, 'stats_df.csv'), index=False)
