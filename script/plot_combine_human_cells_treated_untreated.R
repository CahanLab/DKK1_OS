library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(cowplot)

out_path = "../output/treated_untreated/human_combined/figures/"

##### plot out the treated vs untreated gsea #####
gsea_results = read.csv("../output/treated_untreated/human_combined/Untreated_gsea_results_go_process.csv", row.names = 1)
gsea_results = gsea_results[gsea_results$FDR.q.val < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]

interesting_cat_list = c("interleukin-7-mediated signaling pathway (GO:0038111)", 
                         "regulation of telomere maintenance via telomere lengthening (GO:1904356)", 
                         "negative regulation of gene expression, epigenetic (GO:0045814)", 
                         "DNA methylation (GO:0006306)", 
                         "regulation of mRNA processing (GO:0050684)", 
                         "regulation of interferon-alpha production (GO:0032647)")

gsea_results[gsea_results$FDR.q.val == 0, "FDR.q.val"] = 10^-6
gsea_results = gsea_results[gsea_results$Term %in% interesting_cat_list, ]
gsea_results$logpval = -log10(gsea_results$FDR.q.val)
untreated_gsea_results = gsea_results
untreated_gsea_results$condition = 'Untreated'

# load in the Treated 1 
gsea_results = read.csv("../output/treated_untreated/human_combined/Treated 1_gsea_results_go_process.csv", row.names = 1)
gsea_results = gsea_results[gsea_results$FDR.q.val < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]

interesting_cat_list = c("cytoplasmic translation (GO:0002181)", 
                         "negative regulation of cell cycle G2/M phase transition (GO:1902750)", 
                         "aerobic respiration (GO:0009060)",
                         "Wnt signaling pathway, planar cell polarity pathway (GO:0060071)", 
                         "tumor necrosis factor-mediated signaling pathway (GO:0033209)")



gsea_results[gsea_results$FDR.q.val == 0, "FDR.q.val"] = 10^-6
gsea_results = gsea_results[gsea_results$Term %in% interesting_cat_list, ]
gsea_results$logpval = -log10(gsea_results$FDR.q.val)
treated_1_gsea_results = gsea_results
treated_1_gsea_results$condition = 'Treated 1'

# load in the Treated 2
gsea_results = read.csv("../output/treated_untreated/human_combined/Treated 2_gsea_results_go_process.csv", row.names = 1)
gsea_results = gsea_results[gsea_results$FDR.q.val < 0.05, ]
gsea_results = gsea_results[gsea_results$NES > 0, ]

interesting_cat_list = c("extracellular matrix organization (GO:0030198)", 
                         "regulation of ERBB signaling pathway (GO:1901184)", 
                         "positive regulation of osteoblast differentiation (GO:0045669)", 
                         "Rho protein signal transduction (GO:0007266)", 
                         "sprouting angiogenesis (GO:0002040)")

gsea_results[gsea_results$FDR.q.val == 0, "FDR.q.val"] = 10^-6
gsea_results = gsea_results[gsea_results$Term %in% interesting_cat_list, ]
gsea_results$logpval = -log10(gsea_results$FDR.q.val)
treated_2_gsea_results = gsea_results
treated_2_gsea_results$condition = 'Treated 2'

combined_gsea = rbind(untreated_gsea_results, treated_1_gsea_results, treated_2_gsea_results)
combined_gsea$Term = stringr::str_split_fixed(combined_gsea$Term, pattern = " \\(", n = 2)[, 1]

p <- ggplot(combined_gsea, aes(x=reorder(Term, logpval), y=logpval, fill = condition)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'Set2')) + 
  ylab("-log10(adj p-val)") + 
  xlab("GO Biological Processes") + 
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15), axis.text.y = element_text(size = 15), text = element_text(size = 15)) + 
  facet_grid(
    rows = vars(condition),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 0.2, units = "lines"),
    strip.background = element_blank(), 
  ) + 
  coord_flip()

ggsave(filename = file.path(out_path, "gsea_enrichment.png"), plot = p, width = 10, height = 8)
 
##### table the proportion plots ##### 

samp_tab = read.csv(file.path(out_path, 'treated_untreated_sample_tab.csv'), row.names = 1)

phase_proportion = samp_tab %>% 
                    group_by(., clusters_id) %>%
                    count(., clusters_id, phase) %>% 
                    mutate(proportion = n / sum(n))

p = ggplot(phase_proportion, aes(x = phase, y = proportion, fill = clusters_id)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_brewer(palette = 'Set2') + 
  theme_half_open() + 
  ylab("Cell Proportions") + 
  xlab("Phases") + 
  labs(fill='Clusters')
ggsave(file.path(out_path, "cellcycle_proportion.png"), plot = p, height = 6, width = 6)
