library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

out_path = "../output/treated_untreated/mouse_combined/figures/"

##### plot out the proportion of cells #####
samp_tab = read.csv(file.path(out_path, 'sample_tab.csv'), row.names = 1)

ct_proportion = samp_tab %>% 
  group_by(., condition) %>%
  count(., condition, Cell.Types) %>% 
  mutate(proportion = n / sum(n))

p = ggplot(ct_proportion, aes(x = Cell.Types, y = proportion, fill = condition)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_brewer(palette = 'Set2') + 
  theme_half_open() + 
  ylab("Cell Proportions") + 
  xlab("Cell Types") + 
  labs(fill='Condition') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file.path(out_path, "ct_proportion.png"), plot = p, height = 6, width = 13)
