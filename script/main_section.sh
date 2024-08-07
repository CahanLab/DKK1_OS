# subset out the human and mouse data from the single cell
python subsetting_exploration.py 

# filter out the adta 
python filter_data.py 

# combine the data between ortho and sub data 
python combine_human_cells_treated_untreated.py 

# make prettier plots of the combine data between ortho and sub 
python plot_combine_human_cells_treated_untreated.py

# make prettier plots for the combined data between ortho and sub 
Rscript plot_combine_human_cells_treated_untreated.R

# match the tumor cells to the osteoblast developmental trajectory 
python pySCN_tulane_combine_human_cells_treated_untreated.py 

# look at the mouse cells 
python combine_mouse_cells_treated_untreated.py 

# make prettier plots for treated and untreated mouse cells 
python plot_combine_mouse_cells_treated_untreated.py

