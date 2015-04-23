# Boxplot Figure
R script and example data for boxplot figure

**Figure:**
This functions creates a boxplot for all, effector ans secreted with the old annotation based on the effector colum (EFF) and the secreted column (Yes)

### Usuage ###
Please install this packages:
```R
install.packages("ggplot2") # used for plotting
```

### Generate Figure ###
You can create the files with following code:

```R
# read data
dataset <- read.csv2("data/input_data.txt", sep=";", header=T)
sec_list <- read.csv2("data/secreted_genes_v_gfam.txt", sep=";", header=T)
eff_list <- read.csv2("data/csep_genes_v_gfam.txt", sep=";", header=T) 
source("create_figure.R")

# creates figure
boxplot_figure <- create_figure(dataset, sec_list, eff_list)
print(boxplot_figure)
```