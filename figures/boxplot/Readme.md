# Boxplot Figure
R script and example data for boxplot figure

**Figure:**
This functions creates a boxplot for all, effector ans secreted with the old annotation based on the effector colum (EFF) and the secreted column (Yes)

![Figure example](https://github.com/hzi-bifo/dndsAnalysis/tree/master/project_specific/colletotrichum/figures/boxplot/figure.png)

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

# creates figure
boxplot_figure <- csep_figure(c_patho, c_t, fuNOG_annotation)
print(boxplot_figure)
```