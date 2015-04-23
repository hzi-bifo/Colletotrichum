# CSEP Figure
R script and example data for CSEP figure

Author: philipp.muench@helmholtz-hzi.de

**Figure:**
This functions creates a boxplot for one dataset based on a given fuNOG annotation. 
This boxplod is sorted by descending dN/dS ratio und divided into all, effector and secreded subsets.  

![Alt text](figure.jpeg?raw=true "dN/dS ratio as a function of the log10 p-value")

### Usuage ###
Please install this packages:
```R
install.packages("ggplot2") # used for plotting
```

### Generate Figure ###
You can create the files with following code:

```R
# read data

c_patho <- read.csv2("data/data_cpatho.csv", sep=";", header=T)
c_t <- read.csv2("data/data_ct.csv", sep=";", header=T)
annotation <- read.csv2("data/fuNOG_annotation.txt", sep=";", header=T)
source("create_figure.R")

# creates figure
fuNOG_figure <- create_figure(c_patho, c_t, annotation)
print(csep_figure)
```