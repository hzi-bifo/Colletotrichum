# CSEP Figure
R script and example data for CSEP figure

**Figure:**
dN/dS ratio as a function of the log10 p-value

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

c_patho <- read.csv2("data/set_a.txt", sep=";", header=T)
c_t <- read.csv2("data/set_b.txt", sep=";", header=T)
annotation <- read.csv2("data/fuNOG_annotation.txt", sep=";", header=T)
source("create_figure.R")

# creates figure
csep_figure <- create_figure(c_patho, c_t, fuNOG_annotation)
print(csep_figure)
```