# CSEP Figure
R script and example data for figure 3b

**Figure:**
CSEP gene families with evidence of positive selection (grey line:
p=0.05) in five pathogenic Colletotrichum species (left) compared to five C. tofieldiae isolates (right).
The CSEPs that are shared or unique to each set are highlighted in red and blue, respectively.

![Alt text](figure.jpeg?raw=true "fig3b")

### Usuage ###
Please install this packages:
```R
install.packages("ggplot2")
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
csep_figure <- create_figure(c_patho, c_t, annotation)
print(csep_figure)
```
