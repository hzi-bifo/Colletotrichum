# Figure 3b
R script and example data to generate figure 3b

![Alt text](figure.jpeg?raw=true "fig3b")

*CSEP gene families with evidence of positive selection (grey line:
p=0.05) in five pathogenic Colletotrichum species (left) compared to five C. tofieldiae isolates (right).
The CSEPs that are shared or unique to each set are highlighted in red and blue, respectively.*


### Generate figure ###
Please install this packages:
```R
install.packages("ggplot2")
```

You can create the files with following code:

```R
# read data

c_patho <- read.csv2("data/set_a.txt", sep=";", header=T)
c_t <- read.csv2("data/set_b.txt", sep=";", header=T)
annotation <- read.csv2("data/annotation_full.txt", header=F)
source("create_figure.R")

# creates figure
csep_figure <- create_figure(c_patho, c_t, annotation)
print(csep_figure)
```
