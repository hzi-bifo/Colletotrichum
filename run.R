# check input data

source("violin_plot.R")
source("annotate_csep.R")

# load data
set <- read.csv2("/home/muench/Schreibtisch/github_old/dndsAnalysis/project_specific/colletotrichum/figures/boxplot/data/input_data.txt", sep=";", header=T)

# subset data
set_csep <- set[which(set$eff_new==T),]
set_ssp <- set[which(set$sec_new==T),]
set_non_csep <- set[which(set$eff_new!=T),]

# print mean values based on mean(ratio)
cat(paste("all", mean(as.numeric(as.matrix(set$ratio)), na.rm=T),"sd:",sd(as.numeric(as.matrix(set$ratio)), na.rm=T), "\n"))
cat(paste("csep", mean(as.numeric(as.matrix(set_csep$ratio)), na.rm=T),"sd:", sd(as.numeric(as.matrix(set_csep$ratio)), na.rm=T), "\n"))
cat(paste("ssp", mean(as.numeric(as.matrix(set_ssp$ratio)), na.rm=T),"sd:", sd(as.numeric(as.matrix(set_ssp$ratio)), na.rm=T), "\n"))
cat(paste("non-csep", mean(as.numeric(as.matrix(set_non_csep$ratio)), na.rm=T),"sd:", sd(as.numeric(as.matrix(set_non_csep$ratio)), na.rm=T),"\n"))

# print mean values based on sum(dN) / sum(dS) over sample
cat(paste("all*",sum(as.numeric(as.matrix(set$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(set$sum_pS)), na.rm=T) , "\n"))
cat(paste("csep*", sum(as.numeric(as.matrix(set_csep$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(set_csep$sum_pS)), na.rm=T), "\n"))
cat(paste("ssp*", sum(as.numeric(as.matrix(set_ssp$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(set_ssp$sum_pS)), na.rm=T), "\n"))
cat(paste("non-csep*", sum(as.numeric(as.matrix(set_non_csep$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(set_non_csep$sum_pS)), na.rm=T), "\n"))

# create df for plot
df <- data.frame(ratio=as.numeric(as.matrix(set$ratio)),
                 type=rep("All", length(as.numeric(as.matrix(set$ratio)))), 
                 stringsAsFactors=FALSE) 

df[which(set$sec_new==T),]$type <- "Secreted"
df[which(set$eff_new==T),]$type <- "CSEPs"

figure3 <- violin_plot(df)
figure3

