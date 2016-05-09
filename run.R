# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

options(warn=-1)

source("annotate_csep.R")
source("violin_plot.R")
source("funog_boxplot.R")

# load dndsPipe results
results_file="data/pooled/results/summary/pooled.gap_theshold_0.7.summary.txt_pval.txt"
results_data <- read.table(results_file, sep=";", header=T)

# annotate SSP and CSEP
results_data_annotated <- annotate_csep(results_data)

# generate figure 3
figure3 <- violin_plot(results_data_annotated)
pdf("figure_3.pdf")
print(figure3)
dev.off()

# generate figure 6
annotation <- read.csv2("annotation/fuNOG_annotation.txt", sep=";", header=T)
figures6 <- funog_boxplot(results_data_annotated, annotation, fuNOG_lvl = "lvl2")
pdf("figure_6.pdf", width = 8, height = 7)
print(figures6)
dev.off()

# supplementary table
results_data_annotated$fuNOG_lvl1 <- annotation[match(results_data_annotated$name, annotation$gfam),]$lvl1
results_data_annotated$fuNOG_lvl2 <- annotation[match(results_data_annotated$name, annotation$gfam),]$lvl2
results_data_annotated$fuNOG_lvl3 <- annotation[match(results_data_annotated$name, annotation$gfam),]$lvl3
write.table(results_data_annotated, file="supplementary_table.csv", row.names=F, sep=";", quote=F)

# annotate table
csep_ha <- read.table("annotation/csep_hacquard.csv",se=";", header=F)
gene2gfam <- read.table("annotation/csep_genes_v_gfam.txt", sep=";", head=T)
gene2gfam$ratio <- results_data_annotated[match(gene2gfam$gFam, results_data_annotated$name),]$ratio
csep_ha$dN.dS <- gene2gfam[match(csep_ha$V1, gene2gfam$geneID), ]$ratio
write.table(csep_ha, file="csep_ha_ratio.csv", row.names=F, sep=";", quote=F)

# mean values and fisher test
# test if dnds of CSEPs is significant higher then dnds of rest 
results_data_annotated$ratio <- as.numeric(as.matrix(results_data_annotated$ratio))
csep_subset <- results_data_annotated[which(results_data_annotated$is.csep == T),]
ssp_subset <- results_data_annotated[which(results_data_annotated$is.secreted == T),]
non_csep_subset <- results_data_annotated[which(results_data_annotated$is.csep == F),]
non_ssp_subset <- results_data_annotated[which(results_data_annotated$is.secreted == F),]
testmat_csep <- matrix(c(sum(csep_subset$sum_pN, na.rm=T), sum(non_csep_subset$sum_pN, na.rm=T), sum(csep_subset$sum_pS, na.rm=T), sum(non_csep_subset$sum_pS, na.rm=T)),
                  nrow = 2, dimnames = list(c("CSEPs", "non-CSEPs"),c("sumDn", "sumDs")))
fisher.test(testmat_csep, alternative = "greater") # p-value < 2.2e-16, OR=2.128149
cat(paste("mean all:", mean(results_data_annotated$ratio, na.rm=T), "sd:", sd(results_data_annotated$ratio, na.rm=T), "\n"))
cat(paste("mean CSEPs:", mean(csep_subset$ratio, na.rm=T), "sd:", sd(csep_subset$ratio, na.rm=T), "\n"))
cat(paste("mean Secreted:", mean(ssp_subset$ratio, na.rm=T), "sd:", sd(ssp_subset$ratio, na.rm=T), "\n"))

# test if dnds of Secreted is significant higher then dnds of rest 
testmat_ssp <- matrix(c(sum(ssp_subset$sum_pN, na.rm=T), sum(non_ssp_subset$sum_pN, na.rm=T), sum(ssp_subset$sum_pS, na.rm=T), sum(non_ssp_subset$sum_pS, na.rm=T)),
                  nrow = 2, dimnames = list(c("SSP", "non-SSP"),c("sumDn", "sumDs")))
fisher.test(testmat_ssp, alternative = "greater") 
