#' Process set1. By running this file the dnds datafiles are processed. 
#' this file need following functions
#' dnds_load_files()
#' run_annotation()
#' dnds_test()
#' plot_bar_single()
#' plot_box_fuNOG_single()
#' plot_box_fuNOG_both()
#' plot_density_single()
#' plot_box_cat_single()

#### CONFIG ####
output_folder= "~/Desktop/plot_set1_recalculated2/"

#### load ####
source("DndsAnalysis/R/plot_density_single.R")
source("DndsAnalysis/R/plot_bar_single.R")
source("DndsAnalysis/R/plot_box_cat_single.R")
source("DndsAnalysis/R/plot_box_fuNOG_both.R")
source("DndsAnalysis/R/plot_box_fuNOG_single.R")

#### LOAD DATA ####
cat(paste("load data, this can take some time ...", "\n"))
set1_sites <- dnds_load_files(in_folder =  "~/Desktop/set1_pipe_re/dnds/set1_pipe_re/orig_dim2/",
                                  gap_folder =  "~/Desktop/set1_pipe_re/gap/set1_pipe_re/",
                                  in_colnames = c("all_syn", "all_nonSyn", "tip_syn", "tip_nonSyn", "S", "N", "comparisons"),
                                  dim_method=T,
                                  quietly = TRUE,
                                  comparisons=T,
                                  threshold = 0.7)

#### ANNOTATE DATA ####
cat(paste("annotate data, this can take some time ...", "\n"))
a_set1_sites <- run_annotation(set1_sites) 

#### CLEAN DATA ####
set1_this <- which(is.finite(a_set1_sites$sum_pN) & is.finite(a_set1_sites$sum_pS))
set1 <- a_set1_sites[set1_this,]

#### STATISTICAL ANALYSIS ####
# perform one sided fisher test
cat(paste("test data, this can take some time ...", "\n"))
set1_t_set1 <- dnds_test(set1)

#### CREATE PLOTS
cat(paste("create plots ...","\n"))
# density plot
hist <- plot_density_single(set1_t_set1)

box_old <- plot_box_cat_single(set1_t_set1,new_eff=F)
box <- plot_box_cat_single(set1_t_set1,new_eff=T)
# boxphlot with three categories: all, effector, secreted (old annotation)
# boxplot for fuNOG levels 1-2  TODO: sort the boxplots
box_fuNOG_1_old <- plot_box_fuNOG_single(set1_t_set1, fuNOG_good, "lvl1", remove=T, custom_annotation = F)
box_fuNOG_2_old <- plot_box_fuNOG_single(set1_t_set1, fuNOG_good, "lvl2", remove=T, custom_annotation = F)

box_fuNOG_1 <- plot_box_fuNOG_single(set1_t_set1, fuNOG_good, "lvl1", remove=T, custom_annotation = T)
box_fuNOG_2 <- plot_box_fuNOG_single(set1_t_set1, fuNOG_good, "lvl2", remove=T, custom_annotation = T)

# barplot
bar_desc <- plot_bar_single(set1_t_set1, top_number = 70, order_by = "ratio")

#### SAVE ####
cat(paste("save plots to ",output_folder," ...","\n",sep=""))
pdf(paste(output_folder,"histogram.pdf",sep=""),width=7,height=6); print(hist); dev.off()
pdf(paste(output_folder,"boxplot_sec_eff_all.pdf",sep=""),width=7,height=6); print(box); dev.off()
pdf(paste(output_folder,"boxplot_sec_eff_all_old_annotation.pdf",sep=""),width=7,height=6); print(box_old); dev.off()
pdf(paste(output_folder,"boxplot_fuNOG_lvl1_sec_eff_all.pdf",sep=""),width=9,height=9); print(box_fuNOG_1); dev.off()
pdf(paste(output_folder,"boxplot_fuNOG_lvl2_sec_eff_all.pdf",sep=""),width=9,height=12); print(box_fuNOG_2); dev.off()
pdf(paste(output_folder,"boxplot_fuNOG_lvl1_sec_eff_all_old_annoataion.pdf",sep=""),width=9,height=9); print(box_fuNOG_1_old); dev.off()
pdf(paste(output_folder,"boxplot_fuNOG_lvl2_sec_eff_all_old_annotation.pdf",sep=""),width=9,height=12); print(box_fuNOG_2_old); dev.off()
pdf(paste(output_folder,"barplot_des_by_ratio.pdf",sep=""),width=9,height=12); print(bar_desc); dev.off()




sec_list <- read.csv2("~/Desktop/secreted_genes_v_gfam.txt", sep=";", header=T) # the original file and with ; as seperator
eff_list <- read.csv2("~/Desktop/csep_genes_v_gfam.txt", sep=";", header=T) # the original file and with ; as seperator
sec_list_u <- unique(as.character(as.matrix(sec_list$gFam)))
eff_list_u <- unique(as.character(as.matrix(eff_list$gFam)))
set1_t_set1$eff_new <- FALSE
set1_t_set1$sec_new <- FALSE

set1_t_set1[which(!is.na(match(as.character(as.matrix(set1_t_set1$name)), eff_list_u))),]$eff_new <- TRUE
set1_t_set1[which(!is.na(match(as.character(as.matrix(set1_t_set1$name)), sec_list_u))),]$sec_new <- TRUE
# annotate effector and sec

# create subsets
set1_effector_subset <- set1_t_set1[which(set1_t_set1$eff_new == TRUE),]
set1_noneffector_subset <- set1_t_set1[which(set1_t_set1$eff_new != TRUE),]

set1_secreted_subset <- set1_t_set1[which(set1_t_set1$sec_new == TRUE),]
set1_nonsecreted_subset <- set1_t_set1[which(set1_t_set1$sec_new != TRUE),]

