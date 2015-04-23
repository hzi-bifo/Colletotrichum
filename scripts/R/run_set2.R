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
output_folder="~/Desktop/plot_set2_recalculated/"

#### LOAD DATA ####
cat(paste("load data, this can take some time ...", "\n"))
set2_sites <- dnds_load_files(in_folder = "~/Desktop/new/set2_pipe 2/dnds/set2_pipe/orig_dim/",
                                  gap_folder =  "~/Desktop/new/set2_pipe 3/gap/set2_pipe/",
                                  in_colnames = c("all_syn", "all_nonSyn", "tip_syn", "tip_nonSyn", "S", "N", "comparisons"),
                                  dim_method=T,
                                  quietly = TRUE,
                                  comparisons=T,
                                  threshold = 0.7)
#### ANNOTATE ####
cat(paste("annotate data, this can take some time ...", "\n"))
a_set2_sites <- run_annotation(set2_sites) 

#### CLEAN ####
set2_this <- which(is.finite(a_set2_sites$sum_pN) & is.finite(a_set2_sites$sum_pS))
set2 <- a_set2_sites[set2_this,]

#### STATISTICAL ANALYSIS ####
# fisher test
cat(paste("test data, this can take some time ...", "\n"))
set2_t_set2 <- dnds_test(set2)


#### FIGURES ####
cat(paste("create plots ...","\n"))
# density plot
hist <- plot_density_single(set2_t_set2)

box_old <- plot_box_cat_single(set2_t_set2,new_eff=F)
box <- plot_box_cat_single(set2_t_set2,new_eff=T)
# boxphlot with three categories: all, effector, secreted (old annotation)
# boxplot for fuNOG levels 1-2  TODO: sort the boxplots
box_fuNOG_1_old <- plot_box_fuNOG_single(set2_t_set2, fuNOG_good, "lvl1", remove=T, new_eff = F)
box_fuNOG_2_old <- plot_box_fuNOG_single(set2_t_set2, fuNOG_good, "lvl2", remove=T, new_eff = F)

box_fuNOG_1 <- plot_box_fuNOG_single(set2_t_set2, fuNOG_good, "lvl1", remove=T, custom_annotation = T)
box_fuNOG_2 <- plot_box_fuNOG_single(set2_t_set2, fuNOG_good, "lvl2", remove=T, custom_annotation = T)

# barplot
bar_desc <- plot_bar_single(set2_t_set2, top_number = 70, order_by = "ratio")

#### SAFE ####
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
set2_t_set2$eff_new <- FALSE
set2_t_set2$sec_new <- FALSE

set2_t_set2[which(!is.na(match(as.character(as.matrix(set2_t_set2$name)), eff_list_u))),]$eff_new <- TRUE
set2_t_set2[which(!is.na(match(as.character(as.matrix(set2_t_set2$name)), sec_list_u))),]$sec_new <- TRUE
# annotate effector and sec

# create subsets
set2_effector_subset <- set2_t_set2[which(set2_t_set2$eff_new == TRUE),]
set2_noneffector_subset <- set2_t_set2[which(set2_t_set2$eff_new != TRUE),]
set2_secreted_subset <- set2_t_set2[which(set2_t_set2$sec_new == TRUE),]
set2_nonsecreted_subset <- set2_t_set2[which(set2_t_set2$sec_new != TRUE),]





