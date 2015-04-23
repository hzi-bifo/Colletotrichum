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


#### CONFIGURE ####
output_folder="~/Desktop/plot_pooled_recalculated3/"
dnds_folder="~/Desktop/pooled_pipe_re/dnds/pooled_pipe_re/orig_dim2/"
gap_folder="~/Desktop/pooled_pipe_re/gap/pooled_pipe_re/"

source("DndsAnalysis/R/plot_density_single.R")
source("DndsAnalysis/R/plot_bar_single.R")
source("DndsAnalysis/R/plot_box_cat_single.R")
source("DndsAnalysis/R/plot_box_fuNOG_both.R")
source("DndsAnalysis/R/plot_box_fuNOG_single.R")

#### LOAD DATA ####
cat(paste("load data, this can take some time ...", "\n"))
set_pooled_sites <- dnds_load_files(in_folder =dnds_folder,
                                  gap_folder = gap_folder,
                                  in_colnames = c("all_syn", "all_nonSyn", "tip_syn", "tip_nonSyn", "S", "N", "comparisons"),
                                  dim_method=T,
                                  quietly = TRUE,
                                  comparisons=T,
                                  threshold = 0.7)
#### ANNOTATE ####
cat(paste("annotate data, this can take some time ...", "\n"))
a_set_pooled_sites <- run_annotation(set_pooled_sites) 

#### CLEAN ####
# clean data
set_pooled_this <- which(is.finite(a_set_pooled_sites$sum_pN) & is.finite(a_set_pooled_sites$sum_pS))
set_pooled <- a_set_pooled_sites[set_pooled_this,]

#### TEST ####
# fisher test
cat(paste("test data, this can take some time ...", "\n"))
set_pooled_t_set_pooled <- dnds_test(set_pooled)
    #write.table(set1_t_set1, file=paste(output_folder,"table.csv"), sep=";",quote=F)
  
#### CREATE PLOTS ####
cat(paste("create plots ...","\n"))
# density plot
hist <- plot_density_single(set_pooled_t_set_pooled)
box_old <- plot_box_cat_single(set_pooled_t_set_pooled,new_eff=F)
box <- plot_box_cat_single(set_pooled_t_set_pooled,new_eff=T)
# boxplot with three categories: all, effector, secreted (old annotation)
# boxplot for fuNOG levels 1-2  TODO: sort the boxplots
box_fuNOG_1_old <- plot_box_fuNOG_single(set_pooled_t_set_pooled, fuNOG_good, "lvl1", remove=T, new_eff = F)
box_fuNOG_2_old <- plot_box_fuNOG_single(set_pooled_t_set_pooled, fuNOG_good, "lvl2", remove=T, new_eff = F)

box_fuNOG_1 <- plot_box_fuNOG_single(set_pooled_t_set_pooled, fuNOG_good, "lvl1", remove=T, custom_annotation = T)
box_fuNOG_2 <- plot_box_fuNOG_single(set_pooled_t_set_pooled, fuNOG_good, "lvl2", remove=T, new_eff = T)

# barplot
bar_desc <- plot_bar_single(set_pooled_t_set_pooled, top_number = 70, order_by = "ratio")
   
#### SAVE ####
# save plots as pdf to folder
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
set_pooled_t_set_pooled$eff_new <- FALSE
set_pooled_t_set_pooled$sec_new <- FALSE

set_pooled_t_set_pooled[which(!is.na(match(as.character(as.matrix(set_pooled_t_set_pooled$name)), eff_list_u))),]$eff_new <- TRUE
set_pooled_t_set_pooled[which(!is.na(match(as.character(as.matrix(set_pooled_t_set_pooled$name)), sec_list_u))),]$sec_new <- TRUE
# annotate effector and sec

# create subsets
setpooled_effector_subset <- set_pooled_t_set_pooled[which(set_pooled_t_set_pooled$eff_new == TRUE),]
setpooled_noneffector_subset <- set_pooled_t_set_pooled[which(set_pooled_t_set_pooled$eff_new != TRUE),]

setpooled_secreted_subset <- set_pooled_t_set_pooled[which(set_pooled_t_set_pooled$sec_new == TRUE),]
setpooled_nonsecreted_subset <- set_pooled_t_set_pooled[which(set_pooled_t_set_pooled$sec_new != TRUE),]

