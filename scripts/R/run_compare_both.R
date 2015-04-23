#' Compare set1 and set2. By running this file the dnds datafiles are processed. 
#' this file need following functions
#' dnds_test()
#' plot_box_fuNOG_both()
#' plot_pval()
source("DndsAnalysis/R/plot_pval.R")
# please run run_set1.R and run_set2.R before...

#### CONFIG ####
output_folder <- "~/Desktop/plot_background/"
#### TEST ####
test11 <- set1_t_set1
test22 <- set2_t_set2
test12 <- dnds_test(set1,
                         all_pN = sum(as.numeric(as.matrix(set2$sum_pN)), na.rm=T),
                         all_pS = sum(as.numeric(as.matrix(set2$sum_pS)), na.rm=T))# set 1 to backround of set2

test21 <- dnds_test(set2,
                         all_pN = sum(as.numeric(as.matrix(set1$sum_pN)), na.rm=T),
                         all_pS = sum(as.numeric(as.matrix(set1$sum_pS)), na.rm=T)) # set 2 to background of set1

#### PLOT ####

# create boxplot for comparing set1 and se2
# lvl 1
lvl1_set1vsset2_all <- plot_box_fuNOG_both(test11,test22, fuNOG_good, 
                                           fuNOG_lvl="lvl1", remove = F, custom_annotation = T)
lvl1_set1vsset2_removed <- plot_box_fuNOG_both(test11,test22, fuNOG_good, 
                                               fuNOG_lvl="lvl1", remove = T, , custom_annotation = T)
# lvl 2
lvl2_set1vsset2_all <- plot_box_fuNOG_both(test11,test22, fuNOG_good, 
                                           fuNOG_lvl="lvl2", remove = F, custom_annotation = T)
lvl2_set1vsset2_removed <- plot_box_fuNOG_both(test11,test22, fuNOG_good,
                                               fuNOG_lvl="lvl2", remove = T, custom_annotation = T)

# create boxplot with fdr corrected p-value

pval_lvl1_12 <- plot_pval(test12, fuNOG_good, fuNOG_lvl="lvl1")
pval_lvl2_12 <- plot_pval(test12, fuNOG_good, fuNOG_lvl="lvl2")
#pval_lvl3_12 <- plot_pval(test12, fuNOG_good, fuNOG_lvl="lvl3")

pval_lvl1_21 <- plot_pval(test21, fuNOG_good, fuNOG_lvl="lvl1")
pval_lvl2_21 <- plot_pval(test21, fuNOG_good, fuNOG_lvl="lvl2")
#pval_lvl3_12 <- plot_pval(test12, fuNOG_good, fuNOG_lvl="lvl3")


#### SAVE ####
# lvl1
pdf(paste(output_folder,"set1_w_set2background_lvl1_all.pdf",sep=""),width=10,height=12); print(lvl1_set1vsset2_all); dev.off()
pdf(paste(output_folder,"set1_w_set2background_lvl1_removed.pdf.pdf",sep=""),width=10,height=12); print(lvl1_set1vsset2_removed); dev.off()
# lvl2
pdf(paste(output_folder,"set1_w_set2background_lvl2_all.pdf",sep=""),width=10,height=12); print(lvl2_set1vsset2_all); dev.off()
pdf(paste(output_folder,"set1_w_set2background_lvl2_removed.pdf.pdf",sep=""),width=10,height=12); print(lvl2_set1vsset2_removed); dev.off()

pdf(paste(output_folder,"set1_backgorund_set2_lvl1.pdf",sep=""),width=6,height=6); print(pval_lvl1_12); dev.off()
pdf(paste(output_folder,"set2_backgorund_set1_lvl1.pdf",sep=""),width=6,height=6); print(pval_lvl1_21); dev.off()

pdf(paste(output_folder,"set1_backgorund_set2_lvl2.pdf",sep=""),width=6,height=7); print(pval_lvl2_12); dev.off()
pdf(paste(output_folder,"set2_backgorund_set1_lvl2.pdf",sep=""),width=6,height=7); print(pval_lvl2_21); dev.off()



