#' This functions creates a boxplot for all, effector ans secreted with the old annotation  
#' based on the effector colum (EFF) and the secreted column (Yes)
#' @param set: input dataset
#' @param sec_list: secretion annotation
#' @param eff_list: effector annotation
#' @param new_eff=F take the new effector annotation DEFAULT TURE
#' @keywords dnds, plot
#' @export
#' @examples  
#' create_figure()
create_figure <- function(set, sec_list, eff_list, new_eff=T){
  require(ggplot2)
  
  if(new_eff){
    sec_list_u <- unique(as.character(as.matrix(sec_list$gFam)))
    eff_list_u <- unique(as.character(as.matrix(eff_list$gFam)))
    
    set$eff_new <- FALSE
    set$sec_new <- FALSE
    set[which(!is.na(match(as.character(as.matrix(set$name)), eff_list_u))),]$eff_new <- TRUE
    set[which(!is.na(match(as.character(as.matrix(set$name)), sec_list_u))),]$sec_new <- TRUE
    # annotate effector and sec
    # create subsets
    effector_subset <- set[which(set$eff_new == TRUE),]
    secreted_subset <- set[which(set$sec_new == TRUE),]
    no_subset <- set
    
    # annotate
    no_subset$subset_type <- "all"
    effector_subset$subset_type <- "CSEP"
    secreted_subset$subset_type <- "SSP"
    
  } else {
    # create effector, secreted, and all datasets
    effector_subset <- set[which(set$effector == "EFF"),]
    secreted_subset <- set[which(set$secreted == "Yes"),]
    no_subset <- set
    
    # annotate datasets
    no_subset$subset_type <- "all"
    effector_subset$subset_type <- "CSEP"
    secreted_subset$subset_type <- "SSP"
  }
  
  # create data frame for plotting
  together <- rbind(no_subset, effector_subset, secreted_subset)
  df <- data.frame(as.character(as.matrix(together$desc)),
                   as.numeric(as.matrix(together$ratio)),
                   as.character(as.matrix(together$name)),
                   as.numeric(as.matrix(together$sum_pN)),
                   as.numeric(as.matrix(together$sum_pS)),                  
                   as.numeric(as.matrix(together$fdr)),
                   as.character(as.matrix(together$stars)),
                   as.character(as.matrix(together$subset_type)))
  colnames(df) <- c("desc","ratio","gfam","sum_pN","sum_pS", "fdr", "stars","subset")
  
  # create plot 
  d <- ggplot(df, aes(factor(subset),ratio))
  d <- d + geom_boxplot( outlier.size = 1, fill="grey80") #+ geom_text(data = label.df, label = "***")
  d <- d + theme_bw() 
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  # if(new_eff){
  #  d <- d + ggtitle("New sec/eff annotation")}
  d <- d + xlab("category") +  ylab("ratio")
  
  # test for significance difference between all/effector all/secreted
  
  # test effector vs all
  eff_pN <- sum(as.numeric(as.matrix(effector_subset$sum_pN)),na.rm=T)
  eff_pS <- sum(as.numeric(as.matrix(effector_subset$sum_pS)),na.rm=T)
  all_pN <- sum(as.numeric(as.matrix(no_subset$sum_pN)),na.rm=T) - eff_pN
  all_pS <- sum(as.numeric(as.matrix(no_subset$sum_pS)),na.rm=T) - eff_pS
  test_mat <- matrix(c(eff_pN, eff_pS, 
                       all_pN,all_pS),
                     nrow = 2)
  
  test_eff <- fisher.test(test_mat, alternative = "greater")
  cat(paste("p-value (effector vs. all with one-sided fishers test):", test_eff$p.value, "\n"))
  
  # test secreted vs all ->>> p-value < 2.2e-16 (8.120612e-79)
  eff_pN <- sum(as.numeric(as.matrix(secreted_subset$sum_pN)),na.rm=T)
  eff_pS <- sum(as.numeric(as.matrix(secreted_subset$sum_pS)),na.rm=T)
  all_pN <- sum(as.numeric(as.matrix(no_subset$sum_pN)),na.rm=T) - eff_pN
  all_pS <- sum(as.numeric(as.matrix(no_subset$sum_pS)),na.rm=T) - eff_pS
  test_mat <- matrix(c(eff_pN, eff_pS, 
                       all_pN,all_pS),
                     nrow = 2)
  
  test_sec <- fisher.test(test_mat, alternative = "greater")
  test_sec$p.value
  cat(paste("p-value (secreted vs. all with one-sided fishers test):", test_sec$p.value, "\n"))
  
  return(d)  
}