#' A DndsAnalysis Function
#'
#' This function runs a fishers test with fdr correction
#' @param  
#' @keywords 
#' @export
#' @examples 
#' dnds_test()
dnds_test <- function(dnds_table,
                      pN_count=2, 
                      pS_count=3,
                      all_pN,
                      all_pS){
  
  return_tab <- dnds_table
  return_tab$pvalue <- -1
  return_tab$fdr <- NULL
  
  # if no other background data is provided, compare with the sample
  if(missing(all_pN)){
    all_pN <- sum(as.numeric(as.matrix(dnds_table[,pN_count])), na.rm=T) # total number of pN in sample    
  }
  if(missing(all_pS)){
    all_pS <- sum(as.numeric(as.matrix(dnds_table[,pS_count])), na.rm=T) # total number of pN in sample
  }
  
  i <- 1
  for (i in 1:nrow(dnds_table)){
    
  #  cat(".")
    if (is.finite(dnds_table[i,]$sum_pN ) & is.finite(dnds_table[i,]$sum_pS)){
    test_mat <- matrix(c(as.numeric(as.matrix(dnds_table[i,]$sum_pN)), 
                         as.numeric(as.matrix(dnds_table[i,]$sum_pS)), 
                         as.numeric(as.matrix(all_pN)) - as.numeric(as.matrix(dnds_table[i,]$sum_pN)),
                         as.numeric(as.matrix(all_pS)) - as.numeric(as.matrix(dnds_table[i,]$sum_pS))),
                       nrow = 2,)
                       dimnames =
                         list(c("pN", "pS"),
                              c("gfam", "all"))
  
    test <- fisher.test(test_mat, alternative = "greater")
    return_tab[i,]$pvalue <- test$p.value
   # cat(paste(i,dnds_table[i,]$ratio, test$p.value, "\n"))
    } else {
      return_tab[i,]$pvalue <- -1
    }
  }
  return_tab$fdr <- p.adjust(return_tab$pvalue, method="fdr")
#  stars <- ifelse(return_tab$fdr < 0.001, "***", 
#                  ifelse(return_tab$fdr < 0.01, "**", 
#                         ifelse(return_tab$fdr < 0.05, "*", " ")))  
  return_tab$stars <- "ns"
  return_tab$stars[which(return_tab$fdr < 0.001)] <- "***"
  return_tab$stars[which(return_tab$fdr > 0.001 & return_tab$fdr < 0.01 )] <- "**"
  return_tab$stars[which(return_tab$fdr > 0.01 & return_tab$fdr < 0.05 )] <- "*"

  return(return_tab)
}
