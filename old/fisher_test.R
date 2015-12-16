fisher_test <- function(dnds_summary){
  # input: name;sum_pN;sum_pS;ratio
  # output: name;sum_pN;sum_pS;ratio;pvalue;fdr;stars
  
  # iterate over rows and run one sided fisher test with same sample as background
  options(warn=-1)
  dnds_summary$pvalue <- NULL
  dnds_summary$fdr <- NULL
  dnds_summary$stars <- "ns"
  
  for (i in 1:nrow(dnds_summary)) {
    
    TestMartix <-
      matrix(c(
        dnds_summary$sum_pN[i],
        dnds_summary$sum_pS[i], 
        sum(dnds_summary$sum_pN)-dnds_summary$sum_pN[i],
        sum(dnds_summary$sum_pS)-dnds_summary$sum_pS[i]),
        nrow = 2,
        dimnames = list(type = c("dN", "dS"),
                        sample = c("case", "control")))
    
    dnds_summary$pvalue[i] <- fisher.test(TestMartix, alternative = "greater")$p.value
    
  }
  # fdr correction
  dnds_summary$fdr <- p.adjust(dnds_summary$pvalue, method = "fdr")
  
  dnds_summary$stars  <- ifelse(dnds_summary$fdr < .001, "***",
                                ifelse(dnds_summary$fdr < .01, "** ",
                                       ifelse(dnds_summary$fdr < .05, "* ", " ")))
  return(dnds_summary)
}