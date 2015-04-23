#' This functions creates CSEP figure
#' @author Philipp C. Münch, philipp.muench@helmholtz-hzi.de
#' @param 
#' @keywords dnds
#' @export
#' @examples 
create_figure <- function(patho_set,
                        nonpatho_set,
                        fuNOG_annotation=fuNOG_good){
  require(ggplot2)
  
  # first set
  df_a <- as.data.frame(patho_set)
  df_a$fdr <- as.numeric(as.matrix(df_a$fdr))
  df_a$ratio <- as.numeric(as.matrix(df_a$ratio))
  eff_a <- df_a[which(df_a$eff_new == T),]
  eff_a <- eff_a[with(eff_a, order(fdr,-ratio)), ]
  eff_a$type <- "C*"
  mean_a <- sum(as.numeric(as.matrix(eff_a$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(eff_a$sum_pS)), na.rm=T)
  
  # second set
  df_b <- as.data.frame(nonpatho_set)
  df_b$fdr <- as.numeric(as.matrix(df_b$fdr))
  df_b$ratio <- as.numeric(as.matrix(df_b$ratio))
  eff_b <- df_b[which(df_b$eff_new == T),]
  eff_b <- eff_b[with(eff_b, order(fdr,-ratio)), ]
  eff_b$type <- "Ct"
  mean_b <- sum(as.numeric(as.matrix(eff_b$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(eff_b$sum_pS)), na.rm=T)
  
  # check for unique in patho_set
  eff_a$unique <- FALSE
  i <- 1
  for (element_name in eff_a$name) {
    elem_match <- match(as.character(as.matrix(element_name)), eff_b$name)
    if(is.na(elem_match)){
      eff_a$unique[i] <- TRUE
    }
    i <- i + 1 
  }
  # check for unique in patho_set
  eff_b$unique <- FALSE
  i <- 1
  for (element_name in eff_a$name) {
    elem_match <- match(as.character(as.matrix(element_name)), eff_a$name)
    if(is.na(elem_match)){
      eff_b$unique[i] <- TRUE
    }
    i <- i + 1 
  }
  
  df_both <- rbind(eff_a, eff_b)
  df_both[which(df_both$type =="C*" & df_both$unique==T),]$type <- "C* unique"
  df_both[which(df_both$type =="C*" & df_both$unique!=T),]$type <- "C* non-unique"
  df_both[which(df_both$type =="Ct" & df_both$unique==T),]$type <- "Ct unique"
  df_both[which(df_both$type =="Ct" & df_both$unique!=T),]$type <- "Ct non-unique"

  d <- ggplot(df_both, aes(x= ratio,y=fdr, colour=type))
  d <- d + geom_point()  + theme_bw() + scale_y_log10()
  d <- d +  xlab("dN/dS ratio") + ylab("log10 p-value")
  d <- d + geom_hline(yintercept=0.05, alpha=0.2)

 return(d)
}