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
  eff_a$sample <- "C*"
  eff_a$type <- "non-unique effector"
    mean_a <- sum(as.numeric(as.matrix(eff_a$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(eff_a$sum_pS)), na.rm=T)
  
  # non-effector for background
  rest_a <- df_a[which(df_a$eff_new != T),]
  rest_a <- rest_a[with(rest_a, order(fdr,-ratio)), ]
  rest_a$type <- "non-effector"
  rest_a$sample <- "C*"
  rest_a$unique <- F
  
  # second set
  df_b <- as.data.frame(nonpatho_set)
  df_b$fdr <- as.numeric(as.matrix(df_b$fdr))
  df_b$ratio <- as.numeric(as.matrix(df_b$ratio))
  eff_b <- df_b[which(df_b$eff_new == T),]
  eff_b <- eff_b[with(eff_b, order(fdr,-ratio)), ]
  eff_b$type <- "non-unique effector"
  eff_b$sample <- "Ct"
  mean_b <- sum(as.numeric(as.matrix(eff_b$sum_pN)), na.rm=T)/sum(as.numeric(as.matrix(eff_b$sum_pS)), na.rm=T)

  # non-effector for background
  rest_b <- df_a[which(df_b$eff_new != T),]
  rest_b <- rest_b[with(rest_b, order(fdr,-ratio)), ]
  rest_b$type <- "non-effector"
  rest_b$sample <- "Ct"
  rest_b$unique <- F
  
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
  
  df_both_eff <- rbind(eff_a, eff_b)
  df_both_rest <- rbind(rest_a, rest_b)
  
  df_both_eff[which(df_both_eff$type =="non-unique effector" & df_both_eff$unique==T),]$type <- "unique effector"

  d <- ggplot(df_both_eff, aes(x= ratio,y=fdr, colour=type))
  d <- d + geom_point(data=df_both_rest, alpha=0.1, color="black")
  d <- d + geom_point(alpha=0.8)  + theme_bw() + scale_y_log10()
  d <- d + xlab("dN/dS ratio") + ylab("log10 p-value")
  d <- d + facet_grid(. ~ sample)
  d <- d + geom_hline(yintercept=0.01, alpha=0.2)
  
 return(d)
}