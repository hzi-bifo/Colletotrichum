# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

violin_plot <- function(result_mat){
  require(ggplot2)
  require(scales)
  source("reverselog_trans.R")

  # change order
  
  
  # generate violin plot
  f <- ggplot(result_mat, aes(x=type, y=ratio, colour=type, fill=type)) 
  f <- f + theme_bw() + scale_y_log10() + theme_minimal()
  f <- f + geom_violin(trim=FALSE, fill="gray90", colour="grey50")
  f <- f + geom_boxplot(width=0.1, outlier.colour = "grey60", outlier.size = 3, colour="black", fill="grey30")
  f <- f + labs(x=" ",y="log 10 ratio") 

  return(f)
}