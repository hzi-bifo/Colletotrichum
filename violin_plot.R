# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


violin_plot <- function(result_mat){
  require(ggplot2)
  require(scales)
  source("reverselog_trans.R")

  # generate violin plot
  f <- ggplot(result_mat, aes(x=type, y=ratio, colour=type, fill=type)) 
  f <- f + geom_violin(trim=FALSE, colour="grey50")

  
  f <- f + geom_boxplot(width=0.1, outlier.colour = "grey60", outlier.size = 3, colour="black", fill="grey30")
  f <- f + labs(x=" ",y="log10 dN/dS ratio") 
  f <- f + scale_y_continuous(breaks=c(0.1,0.01,1), trans="log1p")
  f <- f + theme_bw() + theme_minimal()

  f <- f + scale_fill_manual(values=c("#C2C2C2", "#F0DE4B", "#8A91C5"), guide = FALSE)
  return(f)
}