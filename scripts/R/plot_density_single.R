#' A DndsAnalysis Function
#'
#' This functions creates a density plot for one dataset (e.g. set1/set2). With a column "ratio" 
#' @param  
#' @keywords plot
#' @export
#' @examples  
#' plot_density_single()
plot_density_single <- function(set){
  require(ggplot2)
  
  den_set <- as.data.frame(set)
  den_set$ratio <- as.numeric(as.matrix(set$ratio))
  
  d <- ggplot(den_set, aes(ratio)) 
  d <- d + geom_histogram(alpha=0.8,binwidth=.05) + theme_bw()
  d <- d + ggtitle("Histogram")
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + xlab("ratio") +  ylab("number of protein families")
  return(d)  
}
