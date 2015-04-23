#' A DndsAnalysis Function
#'
#' This functions creates a barplot for one dataset (e.g. set1/set2). 
#' It will takes the top_numer significant gene familiy descriptions
#' order_by c("ratio","fdr")
#' @param  
#' @keywords plot
#' @export
#' @examples  
#' plot_bar_single()
plot_bar_single <- function(set, 
                            top_number=60,
                            order_by="ratio",
                            mean_if_mult=T){
  require(ggplot2)
  
  # a function that overwrites the ratio and fdr for gfams with same desc
  meanDu <- function(in_mat){
    out_mat <- in_mat
    # convert to numeric from factor
    #  out_mat[, c(2:4,12,13)] <- sapply(out_mat[, c(2:4,12,13)], as.numeric(as.matrix()))
    
    out_mat$sum_pN = as.numeric(as.character(out_mat$sum_pN))
    out_mat$sum_pS = as.numeric(as.character(out_mat$sum_pS))
    out_mat$ratio = as.numeric(as.character(out_mat$ratio))
    out_mat$fdr = as.numeric(as.character(out_mat$fdr))
    
    
    correct_this <- as.character(as.matrix(unique(in_mat[which(duplicated(in_mat$desc)),]$desc)))
    for(j in correct_this){
      #  cat(paste(j, "\n"))
      matches <- which(as.character(as.matrix(in_mat$desc)) == j)
      #  for (k in matches){
      out_mat[matches,]$ratio <-  median(as.numeric(as.matrix(in_mat[matches,]$ratio)),na.rm=T)
      out_mat[k,]$sum_pN <- median(as.numeric(as.matrix(in_mat[matches,]$sum_pN)),na.rm=T)
      out_mat[k,]$sum_pS <- median(as.numeric(as.matrix(in_mat[matches,]$sum_pS)),na.rm=T)  
      out_mat[k,]$fdr <- median(as.numeric(as.matrix(in_mat[matches,]$fdr)),na.rm=T)  
      
      # }
      
    }
    return(out_mat)
  }
    
  # create subsets for efector, secreted, all
  effector_subset <- set[which(set$effector == "EFF"),]
  secreted_subset <- set[which(set$secreted == "Yes"),]
  no_subset <- set
  no_subset$subset_type <- "all"
  effector_subset$subset_type <- "effector"
  secreted_subset$subset_type <- "secreted"
  
  # create a data frame for plotting
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
  
 
  
  # if a desc exists multiple times in the data frame, we can pool the same desc together and use the mean ratio
  if(mean_if_mult){
    cat(paste("calculate mean values for multiple desc ... "))
    df <- meanDu(df)
    # now we can work with unique desc
    df <- unique(df)
    cat(paste("done!","\n"))
  }  
  
  # order and remove ""
  if(order_by=="ratio"){
    df <- df[with(df, order(-ratio,fdr)), ]      
  } 
 # if(order_by=="fdr"){
#    df <- df[with(df, order(fdr,-ratio)), ]      
#  } 
  df <- df[which(df$desc !="character(0)"),]
  
  # check if this is correctly sorted
  #head(df)
  
  # create plot ordered by fdr
  cat(paste("create plots ... "))
  if(order_by=="ratio"){
    d <- ggplot(df[1:top_number,], aes(reorder(factor(desc), ratio, FUN=median) ,ratio, fill=fdr))
  }
#  if(order_by=="fdr"){
#    d <- ggplot(df[1:top_number,], aes(reorder(factor(desc), fdr, FUN=median) ,ratio, fill=fdr))    
#  }
  d <- d + geom_bar(stat="identity", position="dodge")  #+ geom_text(data = label.df, label = "***")
  d <- d + theme_bw() + coord_flip() #+ scale_fill_manual(values=cbPalette) 
  d <- d + ggtitle(paste("Boxplot ordered by",order_by))
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + xlab("ω ratio") +  ylab("number of protein families")
  cat(paste("done!","\n"))
  return(d)  
}
