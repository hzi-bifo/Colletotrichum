#' A DndsAnalysis Function
#'
#' This functions creates a boxplot for one dataset based on the p-values. 
#' This boxplod is divided into all, effector and secreded subsets.  
#' 
#' @param set A dataset with the columns name ratio, desc, sum_pN, sum_pS, fdr, pvalue, stars
#' @param fuNOG_annotation A path to fuNOG annotation file with columns gfam and fuNOG. 
#' This file are not filtered by accuracy so it must filtered before
#' @param fuNOG_lvl A character with fuNOG lvl of interest. e.g. "lvl1", or "lvl2"
#' @param remove Logical, if ture the fuNOG annotations "NA" and "
#' @keywords dnds
#' @export
#' @examples 
plot_pval <- function(set,
                                fuNOG_annotation=fuNOG_good,
                                fuNOG_lvl="lvl1",
                                remove=T){
  require(ggplot2)
  
  
  df <- as.data.frame(set)
  df$fdr <- as.numeric(as.matrix(df$fdr))
  df$ratio <- as.numeric(as.matrix(df$ratio))

  # add fuNOG annotation
  df$lvl <- NULL
  df$lvl <- fuNOG_annotation[match(df$name, fuNOG_annotation$gfam),][[fuNOG_lvl]] 
    
  # remove
  
  if(remove){
  #  df <- df[which(df$lvl !="" & df$lvl != "Poorly characterized" & df$lvl != "NA"),]
  
  df <- df[which(as.character(as.matrix(df$lvl))!="General function prediction only" 
                     & as.character(as.matrix(df$lvl)) !="Function unknown" 
                     & as.character(as.matrix(df$lvl)) !="" 
                     & as.character(as.matrix(df$lvl)) != "Poorly characterized" 
                     & df$lvl != "NA"),]
  }
  
  # sort
  df <- df[with(df, order(fdr,-ratio)), ]

  # select only significant
  df <- df[which(df$fdr < 0.001),]
  
  #colnames(df) <- c("desc","ratio","sum_pN","sum_pS", "gfam","fdr", "stars","cat","subset")#
  
  # create plot 
  d <- ggplot(df, aes(x= reorder(lvl, -fdr, median, order=TRUE),y=fdr))
  d <- d + geom_boxplot()  + scale_y_log10()
  d <- d + theme_bw() + coord_flip() + 
    xlab(paste("fuNOG categories on", fuNOG_lvl)) +
    ylab("log10 pvalue") #+
   # ggtitle("Significant Ct with C* as background")
  
    return(d)  
  
  df2 <- data.frame(as.character(as.matrix(df$lvl)),
                       as.numeric(as.matrix(df$ratio)),
                    as.numeric(as.matrix(df$sum_pN)),
                    as.numeric(as.matrix(df$sum_pS)),
                       as.numeric(as.matrix(df$fdr)))
  write.table(df2, file="~/Desktop/Ct_Cpatho_background_lvl3.csv", sep=";", row.names=F, quote=F)
  
}





