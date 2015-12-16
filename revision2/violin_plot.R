# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

violin_plot <- function(df){
  require(ggplot2)
  require(scales)
  source("reverselog_trans.R")
  all_subset <- df
  all_subset$type2 <- "All"
  csep_subset <- df[which(df$is.csep == T),]
  csep_subset$type2 <- "CSEPs"
  ssp_subset <- df[which(df$is.secreted == T),]
  ssp_subset$type2 <- "Secreted"
  
  df2 <- rbind(all_subset, ssp_subset, csep_subset)
  # reorder
  df2$type2 <- factor(df2$type2, levels = c('All','Secreted','CSEPs'))
  
 
  # generate violin plot
  f <- ggplot(df2, aes(x=type2, y=ratio, colour=type2, fill=type2)) 
  f <- f + geom_violin(trim=FALSE, colour="grey50")
  f <- f + geom_boxplot(width=0.1, outlier.colour = "grey60", outlier.size = 3, colour="black", fill="grey30")
  f <- f + labs(x=" ",y="log10 dN/dS ratio") 
  f <- f + scale_y_continuous(breaks=c(0.1,0.01,1), trans="log1p")
  f <- f + theme_bw() + theme_minimal()
  f <- f + scale_fill_manual(values=c("#C2C2C2", "#8A91C5", "#F0DE4B"), guide = FALSE)
  f <-f + geom_text(x = c(1, 1), y = c(0, 0) , label = c("***", "***") 
  return(f)
}