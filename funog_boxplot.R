

funog_boxplot <- function(set,
                          fuNOG_annotation=annotation, 
                          fuNOG_lvl="lvl2", 
                          remove=T,
                          sec_path="annotation/secreted_genes_v_gfam.txt",
                          eff_path="annotation/csep_genes_v_gfam.txt"){
  require(ggplot2)
  
    sec_list <- read.csv2(sec_path, sep=";", header=T) # the original file and with ; as seperator
    eff_list <- read.csv2(eff_path, sep=";", header=T) # the original file and with ; as seperator
    sec_list_u <- unique(as.character(as.matrix(sec_list$gFam)))
    eff_list_u <- unique(as.character(as.matrix(eff_list$gFam)))
    
    set$eff_new <- FALSE
    set$sec_new <- FALSE
    set[which(!is.na(match(as.character(as.matrix(set$name)), eff_list_u))),]$eff_new <- TRUE
    set[which(!is.na(match(as.character(as.matrix(set$name)), sec_list_u))),]$sec_new <- TRUE
    # annotate effector and sec
    
    # create subsets
    effector_subset <- set[which(set$eff_new == TRUE),]
    secreted_subset <- set[which(set$sec_new == TRUE),]
    no_subset <- set
    
    # annotate
    no_subset$subset_type <- "All"
    effector_subset$subset_type <- "CSEPs"
    secreted_subset$subset_type <- "Secreted"
    
  
  # create data frame for plotting
  together <- rbind(no_subset, effector_subset, secreted_subset)
  df <- data.frame(as.numeric(as.matrix(together$ratio)),
                   as.character(as.matrix(together$name)),
                   as.numeric(as.matrix(together$sum_pN)),
                   as.numeric(as.matrix(together$sum_pS)),                  
                   as.numeric(as.matrix(together$fdr)),
                   as.character(as.matrix(together$subset_type)))
  colnames(df) <- c("ratio","gfam","sum_pN","sum_pS", "fdr","subset")
  
  
  # reorder
  df$subset <- factor(df$subset, levels = c('All','Secreted','CSEPs'))  
  
  # add fuNOG annotation
  df$lvl <- NULL
  df$lvl <- fuNOG_annotation[match(df$gfam, fuNOG_annotation$gfam),][[fuNOG_lvl]] 
  
  if(remove){
    # remove missing annotation
    df <- df[which(df$lvl !="" & df$lvl != "Poorly characterized" & df$lvl != "NA"),]
  }
  
  # create plot 
  d <- ggplot(df, aes(x= reorder(lvl, as.numeric(as.matrix(ratio)), median, order=TRUE),y=ratio, fill=subset))
  d <- d + geom_boxplot(outlier.size = 1)
  d <- d + scale_fill_manual(name = "Sample", values = c("#C2C2C2", "#8A91C5", "#F0DE4B"))
  d <- d + theme_bw() + coord_flip()      
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + ggtitle(paste("Selection on protein families (fuNOG", fuNOG_lvl,")"))
  d <- d + xlab(" ") +  ylab("dN/dS ratio")
  d <- d + theme(legend.justification=c(1,0), legend.position=c(1,0))
 
  return(d)  
}