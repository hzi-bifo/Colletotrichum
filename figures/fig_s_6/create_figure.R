#' A DndsAnalysis Function
#' This functions creates a boxplot for one dataset based on a given fuNOG annotation. 
#' This boxplod is sorted by descending dN/dS ratio und divided into all, effector and secreded subsets.  
#' @author Philipp C. Münch philipp.muench@helmholtz-hzi.de
#' @param set A dataset with the columns name ratio, desc, sum_pN, sum_pS, fdr, pvalue, stars
#' @param fuNOG_annotation A path to fuNOG annotation file with columns gfam and fuNOG. 
#' @param fuNOG_lvl A character with fuNOG lvl of interest. e.g. "lvl1", or "lvl2"
#' @param remove Logical, if ture the fuNOG annotations "NA" and "
#' @keywords dnds
create_figure <- function(set,
                                  fuNOG_annotation=fuNOG_good, 
                                  fuNOG_lvl="lvl2", 
                                  remove=F,
                                  custom_annotation=T,
                                  sec_path="data/secreted_genes_v_gfam.txt",
                                  eff_path="data/csep_genes_v_gfam.txt"){
  require(ggplot2)
  
  if(custom_annotation){
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
    no_subset$subset_type <- "all"
    effector_subset$subset_type <- "CSEP"
    secreted_subset$subset_type <- "SSP"
    
  } else {  
    # use the ct annotation
    # create effector, secreted, and all datasets
    effector_subset <- set[which(set$effector == "EFF"),]
    secreted_subset <- set[which(set$secreted == "Yes"),]
    no_subset <- set
    
    # annotate datasets
    no_subset$subset_type <- "all"
    effector_subset$subset_type <- "CSEP"
    secreted_subset$subset_type <- "SSP"
  }
  
  # create data frame for plotting
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
  
  # add fuNOG annotation
  df$lvl <- NULL
  df$lvl <- fuNOG_annotation[match(df$gfam, fuNOG_annotation$gfam),][[fuNOG_lvl]] 
  
  if(remove){
    # remove missing annotation
    df <- df[which(df$lvl !="" & df$lvl != "Poorly characterized" & df$lvl != "NA"),]
  }
  
  # create plot 
  d <- ggplot(df, aes(x= reorder(lvl, as.numeric(as.matrix(ratio)), median, order=TRUE),y=ratio, fill=subset))
  #d <- ggplot(df, aes(factor(lvl),ratio, fill=subset))
  d <- d + geom_boxplot(outlier.size = 1)
  d <- d + scale_fill_manual(name = "Sample", values = c("grey", "#FF6666", "#6699FF"))

  
 # d <- d + geom_point(aes(color=subset), alpha=0.1)
  d <- d + theme_bw() + coord_flip()      
  
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + ggtitle(paste("Selection on protein families (fuNOG", fuNOG_lvl,")"))
  d <- d + xlab(" ") +  ylab("dN/dS ratio")
 d <- d + theme(legend.justification=c(1,0), legend.position=c(1,0))
 
  return(d)  
}





