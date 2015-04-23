#' A DndsAnalysis Function
#'
#' This functions creates a boxplot for two dataset based on a given fuNOG annotation. 
#' This boxplod is divided into all, effector and secreded subsets.  
#' 
#' @param set_a A dataset with the columns name ratio, desc, sum_pN, sum_pS, fdr, pvalue, stars
#' @param set_b A dataset with the columns name ratio, desc, sum_pN, sum_pS, fdr, pvalue, stars
#' @param fuNOG_annotation A path to fuNOG annotation file with columns gfam and fuNOG. 
#' This file are not filtered by accuracy so it must filtered before
#' @param fuNOG_lvl A character with fuNOG lvl of interest. e.g. "lvl1", or "lvl2"
#' @param remove Logical, if ture the fuNOG annotations "NA" and "
#' @keywords dnds
#' @export
#' @examples first_graphic <- plot_box_fuNOG_both(set1, set2, fuNOG_good, "lvl1", remove=T, custom_annotation = T)
#' plot_box_fuNOG_both()
plot_box_fuNOG_both <- function(set_a,
                                set_b,
                                fuNOG_annotation,
                                fuNOG_lvl="lvl1",
                                remove=F,
                                custom_annotation=T,
                                sec_path="DndsAnalysis/annotation/secreted_genes_v_gfam.txt",
                                eff_path="DndsAnalysis/annotation/csep_genes_v_gfam.txt"){
  require(ggplot2)
  
  if(custom_annotation){
    ### set_a
    # load secretion annotation
    a_sec_list <- read.csv2(sec_path, sep=";", header=T) # the original file and with ; as seperator
    a_eff_list <- read.csv2(eff_path, sep=";", header=T) # the original file and with ; as seperator
    a_sec_list_u <- unique(as.character(as.matrix(a_sec_list$gFam)))
    a_eff_list_u <- unique(as.character(as.matrix(a_eff_list$gFam)))
    
    set_a$eff_new <- FALSE
    set_a$sec_new <- FALSE
    set_a[which(!is.na(match(as.character(as.matrix(set_a$name)), a_eff_list_u))),]$eff_new <- TRUE
    set_a[which(!is.na(match(as.character(as.matrix(set_a$name)), a_sec_list_u))),]$sec_new <- TRUE
    # annotate effector and sec
    
    # create subsets
    a_effector_subset <- set_a[which(set_a$eff_new == TRUE),]
    a_secreted_subset <- set_a[which(set_a$sec_new == TRUE),]
    a_no_subset <- set_a
    
    # annotate
    a_no_subset$subset_type <- "all"
    a_effector_subset$subset_type <- "CSEP"
    a_secreted_subset$subset_type <- "SSP"
    
    a_no_subset$cat <- "Cpatho"
    a_effector_subset$cat <- "Cpatho"
    a_secreted_subset$cat <- "Cpatho"
    
    ### set_b
    # load secretion annotation
    b_sec_list <- read.csv2("~/Desktop/secreted_genes_v_gfam.txt", sep=";", header=T) # the original file and with ; as seperator
    b_eff_list <- read.csv2("~/Desktop/csep_genes_v_gfam.txt", sep=";", header=T) # the original file and with ; as seperator
    b_sec_list_u <- unique(as.character(as.matrix(b_sec_list$gFam)))
    b_eff_list_u <- unique(as.character(as.matrix(b_eff_list$gFam)))
    
    set_b$eff_new <- FALSE
    set_b$sec_new <- FALSE
    set_b[which(!is.na(match(as.character(as.matrix(set_b$name)), b_eff_list_u))),]$eff_new <- TRUE
    set_b[which(!is.na(match(as.character(as.matrix(set_b$name)), b_sec_list_u))),]$sec_new <- TRUE
    # annotate effector and sec
    
    # create subsets
    b_effector_subset <- set_b[which(set_b$eff_new == TRUE),]
    b_secreted_subset <- set_b[which(set_b$sec_new == TRUE),]
    b_no_subset <- set_b
    
    # annotate
    b_no_subset$subset_type <- "all"
    b_effector_subset$subset_type <- "CSEP"
    b_secreted_subset$subset_type <- "SSP"
    
    b_no_subset$cat <- "Ct"
    b_effector_subset$cat <- "Ct"
    b_secreted_subset$cat <- "Ct" 
  }
  
  # create data frame for plotting
  
  together <- rbind(a_no_subset, a_effector_subset, a_secreted_subset,
                    b_no_subset, b_effector_subset, b_secreted_subset)
  
  df <- data.frame(as.character(as.matrix(together$desc)),
                   as.numeric(as.matrix(together$ratio)),
                   as.numeric(as.matrix(together$sum_pN)),
                   as.numeric(as.matrix(together$sum_pS)),
                   as.character(as.matrix(together$name)),
                   as.numeric(as.matrix(together$fdr)),
                   as.character(as.matrix(together$stars)),
                   as.character(as.matrix(together$cat)),
                   as.character(as.matrix(together$subset_type)))
  colnames(df) <- c("desc","ratio","sum_pN","sum_pS", "gfam","fdr", "stars","cat","subset")#
  
  
  # add fuNOG annotation
  df$lvl <- NULL
  df$lvl <- fuNOG_annotation[match(df$gfam, fuNOG_annotation$gfam),][[fuNOG_lvl]] 
  
  if(remove){
    # remove missing annotation
    df <- df[which(df$lvl !="" & df$lvl != "Poorly characterized" & df$lvl != "NA"),]
  }
  
  # create plot 
  d <- ggplot(df, aes(x= reorder(lvl, -ratio, median, order=TRUE),y=ratio, fill=subset))
  d <- d + geom_boxplot(outlier.size = 1) #+ scale_fill_manual(values=cbPalette)
  d <- d + theme_bw() + coord_flip() + facet_grid(. ~ cat)      
  d <- d + xlab(paste("fuNOG categories on", fuNOG_lvl)) +  ylab("ratio")
  #d <- ggplot(df, aes(factor(lvl),ratio, fill=subset))
  #  d <- d + geom_boxplot( ) #+ geom_text(data = label.df, label = "***")
  d <- d + theme_bw() + coord_flip()      
  if(new_eff){
   # d <- d + ggtitle("New sec/eff annotation")
  }
  
  #d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  #d <- d + ggtitle(paste("fuNOG", fuNOG_lvl))
  #d <- d + xlab(paste("fuNOG categories on", fuNOG_lvl)) +  ylab("ω ratio")
  
  # here we can test all against all and cat() if we find something significant...  
  return(d)  
}



