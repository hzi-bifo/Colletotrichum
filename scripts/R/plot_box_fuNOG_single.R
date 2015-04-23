#' A DndsAnalysis Function
#'
#' This functions creates a boxplot for one dataset based on a given fuNOG annotation. 
#' This boxplod is sorted by descending dN/dS ratio und divided into all, effector and secreded subsets.  
#' 
#' @param set A dataset with the columns name ratio, desc, sum_pN, sum_pS, fdr, pvalue, stars
#' @param fuNOG_annotation A path to fuNOG annotation file with columns gfam and fuNOG. 
#' This file are not filtered by accuracy so it must filtered before
#' @param fuNOG_lvl A character with fuNOG lvl of interest. e.g. "lvl1", or "lvl2"
#' @param remove Logical, if ture the fuNOG annotations "NA" and "
#' @keywords dnds
#' @export
#' @examples first_graphic <- plot_box_fuNOG_single(set1, fuNOG_good, "lvl1", remove=T, custom_annotation = T)
#' plot_box_fuNOG_single()
plot_box_fuNOG_single <- function(set,
                                  fuNOG_annotation=fuNOG_good, 
                                  fuNOG_lvl="lvl2", 
                                  remove=F,
                                  custom_annotation=T,
                                  sec_path="DndsAnalysis/annotation/secreted_genes_v_gfam.txt",
                                  eff_path="DndsAnalysis/annotation/csep_genes_v_gfam.txt"){
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
  d <- ggplot(df, aes(x= reorder(lvl, -ratio, median, order=TRUE),y=ratio, fill=subset))
  #d <- ggplot(df, aes(factor(lvl),ratio, fill=subset))
  d <- d + geom_boxplot( outlier.size = 1) #+ geom_text(data = label.df, label = "***")
  d <- d + theme_bw() + coord_flip()      
  if(new_eff){
    d <- d + ggtitle("New sec/eff annotation")
  }
  
  d <- d + theme(plot.title = element_text(size=15, vjust=1, lineheight=0.6))
  d <- d + ggtitle(paste("fuNOG", fuNOG_lvl))
  d <- d + xlab(paste("fuNOG categories on", fuNOG_lvl)) +  ylab("ratio")
  
  
  # test
  
  # test CSEP(CPAS) vs all(CPAS)
#  case_pN <-sum(as.numeric(as.matrix(df[which(df$lvl=="Cellular processes and signaling" & df$subset=="CSEP"),]$sum_pN)))
#  case_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Cellular processes and signaling" & df$subset=="CSEP"),]$sum_pS)))
  
#  all_pN <- sum(as.numeric(as.matrix(df[which(df$lvl=="Cellular processes and signaling" & df$subset=="all"),]$sum_pN)))
#  all_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Cellular processes and signaling" & df$subset=="all"),]$sum_pS)))

#  test_mat <- matrix(c(case_pN, case_pS, 
#                       all_pN,all_pS),
#                     nrow = 2)
  
 # test_case <- fisher.test(test_mat, alternative = "greater")
  
  # test CSEP(INFO) vs all(INFO)
#  case_pN <-sum(as.numeric(as.matrix(df[which(df$lvl=="Information storage and processing" & df$subset=="CSEP"),]$sum_pN)))
#  case_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Information storage and precessing" & df$subset=="CSEP"),]$sum_pS)))
  
#  all_pN <- sum(as.numeric(as.matrix(df[which(df$lvl=="Information storage and precessing" & df$subset=="all"),]$sum_pN)))
#  all_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Information storage and precessing" & df$subset=="all"),]$sum_pS)))
  
#  test_mat <- matrix(c(case_pN, case_pS, 
#                       all_pN,all_pS),
#                     nrow = 2)
  
#  test_case <- fisher.test(test_mat, alternative = "greater")


# test CSEP(defense) vs all(defense)
  #case_pN <-sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="CSEP"),]$sum_pN)))
  #case_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="CSEP"),]$sum_pS)))

 # all_pN <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="all"),]$sum_pN)))
#  all_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="all"),]$sum_pS)))

#  test_mat <- matrix(c(case_pN, case_pS, 
#                       all_pN,all_pS),
#                     nrow = 2)

 # test_case <- fisher.test(test_mat, alternative = "greater")



# test CSEP(turnover) vs all(turnover)
#case_pN <-sum(as.numeric(as.matrix(df[which(df$lvl=="Posttranslational modification, protein turnover, chaperones" & df$subset=="CSEP"),]$sum_pN)))
#case_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="CSEP"),]$sum_pS)))

#all_pN <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="all"),]$sum_pN)))
#all_pS <- sum(as.numeric(as.matrix(df[which(df$lvl=="Defense mechanisms" & df$subset=="all"),]$sum_pS)))#
#
#test_mat <- matrix(c(case_pN, case_pS, 
#                     all_pN,all_pS),
#                   nrow = 2)
#
#test_case <- fisher.test(test_mat, alternative = "greater")

  return(d)  
}





