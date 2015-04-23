#' A DndsAnalysis Function
#'
#' This function allows you to import dnds pathway results
#' @param Path to folder where the .DndsRatio.txt file are stored
#' @param Colnames of input data. Default: "Sd|Nd|S|N|comparisons|pS|pN|pNpS"
#' @param Path to folder where the .gap.txt files are stored
#' @param File ending of gap files. Default: ".gap.txt"
#' @param Threshold: All positions with gap value above this threshold will be removed. Default: 1
#' @param Use JC to calculate dnds
#' @param Numer of characters to remove from filename to get the protein identifyier. Default: 14
#' @param Verbose mode
#' @keywords dnds
#' @export
#' @examples dnds_run1 <- dnds_load_files()
#' dnds_load_files()

dnds_load_files <- function(in_folder = "/Volumes/meerkat.cs.uni-duesseldorf.de/metagenomics/projects/colletotrichum/dnds/pipeline/0001_run1/dnds/0001/2/",
                            in_colnames = c("Sd", "Nd", "S", "N", "comparisons", "pS", "pN","pNpS"),
                            gap_folder = "/Volumes/meerkat.cs.uni-duesseldorf.de/metagenomics/projects/colletotrichum/dnds/pipeline/0001_run1/gap/0001/",
                            gap_ending = ".gap.txt",
                            original=F,
                            threshold = 1,
                            jc = FALSE,
                            filename_ending = 14,
                            quietly = FALSE){
  # check arguments
  if(!is.character(in_folder))
    stop("Error: input path is not a character")
  if(!is.character(gap_folder))
    stop("Error: gap path is not a character")
  if(!is.numeric(threshold))
    stop("Error: Threshold is not numeric")
  if(threshold>1 || threshold<0)
    stop("Error: Threshold must be between 0 and 1")
  
  # process arguments
  file_list <- list.files(in_folder)
  
  # create results matrix
  mat = data.frame(matrix(vector(), length(file_list), 4, 
                          dimnames=list(c(), c("name", "sum_pN", "sum_pS","ratio"))), stringsAsFactors=F) # name, sum(dN), sum(dS), sum(dN/dS) 
  # iterate over every file
  i <- 1
  for(file in file_list){
    # exctract family identifier
    fam_name <- substr(as.character(file), 1, nchar(as.character(file))-filename_ending)
    if(!quietly){
      cat(paste("Load",fam_name,"\n"))
    }
    file_data <- read.table(paste(in_folder,file,sep=""), header=T)
    if (threshold < 1){
      # exclude dnds values based on gap proportion
      gap_data <- read.table(paste(gap_folder,fam_name,gap_ending,sep=""), header=F)
      file_data <- file_data[which(gap_data<threshold),]
    }
    
    if(original){
      # remove NA and NaN and Inf
      file_data <- file_data[complete.cases(file_data), , drop=FALSE]
      
      # calculate output values
      sum_pN <- sum(file_data$all_nonSyn) + sum(file_data$tip_nonSyn)
      sum_pS <- sum(file_data$all_syn) + sum(file_data$tip_syn)
    }  
    
    else  {
      # remove NA and NaN and Inf
      file_data <- file_data[complete.cases(file_data), , drop=FALSE]
      file_data <- file_data[is.finite(file_data$pN) & is.finite(file_data$pS),]
      
      # calculate output values
      sum_pN <- sum(file_data$pN)
      sum_pS <- sum(file_data$pS)
    }
    # calculate ratio
    if (jc){
      if(original){
        stop("JC calculation for original data not possible")
      } else {
        # this is not working !!!
        #  ratio <- sum(( (-3*log(1-(4*file_data$pN/3)))/4 ) / ( (-3 * log(1-(4*file_data$pS/3)))/4 ), na.rm=T) # why this is inf??
      }
    } else { # no JC
      if (sum_pS > 10){
        ratio <- sum_pN / sum_pS
      } else {
        ratio <- NA
      }
    }
    
    # add row to output matrix
    mat$name[i] <- fam_name
    mat$sum_pN[i] <- sum_pN
    mat$sum_pS[i] <- sum_pS
    mat$ratio[i] <- ratio
    i <- i + 1   
  }
  return(mat)
}
