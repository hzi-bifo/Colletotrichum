dnds_load_files <- function(in_folder = "",
                            in_colnames = c("Sd", "Nd", "S", "N", "comparisons", "pS", "pN","pNpS"),
                            gap_folder = "",
                            gap_ending = ".gap.txt",
                            original=F,
                            min_method_sites=F,
                            comparisons=F,
                            dim_method=F,
                            threshold = 1,
                            jc = FALSE,
                            filename_ending = 14,
                            quietly = TRUE){
  # check arguments
  if(!is.character(in_folder))
    stop("Error: input path is not a character")
  if(!is.character(gap_folder))
    stop("Error: gap path is not a character")
  if(!is.numeric(threshold))
    stop("Error: Threshold is not numeric")
  if(threshold>1 || threshold<0)
    stop("Error: Threshold must be between 0 and 1")
  if(original && min_method_sites && dim_method)
    stop("Error: Please use either original, dim or min method")
  # process arguments
  file_list <- list.files(in_folder)
  #file_list_gap <- list.files(gap_folder)
  # only process files with gap and dnds information
  #dnds_intlist <- strsplit(file_list, ".", fixed=T)
  #filename_dnds <- matrix(unlist(dnds_intlist),ncol=3,byrow=TRUE)[,1]
  #gap_intlist <- strsplit(file_list_gap, ".", fixed=T)
  #filename_gap <- matrix(unlist(gap_intlist),ncol=3,byrow=TRUE)[,1]  
  #onlythis <- intersect(filename_gap, filename_dnds)
  
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
        if (file.exists(paste(gap_folder,fam_name,gap_ending,sep=""))){
          gap_data <- read.table(paste(gap_folder,fam_name,gap_ending,sep=""), header=F)
        } else {
        gap_data <- rep(0,length(file_data))
        }
        
        file_data <- file_data[which(gap_data<threshold),]
      }
      if(original){
        # remove NA and NaN and Inf
        file_data <- file_data[complete.cases(file_data), , drop=FALSE]      
        # calculate output values
        sum_pN <- sum(file_data$all_nonSyn) + sum(file_data$tip_nonSyn)
        sum_pS <- sum(file_data$all_syn) + sum(file_data$tip_syn)
      } else if(min_method_sites) { # for the minimal pathway method with site correction 
        # take mean value for Sd and S based on number of comparisons
        file_data$pS <- file_data$Sd / file_data$S
        file_data$pN <- file_data$Nd / file_data$N
        # sum without NA and without Inf
        sum_pN <- sum(file_data[which(is.finite(file_data$pN)),]$pN, na.rm=T)
        sum_pS <- sum(file_data[which(is.finite(file_data$pS)),]$pS, na.rm=T) 
      } else if(dim_method){
        # remove NA and NaN and Inf
        file_data <- file_data[complete.cases(file_data), , drop=FALSE]
        if(comparisons){
          # calculate output values
          temp_count_n <- sum(file_data$all_nonSyn,na.rm=T) + sum(file_data$tip_nonSyn, na.rm=T) 
          temp_count_s <- sum(file_data$all_syn, na.rm=T) + sum(file_data$tip_syn, na.rm=T)
          site_count_n <- sum(file_data$N, na.rm=T)/sum(file_data$comparisons)
          site_count_s <- sum(file_data$S, na.rm=T)/sum(file_data$comparisons)
          sum_pN <- temp_count_n / site_count_n
          sum_pS <- temp_count_s / site_count_s
        } else { 
          # calculate output values
          temp_pN <- sum(file_data$all_nonSyn) + sum(file_data$tip_nonSyn)
          temp_pS <- sum(file_data$all_syn) + sum(file_data$tip_syn)
          file_data$pS <- temp_pS / file_data$S
          file_data$pN <- temp_pN / file_data$N
          sum_pN <- sum(file_data[which(is.finite(file_data$pN)),]$pN, na.rm=T)
          sum_pS <- sum(file_data[which(is.finite(file_data$pS)),]$pS, na.rm=T)
        } 
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
        stop("JC calculation for original data not implemented")
      } else {
    #  ratio <- sum(( (-3*log(1-(4*file_data$pN/3)))/4 ) / ( (-3 * log(1-(4*file_data$pS/3)))/4 ), na.rm=T) 
      }
    } else { # no JC
      if(is.finite(sum_pS)){ 
      if (sum_pS > 10){
        ratio <- sum_pN / sum_pS
      } else {
        ratio <- NA
      }
      }
    }
    
    # add row to output matrix
    mat$name[i] <- as.character(fam_name)
    mat$sum_pN[i] <- sum_pN
    mat$sum_pS[i] <- sum_pS
    mat$ratio[i] <- ratio
    i <- i + 1   
  }
  return(mat)
}
