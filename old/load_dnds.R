# This function allows you to import dnds results

load_dnds <- function(in_folder = "dnds/",
                      filename_ending = 18,
                      in_colnames = c("pos", "Sd","Nd", "Sd_tip", "Nd_tip", "S", "N", "comparisons"),
                      gap_folder = "gap/",
                      gap_ending = ".gap.txt",
                      threshold = 1,
                      sum_threshold = 10,
                      quite= FALSE){
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
  for(file in file_list)
  {
    # exctract family identifier
    fam_name <- substr(as.character(file), 1, nchar(as.character(file))-filename_ending)
    if(!quite){
      cat(paste("Load",fam_name,"\n"))
    }
    file_data <- read.table(paste(in_folder,file,sep=""), header=T)
    if (threshold < 1){
      # exclude dnds values based on gap proportion
      if (file.exists(paste(gap_folder,fam_name,gap_ending,sep=""))){
        gap_data <- read.table(paste(gap_folder,fam_name,gap_ending,sep=""), header=F)        
      } else { gap_data <- rep("0",nrow(file_data))}
      file_data <- file_data[which(gap_data<threshold),]
    }
    
    # minimal pathway method with correction for multiple substitutions
    # take mean value for Sd and S based on number of comparisons
    file_data$pS <- (file_data$Sd + file_data$Sd_tip) / file_data$S
    file_data$pN <- (file_data$Nd + file_data$Nd_tip) / file_data$N
    
    # sum without NA and without Inf
    sum_pN <- sum(file_data[which(is.finite(file_data$pN)),]$pN, na.rm=T)
    sum_pS <- sum(file_data[which(is.finite(file_data$pS)),]$pS, na.rm=T)     
    
    if(is.finite(sum_pS)){ 
      if (sum_pS > sum_threshold){
        ratio <- sum_pN / sum_pS
      } else {
        ratio <- NA
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