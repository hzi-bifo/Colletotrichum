#' A DndsAnalysis Function
#'
#' This function assign a list of peg id to given gfam_id. To reduce runtime please supply mapping file as parameter
#' @param  
#' @keywords annotation
#' @export
#' @examples 
#' gfam2peg()
gfam2peg <- function(gfam,
                     mapping,
                     species_number = "1001"){
  # check arguments
  if(!is.character(gfam))
    stop("Error: gfam must be a character")
  if(missing(mapping)){
    cat("No mapping data specified, Read default mapping matrix. This will take some time. \n")
    mapping <- read.csv2("DndsAnalysis/annotation/peg_table.csv",
                        sep=";", header=F)    
  }
  #cat(paste("translate",gfam,"...\n"))
  # lookup gfam in mapping file
  row_match <- mapping[which(mapping[,1] == gfam),]
  row_match <- row_match[which(nchar(as.matrix(row_match))>0)] # get non-empty ids
  
  # if there is a annotation on 1001. then take a random 1001 annotation. If not take a random annotation
  peg_1001 <- grep(paste("^",species_number, sep=""),as.character(as.matrix(row_match)))
  if (length(peg_1001)>0){
    # take a random one
    peg <- row_match[peg_1001[sample(1:length(peg_1001),1)]]
  } else
    # take random one
    peg <- row_match[sample(2:length(row_match),1)]
  
  return(as.character(as.matrix(peg)))
}
