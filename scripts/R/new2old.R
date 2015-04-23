#' A DndsAnalysis Function
#'
#' This function assign a old prot_id to a new prot_id
#' @param  
#' @keywords annotation
#' @export
#' @examples new2old("0001.peg.1") 
#' new2old()
new2old <- function(new_id, mapping){
  # check arguments
  if(!is.character(new_id))
    stop("Error: new_id must be a character")
  if(missing(mapping)){
    cat("No mapping data specified, Read default mapping matrix. This will take some time. \n")
    mapping <- dnds_load_mapping()    
  }
  #cat(paste("translate",new_id,"...\n"))

  # lookup gfam in mapping file
  old_id <- mapping$old_prot_ids[which(mapping$new_prot_ids == new_id)]

  return(as.character(as.matrix(old_id)))
}
