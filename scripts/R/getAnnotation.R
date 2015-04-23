#' A DndsAnalysis Function
#'
#' This function returns a annotation based on the old_id. 
#' @param old_id
#' @param column for annotation e.g. "Description" to get the Description colum or "secreted" to get the secreted column
#' @param the annotation matrix. If not specified it will load a matrix.
#' @keywords annotation
#' @export
#' @examples getAnnotation("CT04_00002", "Eff_NGP") 
#' getAnnotation()
getAnnotation <- function(old_id="CT04_00002", column="Description", annotation){
  # check arguments
  if(!is.character(new_id))
    stop("Error: old_id must be a character")
  if(missing(annotation)){
    cat("No annotation data specified, Read default mapping matrix. This will take some time. \n")
    annotation <- dnds_load_annotation()    
  }
 # cat(paste("annotate",column,"of",old_id,"...\n"))
  
  # lookup gfam in mapping file
  row <- which(annotation$GeneID == old_id)
  annot <- annotation[column][row,]
  return(as.character(as.matrix(annot)))
}
