# load Annotation file

#' A DndsAnalysis Function
#'
#' This function loads annotation data from the server
#' @param  path to annotation file in csv format
#' @keywords mapping annotation
#' @export
#' @examples annotation <- dnds_load_annotation()
#' dnds_load_annotation()

dnds_load_annotation <- function(annotation_file="DndsAnalysis/annotation/Ctv4_flt_genes_annot_new.csv"){
  annotation <- read.csv2(annotation_file, header=T, sep=";")
  cat(paste("load annotation from", annotation_file, "\n"))
  return(annotation)
}
