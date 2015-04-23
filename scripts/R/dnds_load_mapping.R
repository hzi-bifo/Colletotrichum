# load Annotation file

#' A DndsAnalysis Function
#'
#' This function loads all relevant annotation data from the server
#' @param  path to protein_names_mapping in csv format
#' @keywords mapping annotation
#' @export
#' @examples mapping <- dnds_load_mapping()
#' dnds_load_mapping()

dnds_load_mapping <- function(protein_names_mapping="DndsAnalysis/annotation/protein_names_mapping_new.csv"){
 pn_mapping <- read.csv2(protein_names_mapping, sep=";")
 cat(paste("load mapping file from", protein_names_mapping, "\n"))
 return(pn_mapping)
}
