#' A DndsAnalysis Function
#'
#' This function run the annotation step. It will annotate the dnds table based on the gfam
#' @param dnds table created with dnds_load_files() function. On default it will create a example table
#' @keywords annotation
#' @export
#' @examples 
#' dnds_load_mapping()
run_annotation <- function(dnds_table){
  
  # check parameters
  if(missing(dnds_table)){
    cat("No dnds table is specified. Plese run dnds_load_files() to create a dnds table. \n")
    cat("Load some test table... \n")
    dnds_table <- dnds_load_files(quietly=T) # sometimes you have to run this twice
    dnds_table <- dnds_load_files(quietly=T) # ...
  }
  
  # load relevant files 
 if(!exists("new_old_mapping")){new_old_mapping <- dnds_load_mapping()}
 if(!exists("annotation")){annotation <- dnds_load_annotation()}
 if(!exists("mapping_peg")){mapping_peg <-  read.csv2("~/Desktop/temp/peg_table.csv",sep=";", header=F)  }
  
  # assign gfam to peg using peg_table
  peg_col <- sapply(dnds_table$name, function(x) {
    peg <-gfam2peg(as.character(as.matrix(x)), mapping_peg)
    return(peg)
    })
  dnds_table$peg <- peg_col 
  
  # assign old prot_ids to new_prot_ids using mapping talbe
  old_id_col <- sapply(dnds_table$peg, function(x) {
    old <- new2old(as.character(as.matrix(x)), new_old_mapping)
    return(old)
  })
  dnds_table$old <- old_id_col 
  
  # assign description and effector information 
  annot_desc <- sapply(dnds_table$old, function(x) {
    old <- getAnnotation(as.character(as.matrix(x)), "Description" ,annotation)
    return(old)
  })
  dnds_table$desc <- annot_desc
  # annotate effector
  annot_effector <- sapply(dnds_table$old, function(x) {
    old <- getAnnotation(as.character(as.matrix(x)), "Cg_Eff" ,annotation)
    return(old)
  })
  dnds_table$effector <- annot_effector
  # annotate secreted
  annot_secreted <- sapply(dnds_table$old, function(x) {
    old <- getAnnotation(as.character(as.matrix(x)), "secreted" ,annotation)
    return(old)
  })
  dnds_table$secreted <- annot_secreted
  # annotate B2GO_BestHitSpecies CSEP
  annot_B2GO_species <- sapply(dnds_table$old, function(x) {
    old <- getAnnotation(as.character(as.matrix(x)), "Cg_bestHit" ,annotation)
    return(old)
  })
  dnds_table$Cg_bestHit <- annot_B2GO_species

 annot_CSEP <- sapply(dnds_table$old, function(x) {
   old <- getAnnotation(as.character(as.matrix(x)), "CSEP" ,annotation)
   return(old)
 })
 dnds_table$CSEP <- annot_CSEP
 
 return(dnds_table) 
}

