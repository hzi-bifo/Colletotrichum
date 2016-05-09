# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

annotate_csep <- function(result_mat){
  
  # load list of CSEP gene families
  csep <- read.table("annotation/csep_list.txt", header=F)
  
  # load SSP gene families
  ssp <- read.table("annotation/ssp_list.txt", header=F)
  
  # match gene fam ids
  csep_match <- match(as.character(as.matrix(csep)), as.character(as.matrix(result_mat$name)))
  
  # match gene fam ids
  ssp_match <- match(as.character(as.matrix(ssp)), as.character(as.matrix(result_mat$name)))
  
  # subset result_mat csep
  result_csep <- result_mat[csep_match[!is.na(csep_match)],]
  
  # subset result_mat non-csep
  result_non_csep <- result_mat[-csep_match[!is.na(csep_match)],]
  
  # subset result_mat non-csep non ssp
  result_non_csep <- result_mat[-csep_match[!is.na(csep_match)],]
  result_non_csep_ssp <- result_non_csep[-ssp_match[!is.na(ssp_match)],]
  
  # subset result_mat ssp
  result_ssp<- result_mat[ssp_match[!is.na(ssp_match)],]
  
  # mark ssp in matrix
  result_mat$type <- "All"
  result_mat[ssp_match[!is.na(ssp_match)],]$type <- "Secreted"
  
  # mark csep in matrix
  result_mat[csep_match[!is.na(csep_match)],]$type <- "CSEPs"
  
  # detailed annotation
  result_mat$is.secreted <- FALSE
  result_mat$is.csep <- FALSE
  result_mat[which(!is.na(match(result_mat$name, csep$V1))),]$is.csep <- "TRUE"
  result_mat[which(!is.na(match(result_mat$name, ssp$V1))),]$is.secreted <- "TRUE"
  return(result_mat)
}