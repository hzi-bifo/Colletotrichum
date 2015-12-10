# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

annotate_csep <- function(result_mat){
  
  # load list of CSEP gene families
  csep <- read.table("annotation/gfam_csep.txt", header=F)
  
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
  result_mat$type <- "gene families"
  result_mat[ssp_match[!is.na(ssp_match)],]$type <- "SSP"
  
  # mark csep in matrix
  result_mat[csep_match[!is.na(csep_match)],]$type <- "CSEP"
  
  # fisher test  CSEP vs rest
  testmat <- matrix(c(sum(result_csep$sum_pN, na.rm=T), sum(result_non_csep$sum_pN, na.rm=T), sum(result_csep$sum_pS, na.rm=T), sum(result_non_csep$sum_pS, na.rm=T)),
                    nrow = 2, dimnames = list(c("CSEP", "non-SSP-CSEP"),c("sumDn", "sumDs")))
  fisher.test(testmat, alternative = "greater") # p-value < 2.2e-16, OR=2.128149
  
  # fisher test  SSP vs rest
  testmat <- matrix(c(sum(result_ssp$sum_pN, na.rm=T), sum(result_non_csep_ssp$sum_pN, na.rm=T), sum(result_ssp$sum_pS, na.rm=T), sum(result_non_csep_ssp$sum_pS, na.rm=T)),
                    nrow = 2, dimnames = list(c("SSP", "non-SSP+CSEP"),c("sumDn", "sumDs")))
  fisher.test(testmat, alternative = "greater") # p-value < 2.2e-16, OR=2.128149
  
  # add annotation
  suptable$gfam <- "NA"
  suptable <- read.table("annotation/csep.txt",se=";", header=T)
  annottable <- read.table("annotation/gfam2id.txt", sep=";")
  supmatch <- match(as.character(as.matrix(suptable$GeneID)), as.character(as.matrix(annottable$V2)))
  suptable$gfam<- annottable[supmatch,]$V1
  supmatch2 <- match(as.character(as.matrix(suptable$gfam)), as.character(as.matrix(result_mat)))
  suptable$dN.dS <- result_mat[supmatch2,]$ratio
  write.table(suptable, file = "annotated_table.txt", quote=F, row.names=F, sep=";")
  return(result_mat)
}