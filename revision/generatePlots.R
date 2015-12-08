# script by Philipp C. Münch
# philipp.muench@helmholtz-hzi.de

require(ggplot2)
require(scales)
source("reverselog_trans.R")

# load list of CSEP gene families
csep <- read.table("data/annotation/gfam_csep.txt", header=F)

# load dnds pipeline output
result_mat <- read.table("results/summary/pooled.gap_theshold_0.7.summary.txt_pval.txt", sep=";", header=T)

# match gene fam ids
csep_match <- match(as.character(as.matrix(csep)), as.character(as.matrix(result_mat$name)))

# subset result_mat csep
result_csep <- result_mat[csep_match[!is.na(csep_match)],]

# subset result_mat non-csep
result_non_csep <- result_mat[-csep_match[!is.na(csep_match)],]

# mark csep in matrix
result_mat$type <- "non-CSEP"
result_mat[csep_match[!is.na(csep_match)],]$type <- "CSEP"


# fisher test
# testmat <- matrix(c(sum(result_csep$sum_pN, na.rm=T), sum(result_non_csep$sum_pN, na.rm=T), sum(result_csep$sum_pS, na.rm=T), sum(result_non_csep$sum_pS, na.rm=T)),
                    nrow = 2, dimnames = list(c("CSEP", "non-CSEP"),c("sumDn", "sumDs")))

fisher.test(testmat, alternative = "greater") # p-value < 2.2e-16, OR=2.128149

# add annotation
suptable$gfam <- "NA"
suptable <- read.table("data/annotation/csep.txt",se=";", header=T)
annottable <- read.table("data/annotation/gfam2id.txt", sep=";")
supmatch <- match(as.character(as.matrix(suptable$GeneID)), as.character(as.matrix(annottable$V2)))
suptable$gfam<- annottable[supmatch,]$V1

supmatch2 <- match(as.character(as.matrix(suptable$gfam)), as.character(as.matrix(result_mat)))
suptable$dN.dS <- result_mat[supmatch2,]$ratio
write.table(suptable, file = "annotated_table.txt", quote=F, row.names=F, sep=";")

# start plot
d <- ggplot(result_mat, aes(x=fdr,y=ratio, colour=type))

# scale und theme
d <- d + scale_x_continuous(trans=reverselog_trans(10), limits=c(1,1*10^(-168)))
d <- d + geom_point(alpha=0.4)  + theme_bw() 


e <- ggplot(result_mat, aes(x=type, y=ratio, colour=type)) + theme_bw() + scale_y_log10()
e <- e + geom_boxplot(notch = TRUE, outlier.colour = "black", outlier.size = 1) #+ geom_jitter(alpha=0.01)

e <- ggplot(result_mat, aes(x=type, y=ratio, colour=type, fill=type)) + theme_bw() + scale_y_log10() + theme_minimal()
e <- e + geom_violin(trim=FALSE, fill="gray90", colour="black") +  geom_boxplot(width=0.1, outlier.colour = "darkred", outlier.size = 1, colour="black", fill="grey30")  #+ geom_jitter(alpha=0.01)