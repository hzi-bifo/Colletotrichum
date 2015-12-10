options(warn=-1)

dnds_folder="data/dnds/"
gap_folder="data/gap/"

source("load_dnds.R")
# sed -i -e 's/all_syn\tall_nonSyn\ttip_syn\ttip_nonSyn/pos\tSd\tNd\tSd_tip\tNd_tip/g' dnds/gfam_*


# process gfam output of dnds pipeline
pooled <- load_dnds(in_folder = dnds_folder,
                                gap_folder = gap_folder,
                                threshold = 0.7,
                                filename_ending = 14)

write.table(pooled, file="pooled_0.7.txt",row.names=F, col.names=F, quote=F, sep=";")

# statistical test

pooled_pval <- fisher_test(pooled3)

write.table(pooled, file="pooled_0.7_pval.txt",row.names=F, col.names=F, quote=F, sep=";")

