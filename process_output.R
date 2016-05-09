output_folder="results/"
dnds_folder="data/dnds/"
gap_folder="data/gap/"

source("dnds_load_files.R")

# process gfam output of dnds pipeline
pooled <- dnds_load_files(in_folder =dnds_folder,
                                    gap_folder = gap_folder,
                                    in_colnames = c("all_syn", "all_nonSyn", "tip_syn", "tip_nonSyn", "S", "N", "comparisons"),
                                    dim_method=T,
                                    quietly = TRUE,
                                    comparisons=T,
                                    threshold = 0.7)

write.table(pooled, file="pooled_0.7.txt",row.names=F, col.names=F, quote=F, sep=";")