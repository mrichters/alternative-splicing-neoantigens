# 09.13.20

library(tidyverse)
library(WriteXLS)

combine_junctions <- function(cohort, tumor) {
    tsv_files <- list.files(cohort)
    tumor_files <- tsv_files[grepl(tumor, tsv_files)]
    df <- c()
    for ( f in tumor_files ) {
        name <- gsub(".cse_identify_filtered_compare_default.tsv", "", f)
        t <- read_tsv(paste(cohort, f, sep='/'), col_names=T)
        df <- rbind(df, data.frame("Sample" = name, t))
    }
    j <- filter(df, score >= 50 & anchor %in% c("D", "A", "NDA"))
    j <- j[order(j$variant_info),]
    WriteXLS(j, paste0("filtered_regtools/", tumor, "_filtered_junctions.xlsx"))
    #return(j)
}

combine_junctions_samples <- function(cohort, sample) {
    #tsv_files <- list.files(cohort)
    #tumor_files <- tsv_files[grepl(sample, tsv_files)]
    t <- read_tsv(paste(cohort, sample, sep='/'), col_names=T)
    df <- rbind(df, data.frame("Sample" = sample, t))
    j <- filter(df, score >= 50 & anchor %in% c("D", "A", "NDA"))
    j <- j[order(j$variant_info),]
    WriteXLS(j, paste0("filtered_regtools_sample/", tumor, "_filtered_junctions.xlsx"))
    #return(j)
}

setwd("~/Desktop/Alt_Splicing/new_samples/")

tsv_files <- list.files("gbm/regtools/")

#tumors <- unique(sapply(tsv_files, function(x) unlist(strsplit(gsub(".cse_identify_filtered_compare_default.tsv", "", x), '-'))[1]))
samples <- sapply(tsv_files, function(x) gsub(".cse_identify_filtered_compare_default.tsv", "", x))

for (t in tsv_files[1]) {
    print(t)
    combine_junctions_samples("gbm/regtools", t)
}
