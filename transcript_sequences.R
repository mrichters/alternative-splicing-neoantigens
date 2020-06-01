# load packages
library(polyester)
library(ORFik)
library(tidyverse)
library(biomaRt)
library(BSgenome)
library(stringdist)
library(zipR)
library(seqinr)
library(Biostrings)

# setwd
setwd("~/Desktop/Alt_Splicing/runthrough_pipeline/pipeline_results/")

# command line arguments
args = commandArgs(trailingOnly = T)
gtf_junctions <- "junctions.gtf" #args[1]
fasta <- "transcripts.fa" #args[2]

get_tumor_sequences <- function(fasta, junction_name, gene_name) {
    seqs <- readDNAStringSet(fasta, "fasta")
    junction <- seqs[[junction_name]]
    transcript_seq <- DNAStringSet(junction)
    names(transcript_seq) <- junction_name
    # translated into all 3 reading frames and all ORFs found (in total)
    pos <- findORFs(transcript_seq, startCodon = startDefinition(6), longestORF = FALSE, minimumLength = 0)
    orf_seqs <- getSeq(seqs, GRanges(junction_name, unlist(pos), strand = '+'))
    # create aa sequences df
    orf_seqs_aa <- as.data.frame(Biostrings::translate(orf_seqs))
    colnames(orf_seqs_aa) <- "peptide"
    orf_seqs_aa$peptide_length <- apply(orf_seqs_aa, 1, FUN = function(x) nchar(x[1]))
    orf_seqs_aa <- filter(orf_seqs_aa, peptide_length > 20)
    orf_seqs_aa$tumor_id <- junction_name
    orf_seqs_aa$gene <- gene_name
    orf_seqs_aa <- orf_seqs_aa[rev(order(orf_seqs_aa$peptide_length)), ]
    return(orf_seqs_aa)
}

# read in input junctions file
junctions <- as.data.frame(read_tsv(gtf_junctions, col_names = c('seqname', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'attributes')))
# subset with only transcript lines
junction_transcripts <- filter(junctions, feature == 'transcript')
final_tumor_df <- c()
# loop over transcripts 
for ( row in 1:nrow(junction_transcripts) ) {
    # get strand info - shouldn't be needed - strand was specified in getfasta
    stranded <- junction_transcripts[row, "strand"]
    # list all attributes
    attributes <- unlist(strsplit(str_replace_all(junction_transcripts[row, "attributes"], fixed(" "), ""), ";"))
    # get gene name from attributes
    gene_name <- unlist(strsplit(attributes[grep("gene_name", attributes)], '"'))[2]
    # get transcript name from attributes
    junction_name <- unlist(strsplit(attributes[grep("transcript_id", attributes)], '"'))[2]
    print(c(stranded, gene_name, junction_name))
    # tumor df
    tumor_df <- get_tumor_sequences(fasta, junction_name, gene_name)
    final_tumor_df <- rbind(final_tumor_df, tumor_df)
}
write_tsv(final_tumor_df, "tumor_sequences.tsv")