"""
Usage: python3 splicing_workflow.py <gtf_file> <regtools_file>
"""

import os
import sys
import re
import pandas as pd
import subprocess
from gtfparse import read_gtf
from pybiomart import Server, Dataset
import ensembl_rest
from time import sleep

def identify_transcripts(gtf_file, regtools_file):
    """Filter GTF tsv file to find tumor junction coordinates from regtools

    Args:
        gtf_file (string): path to gtf file
        regtools_file (string): path to (filtered) regtools output excel file

    Returns:
        junctions.gtf (file): filtered gtf only containing transcripts that correspond to regtools junctions
        transcripts.fa (file): fasta file with coding transcript sequences corresponding to junctions.gtf file
        gtf_transcripts["gene"] (df): list of altered genes
    """
    # convert gtf to pd df
    gtf_data = read_gtf(gtf_file)
    # read in regtools significant junctions as pd df
    junctions = pd.read_excel(regtools_file, sheet_name = "Sheet1", usecols = ['chrom', 'start', 'end', 'variant_info', 'strand', 'anchor', 'genes'])
    total_transcripts = {}
    for row_index,row in junctions.iterrows():
        start_exons = gtf_data.loc[(gtf_data['end'] == row["start"]) & (gtf_data["seqname"] == row["chrom"])]
        stop_exons = gtf_data.loc[(gtf_data['start'] == row["end"]) & (gtf_data["seqname"] == row["chrom"])]
        transcript_dict = {t:row["genes"] for t in list(start_exons['transcript_id']) if t in list(stop_exons['transcript_id'])}
        total_transcripts.update(transcript_dict)
    gtf_transcripts = gtf_data.loc[(gtf_data["transcript_id"].isin(total_transcripts.keys())) & (gtf_data["feature"] == "transcript")]
    gtf_transcripts["gene"] = [total_transcripts[x] for x in gtf_transcripts["transcript_id"]]
    # write subsetted gtf file
    write_file = open(f'{working_dir}/junctions.gtf', "w")
    for item, gene in zip(list(gtf_transcripts["transcript_id"]), list(gtf_transcripts["gene"])):
        for line in open(gtf_file).readlines():
            if re.search(item, line):
                new_line = line.strip() + f' gene_name "{gene}";\n'
                write_file.write(new_line)
    # create BED12 file from junctions.gtf
    subprocess.Popen("gtfToGenePred junctions.gtf test.genePhred", shell=True)
    subprocess.Popen("genePredToBed test.genePhred results.bed", shell=True)
    # create fasta file with transcript (exon-only) sequences
    subprocess.Popen("bedtools getfasta -fi ~/Documents/ref_fasta/GRCh38.d1.vd1.fa -fo transcripts.fa -bed results.bed -split -name", shell=True)
    return gtf_transcripts

def get_ref_proteins(gene_list):
    """Get wild type protein sequences for each gene with an alternative junction

    Args:
        gene_list (list): list of gene symbols

    Returns:
        final_gene_df (df): reference protein df with sequence, length, and ID
    """
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    total_gene_df = pd.DataFrame()
    gene_info = dataset.query(attributes=["external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "chromosome_name", "start_position", "end_position", "strand", "transcript_start", "transcript_end", "transcription_start_site", "transcript_length", "transcript_tsl", "transcript_biotype"])
    for gene in gene_list:
        gene_df = gene_info.loc[gene_info["Gene name"] == gene]
        gene_df["tsl"] = [int("".join(re.findall(r'\d', x.split(' ')[0]))) for x in list(gene_df["Transcript support level (TSL)"])]
        gene_df = gene_df[(gene_df["Transcript type"] == "protein_coding") & (gene_df["tsl"].isin([1,2]))]
        gene_df["protein sequence"] = [ensembl_rest.sequence_id(x)["seq"] for x in list(gene_df["Protein stable ID"])]
        gene_df["protein length"] = [len(x) for x in list(gene_df["protein sequence"])]
        final_gene_df = gene_df[["Protein stable ID", "protein sequence", "protein length"]]
        final_gene_df["gene"] = gene
        total_gene_df = total_gene_df.append(final_gene_df, ignore_index=True)
    total_gene_df.to_csv("protein_sequences.tsv", sep='\t', index=False)
    return total_gene_df

def get_tumor_proteins(gtf_file, fasta_file):
    subprocess.Popen(f'Rscript transcript_sequences.R {gtf_file} {fasta_file}', shell=True)

def match_sequences(ref_protein_df, tumor_sequences_file):
    """Find sequences that match in reference and tumor plus alteration for neoantigen prediction

    Args:
        final_gene_df (df): reference protein df from get_ref_proteins()
        tumor_sequences_file (string): output file from Rscript command in get_tumor_proteins()

    Returns:
        match_df (df): df with full tumor sequences (ref plus alteration)  
    """
    # load dfs
    reference = ref_protein_df
    tumor = pd.read_table("tumor_sequences.tsv")
    # loop over sequences from each transcript per gene
    transcripts = tumor["tumor_id"].unique().tolist()
    # includes entire seq (1), mismatched portion (2) + fasta junction expansion (for all lengths) (2)
    completed_seq_dict = {}
    for trans in transcripts:
        # transcript-specific tumor sequences
        tumor_seqs = tumor.loc[tumor["tumor_id"] == trans]
        # only ref sequences that match gene from transcript
        ref_seqs = reference.loc[reference["gene"] == "".join(tumor_seqs["gene"].unique().tolist())]
        total_seq_counter = 0
        total_seq_dict = {}
        for p_index, p_row in ref_seqs.iterrows():
            p = ref_seqs.loc[p_index, "protein sequence"]
            for t_index, t_row in tumor_seqs.iterrows():
                t = tumor.loc[t_index, "peptide"] #tumor_seqs.loc[t_index, "peptide"]
                match_seq = ''
                match_seq_counter = 1
                match_seq_dict = {}
                mismatch_seq = ''
                mismatch_seq_counter = 1        
                mismatch_seq_dict = {}
                total_seq = ''
                if p[0:5] == t[0:5]:
                    #print("initial match")
                    for i, (ref,alt) in enumerate(zip(p,t)):
                        if i <= 5:
                            # initialize match_seq and add to match (since I already checked that 0:5 matched)
                            match_seq += alt
                        elif i > 5:
                            if ref == alt:
                                if p[i-1] == t[i-1]:
                                    match_seq += alt
                                if p[i-1] != t[i-1]:
                                    # initialize new match_seq
                                    match_seq = ''
                                    match_seq += alt
                                    if mismatch_seq:
                                        # since now mismatch, add previous mismatch_seq to dictionary
                                        mismatch_seq_dict[f'mismatch{mismatch_seq_counter}'] = mismatch_seq
                                        total_seq += mismatch_seq
                                        mismatch_seq_counter += 1
                            elif ref != alt:
                                if p[i-1] == t[i-1]:
                                    # initialize mismatch_seq
                                    mismatch_seq = ''
                                    mismatch_seq += alt
                                    # since now mismatch, add previous match_seq to dictionary
                                    match_seq_dict[f'match{match_seq_counter}'] = match_seq
                                    total_seq += match_seq
                                    match_seq_counter += 1
                                if p[i-1] != t[i-1]:
                                    mismatch_seq += alt
                                    if alt == '*':
                                        mismatch_seq_dict[f'mismatch{mismatch_seq_counter}'] = mismatch_seq
                                        total_seq += mismatch_seq
                    #print(total_seq)
                    #fasta_expansion = {"classI": }
                    total_seq_counter += 1
                    total_seq_dict[f'{trans}_{total_seq_counter}'] = {"total_seq": total_seq, "mismatch": mismatch_seq}
                    print(total_seq_dict)
                    completed_seq_dict.update(total_seq_dict)
                    total_seq = ''
                    total_seq_counter = 0
                    match_seq = ''
                    match_seq_dict = {}
                    match_seq_counter = 1
                    mismatch_seq = ''
                    mismatch_seq_dict = {}
                    mismatch_seq_counter = 1
    return completed_seq_dict
                
# get working dir
working_dir = os.getcwd()

# input files
gtf_file = '/Users/mrichters/Desktop/Alt_Splicing/runthrough_pipeline/BrMET019-1.denovo.transcripts.gtf'
regtools_file = '/Users/mrichters/Desktop/Alt_Splicing/runthrough_pipeline/all_splicing_variants_default_format.bed_out.xlsx'

junctions_gtf = identify_transcripts(gtf_file, regtools_file)

ref_protein = get_ref_proteins(junctions_gtf["gene"].unique().tolist())

get_tumor_proteins("junctions.gtf", "transcripts.fa")

match_dict = match_sequences(ref_protein, "tumor_sequences.tsv")


