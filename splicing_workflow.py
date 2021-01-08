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

def identify_transcripts(gtf_file, regtools_file, sample_name):
    """Filter GTF tsv file to find tumor junction coordinates from regtools

    Args:
        gtf_file (string): path to gtf file
        regtools_file (string): path to (filtered) regtools output excel file

    Returns:
        junctions.gtf (file): filtered gtf only containing transcripts that correspond to regtools junctions
        transcripts.fa (file): fasta file with coding transcript sequences corresponding to junctions.gtf file
        gtf_transcripts (df): list of altered transcripts
    """
    # convert gtf to pd df
    gtf_data = read_gtf(gtf_file)
    # read in regtools significant junctions as pd df
    junctions = pd.read_excel(regtools_file, sheet_name='Sheet1')
    junctions = junctions.loc[junctions['Sample'] == sample_name]
    print(junctions)
    total_transcripts = {}
    for row_index,row in junctions.iterrows():
        start_exons = gtf_data.loc[(gtf_data['end'] == row["start"]) & (gtf_data["seqname"] == row["chrom"])]
        stop_exons = gtf_data.loc[(gtf_data['start'] == row["end"]) & (gtf_data["seqname"] == row["chrom"])]
        print(start_exons, stop_exons)
        transcript_dict = {t:row["gene_names"] for t in list(start_exons['transcript_id']) if t in list(stop_exons['transcript_id'])}
        total_transcripts.update(transcript_dict)
    gtf_transcripts = gtf_data.loc[(gtf_data["transcript_id"].isin(total_transcripts.keys())) & (gtf_data["feature"] == "transcript")]
    gtf_transcripts["gene"] = [total_transcripts[x] for x in gtf_transcripts["transcript_id"]]
    # filter gtf_transcripts
    if "transcript_type" in gtf_transcripts.columns:
        gtf_transcripts = gtf_transcripts.loc[gtf_transcripts["transcript_type"] == "protein_coding"]
    # write subsetted gtf file
    write_file = open("junctions.gtf", "w")
    for item, gene in zip(list(gtf_transcripts["transcript_id"]), list(gtf_transcripts["gene"])):
        for line in open(gtf_file).readlines():
            if re.search(item, line):
                new_line = line.strip() + f' gene_name "{gene}";\n'
                write_file.write(new_line)
    write_file.close()
    # create BED12 file from junctions.gtf
    # create fasta file with transcript (exon-only) sequences
    subprocess.Popen("gtfToGenePred junctions.gtf test.genePhred && genePredToBed test.genePhred results.bed && bedtools getfasta -fi ~/Documents/ref_fasta/GRCh38.d1.vd1.fa -fo transcripts.fa -bed results.bed -split -name -s && sed -i.bak 's/(-)//;s/(+)//' transcripts.fa", shell=True)
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
        # filter out NaN values
        gene_df = gene_df.loc[gene_df["Transcript support level (TSL)"].astype(str).str.contains("tsl")]
        gene_df["tsl"] = [re.search(r'\d', x).group() for x in list(gene_df["Transcript support level (TSL)"])]
        # filter by protein coding and TSL == 1,2
        gene_df = gene_df[(gene_df["Transcript type"] == "protein_coding") & (gene_df["tsl"].isin(["1","2"]))]
        gene_df["protein sequence"] = [ensembl_rest.sequence_id(x)["seq"] for x in list(gene_df["Protein stable ID"])]
        gene_df["protein length"] = [len(x) for x in list(gene_df["protein sequence"])]
        final_gene_df = gene_df[["Protein stable ID", "protein sequence", "protein length"]]
        final_gene_df["gene"] = gene
        total_gene_df = total_gene_df.append(final_gene_df, ignore_index=True)
    total_gene_df.to_csv("protein_sequences.tsv", sep='\t', index=False)
    return total_gene_df

def get_tumor_proteins(gtf_file, fasta_file):
    """Run Rscript that uses biomaRt to find tumor protein sequences from junctions.gtf

    Args:
        gtf_file (file): junctions.gtf created in identify_transcripts()
        fasta_file (file): transcripts.fa created in identify_transcripts()

    Returns:
        tumor_sequences.tsv (file): tsv file with tumor sequence, sequence length, and relevant transcript/gene
    """
    subprocess.Popen(f'Rscript {original_dir}/transcript_sequences.R {gtf_file} {fasta_file}', shell=True)
    print(f'Rscript {original_dir}/transcript_sequences.R {gtf_file} {fasta_file}')

def match_sequences(ref_protein_df, tumor_sequences_file):
    """Find sequences that match in reference and tumor plus alteration for neoantigen prediction

    Args:
        final_gene_df (df): reference protein df from get_ref_proteins()
        tumor_sequences_file (string): output file from Rscript command in get_tumor_proteins()

    Returns:
        match_df (df): df with full tumor sequences (ref plus alteration)  
    """
    # load dfs
    reference = ref_protein
    tumor = pd.read_table("tumor_sequences.tsv")
    # loop over sequences from each transcript per gene
    transcripts = tumor["tumor_id"].unique().tolist()
    # includes entire seq (1), mismatched portion (2) + fasta junction expansion (for all lengths) (2)
    completed_seq_dict = {}
    for trans in transcripts:
        # transcript-specific tumor sequences
        tumor_seqs = tumor.loc[tumor["tumor_id"] == trans]
        # only ref sequences that match gene from transcript
        # catching the error of more than one gene in tumor_seqs - impossible though?
        ref_seqs = reference.loc[reference["gene"] == "".join(tumor_seqs["gene"].unique().tolist())]
        for p_index, p_row in ref_seqs.iterrows():
            p = ref_seqs.loc[p_index, "protein sequence"]
            for t_index, t_row in tumor_seqs.iterrows():
                t = tumor.loc[t_index, "peptide"]
                #match_seq = ''; match_seq_counter = 1; match_seq_dict = {}
                #mismatch_seq = ''; mismatch_seq_counter = 1; mismatch_seq_dict = {}
                #total_seq = ''
                total_dict = {}
                for i, (ref,alt) in enumerate(zip(p,t)):
                    # IF BASES MATCH
                    # add to first match seq - if ref == alt and ref-1 == alt-1
                    if ref == alt:
                        print('yes')
                        # if first base, add match to match seq
                        if i == 0:
                            match_seq += alt               
                        # if NOT first base
                        elif i > 0:
                            # this base and previous base must match to add to match seq
                            if p[i-1] == t[i-1]:
                                match_seq += alt
                    if ref != alt:
                        print('no')
                        # WHY? if at base right before stop codon and whole sequence is a complete match, add to total_seq
                        #elif i == len(p)-1 and t[i] == '*':
                        #        total_seq += match_seq
                            # if this base matches but previous base does not, assuming that this is SECOND match seq beginning
                            # if p[i-1] != t[i-1]:
                            #     # initialize new match_seq
                            #     match_seq = ''
                            #     match_seq += alt
                            #     # add mismatch seq back to total seq if this is not the first match seq
                            #     if mismatch_seq:
                            #         # since now mismatch, add previous mismatch_seq to dictionary
                            #         total_dict["end"] = i
                            #         total_dict["mismatch_sequence"] = mismatch_seq
                            #         total_seq += mismatch_seq
                            #         mismatch_seq_counter += 1
                    # IF BASES DON'T MATCH
                    # starting a mismatch seq
                    elif ref != alt:
                        if p[i-1] == t[i-1]:
                            # initialize mismatch_seq
                            mismatch_seq = ''
                            mismatch_seq += alt
                            total_dict["start"] = i
                            # since now mismatch, add previous match_seq to dictionary
                            match_seq_dict[f'match{match_seq_counter}'] = match_seq
                            total_seq += match_seq
                            match_seq_counter += 1
                        # continue to add to mismatch seq
                        if p[i-1] != t[i-1]:
                            mismatch_seq += alt
                            # if mismatch ends with frameshift, add mismatch seq to total_seq
                            if alt == '*':
                                total_dict["end"] = i
                                total_dict["mismatch_sequence"] = mismatch_seq
                                total_seq += mismatch_seq
            if total_dict:
                total_dict["transcript_id"] = [trans]
            print(total_dict)
            if total_seq in completed_seq_dict.keys():
                completed_seq_dict[total_seq]["transcript_id"].append(trans)
            else:    
                completed_seq_dict[total_seq] = total_dict
            total_seq = ''
            match_seq = ''; match_seq_dict = {}; match_seq_counter = 1
            mismatch_seq = ''; mismatch_dict = {}; mismatch_seq_counter = 1
    final_dict = {"total_sequence": []}
    for k in completed_seq_dict.keys():
        final_dict["total_sequence"].append(k)
        value_dict = completed_seq_dict[k]
        for key in value_dict.keys():
            if key in final_dict.keys():
                final_dict[key].append(value_dict[key])
            else: 
                final_dict[key] = [value_dict[key]]
    final_dict["transcript_id"] = [",".join(item) if len(item) > 1 and isinstance(item, list) else str(item[0]) for item in final_dict["transcript_id"]]
    final_df = pd.DataFrame(final_dict)  
    return final_df

def create_fasta(match_df, epitope_lengths=[8,9,10,11,15]):
    """Create fasta file with altered transcript sequences for input into pVACbind

    Args:
        match_df (df): output from match_sequences that contains the total transcript sequences and mismatch portion that will be used to create fasta sequence
        epitope_lengths (list): desired epitope lengths for pVACbind

    Returns:
        pvacbind_sequences.fa (file): fasta file for input into pVACbind 
    """
    write_file = open("pvacbind_sequences.fa", "w")
    for row_index, row in match_df.iterrows():
        fasta_seqs = [match_df.loc[row_index, "total_sequence"][match_df.loc[row_index, "start"] - l:match_df.loc[row_index, "end"]] for l in epitope_lengths]
        fasta_headers = [f'>{match_df.loc[row_index, "transcript_id"]} length{l}' for l in epitope_lengths]
        for x,y in zip(fasta_seqs, fasta_headers):
            new_line = f'{y}\n{x}\n\n'
            ### create df and filter out redundant seqs - then write file ###
            write_file.write(new_line)
    write_file.close()

# get working dir
original_dir = os.getcwd()
working_dir = '/Users/mrichters/Desktop/Alt_Splicing/new_samples/'
results_dir = working_dir + 'results'
os.chdir(working_dir)

# input files
gtf_file = working_dir + 'stringtie/BrMET019-2.stringtie_transcripts.gtf'
regtools_file = working_dir + 'brmet/filtered_regtools/BrMET019_filtered_junctions.xlsx'
#gtf_file = working_dir + 'gencode.v22.annotation.gtf'
#regtools_file = working_dir + 'RNF145_sampleline.tsv'

junctions_gtf = identify_transcripts(gtf_file, regtools_file, "BrMET019-2")

ref_protein = get_ref_proteins(junctions_gtf["gene"].unique().tolist())

get_tumor_proteins("junctions.gtf", "transcripts.fa")

match_df = match_sequences(ref_protein, "tumor_sequences.tsv")

create_fasta(match_df)

# read in ref and tumor peptide seqs
reference = ref_protein
tumor = pd.read_table("tumor_sequences.tsv")
# loop over sequences from each transcript per gene
transcripts = tumor["tumor_id"].unique().tolist()

for trans in transcripts:
    # transcript-specific tumor sequences
    tumor_seqs = tumor.loc[tumor["tumor_id"] == trans]
    # only ref sequences that match gene from transcript
    # catching the error of more than one gene in tumor_seqs - impossible though?
    ref_seqs = reference.loc[reference["gene"] == "".join(tumor_seqs["gene"].unique().tolist())]
    # iterating over ref protein sequences (1 whole sequence at a time)
    for p_index, p_row in ref_seqs.iterrows():
        p = ref_seqs.loc[p_index, "protein sequence"]
        # iterating over tumor peptide sequences (1 whole sequence at a time)
        for t_index, t_row in tumor_seqs.iterrows():
            t = tumor.loc[t_index, "peptide"]
            print(p_index, t_index)
            # iterating over each p for each t
            # iterating over each base of p for each t
            for index in range(len(p)):
                #print(index)
                match_seq = ''; match_seq_counter = 0; match_seq_dict = {}
                mismatch_seq = ''; mismatch_seq_counter = 0; mismatch_seq_dict = {}
                total_seq = ''; total_dict = {}
                if p[0:5] == t[0:5]:
                    first_match = True
                    print("initial match", p[0:5])
                else:
                    first_match = False
                for i, (ref,alt) in enumerate(zip(p[index:len(p)],t)):
                    # starting a match seq
                    if i <= 5 and first_match == True:
                        # initialize match_seq and add to match (since I already checked that 0:5 matched)
                        match_seq += alt
                    elif i > 5:
                        if ref == alt:
                            # add to first match seq
                            if p[i-1] == t[i-1]:
                                match_seq += alt
                            # start new match seq
                            if p[i-1] != t[i-1]:
                                # initialize new match_seq
                                match_seq_counter += 1
                                match_seq_dict[f'match{match_seq_counter}'] = match_seq
                                match_seq = ''
                                match_seq += alt
                                # add mismatch seq back to total seq if this is not the first match seq
                                # if mismatch_seq:
                                    # since now mismatch, add previous mismatch_seq to dictionary
                                    # total_dict["end"] = i
                                    # total_dict["mismatch_sequence"] = mismatch_seq
                                    # total_seq += mismatch_seq
                                    # mismatch_seq_counter += 1
                    