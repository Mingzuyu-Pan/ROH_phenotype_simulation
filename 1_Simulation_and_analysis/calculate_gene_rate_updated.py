#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 10:55:02 2025

@author: mzyp
"""

import pandas as pd

gtf_file = "Homo_sapiens.GRCh37.87.chr.gtf"
gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

# To filter the exon data we want from the raw data
exon_df = gtf_df
exon_df = gtf_df[gtf_df["feature"] == "exon"]

# This three columns are not used in my calculation
exon_df = exon_df.drop(columns=["source", "score", "frame"])


exon_df["gene_id"] = exon_df["attribute"].str.extract(r'gene_id "([^"]+)"')
exon_df["gene_version"] = exon_df["attribute"].str.extract(r'gene_version "([^"]+)"')
exon_df["transcript_id"] = exon_df["attribute"].str.extract(r'transcript_id "([^"]+)"')
exon_df["transcript_version"] = exon_df["attribute"].str.extract(r'transcript_version "([^"]+)"')
exon_df["exon_number"] = exon_df["attribute"].str.extract(r'exon_number "([^"]+)"')
exon_df["gene_name"] = exon_df["attribute"].str.extract(r'gene_name "([^"]+)"')
exon_df["gene_source"] = exon_df["attribute"].str.extract(r'gene_source "([^"]+)"')
exon_df["gene_biotype"] = exon_df["attribute"].str.extract(r'gene_biotype "([^"]+)"')
exon_df["transcript_name"] = exon_df["attribute"].str.extract(r'transcript_name "([^"]+)"')
exon_df["transcript_source"] = exon_df["attribute"].str.extract(r'transcript_source "([^"]+)"')
exon_df["transcript_biotype"] = exon_df["attribute"].str.extract(r'transcript_biotype "([^"]+)"')
exon_df["exon_id"] = exon_df["attribute"].str.extract(r'exon_id "([^"]+)"')
exon_df["exon_version"] = exon_df["attribute"].str.extract(r'exon_version "([^"]+)"')
exon_df["tag"] = exon_df["attribute"].str.extract(r'tag "([^"]+)"')
exon_df = exon_df.drop(columns=["attribute"])

# This step is to make sure the transcripts we chose come from the protein coding gene
# Also to filter sudo gene.
exon_df = exon_df[exon_df["gene_biotype"] == "protein_coding"]
exon_df = exon_df[exon_df["transcript_biotype"] == "protein_coding"]

exon_df_filtered = exon_df.drop(columns=["tag","transcript_biotype","gene_biotype","gene_source","transcript_source"])

# to make sure the calculation is correct and not influenced by the type of data
exon_df_filtered["end"] = exon_df_filtered["end"].astype('int')
exon_df_filtered["start"] = exon_df_filtered["start"].astype('int')
exon_df_filtered["exon_length"] = exon_df_filtered["end"] - exon_df_filtered["start"] + 1

# Same. Prepared for the following analysis.
exon_df_filtered['gene_name'] = exon_df_filtered['gene_name'].astype('string')
exon_df_filtered['transcript_id'] = exon_df_filtered['transcript_id'].astype('string')
exon_df_filtered['transcript_version'] = exon_df_filtered['transcript_version'].astype('string')
exon_df_filtered['seqname'] = exon_df_filtered['seqname'].astype('string')


def merge_intervals_length(group):
    # ascending by default
    sorted_group = group.sort_values("start")
    # Note this is a list
    merged_intervals = []

    # To access each row in a convenient way
    for idx, row in sorted_group.iterrows():
        current_start = row["start"]
        current_end = row["end"]
        if not merged_intervals:
            merged_intervals.append([current_start, current_end])
        else:
            last_interval = merged_intervals[-1]
            if current_start <= last_interval[1]:
                last_interval[1] = max(last_interval[1], current_end)
            else:
                merged_intervals.append([current_start, current_end])
    
    total_length = sum(interval[1] - interval[0] + 1 for interval in merged_intervals)
    
    return pd.Series({
        "exon_length": total_length,
        "start": sorted_group["start"].min(),
        "end": sorted_group["end"].max()
    })


exon_length_summary = exon_df_filtered.groupby(
    ["gene_name", "seqname"]
).apply(merge_intervals_length).reset_index()



longest_exon = exon_length_summary

# We are preparing for following analysis now.
# To make sure we only want the transcripts from autosomal chromosomes.
print(longest_exon["seqname"].unique())
longest_exon["seqname"] = longest_exon["seqname"].astype(str)
print(longest_exon["seqname"].unique())
autosomes_to_keep = [str(i) for i in range(1, 23)]
longest_exon = longest_exon[longest_exon['seqname'].isin(autosomes_to_keep)]
print(longest_exon["seqname"].unique())


# To open the gene lists.
with open('universe_list.txt', 'r') as file:
    universe_list = file.read().splitlines()

with open('dominant_list.txt', 'r') as file:
    dominant_list = file.read().splitlines()
    
with open('recessive_list.txt', 'r') as file:
    recessive_list = file.read().splitlines()    

# Set up the data type to prepare for the following analysis and calculation.
are_all_strings = all(isinstance(gene, str) for gene in dominant_list)
print("Are they all strings?(dominant list)", are_all_strings)
are_all_strings_1 = all(isinstance(gene, str) for gene in recessive_list)
print("Are they all strings?(recessive list)", are_all_strings_1)
are_all_strings_2 = all(isinstance(gene, str) for gene in universe_list)
print("Are they all strings?(universe list)", are_all_strings_2)

print("The data type of 'gene_name' in longest_exon:", longest_exon['gene_name'].dtype)
longest_exon['gene_name'] = longest_exon['gene_name'].astype('string')
print("The data type of 'gene_name' in longest_exon:", longest_exon['gene_name'].dtype)
    
# To check if there are any genes in dominant and recessive but not in universe.
extra_dominant = [gene for gene in dominant_list if gene not in universe_list]
extra_recessive = [gene for gene in recessive_list if gene not in universe_list]

# The answer is yes. Then let's add them.
universe_list.extend(extra_dominant)
universe_list.extend(extra_recessive)

# To make sure there are not repeated genes in each gene list.
universe_list = list(set(universe_list))
dominant_list = list(set(dominant_list))
recessive_list = list(set(recessive_list))

# To check if there is any gene belonging to both dominant and recessive.
intersection = [gene for gene in dominant_list if gene in recessive_list] 
# And let's remove it
universe_list = [gene for gene in universe_list if gene not in intersection]
dominant_list = [gene for gene in dominant_list if gene not in intersection]
recessive_list = [gene for gene in recessive_list if gene not in intersection]
# Then we can define the additive list
additive_list = [gene for gene in universe_list if gene not in dominant_list and gene not in recessive_list]

# To check if the gene name in longest_exon exists in the corresponding gene list.
universe_gene = longest_exon[longest_exon['gene_name'].isin(universe_list)]
dominant_gene = longest_exon[longest_exon['gene_name'].isin(dominant_list)]
recessive_gene = longest_exon[longest_exon['gene_name'].isin(recessive_list)]
additive_gene = longest_exon[longest_exon['gene_name'].isin(additive_list)]

# Calculate the total length of longest transcripts.
total_dominant_exon_length = dominant_gene["exon_length"].sum()
total_recessive_exon_length = recessive_gene["exon_length"].sum()
total_additive_exon_length = additive_gene["exon_length"].sum()
total_universe_exon_length = universe_gene["exon_length"].sum()

#Calculate the final ratio.
all_chromosome_total = total_additive_exon_length+total_dominant_exon_length+total_recessive_exon_length
all_additive_ratio = total_additive_exon_length/all_chromosome_total
all_dominant_ratio = total_dominant_exon_length/all_chromosome_total
all_recessive_ratio = total_recessive_exon_length/all_chromosome_total


print(f"for all the genome: dominant:additive:recessive={all_dominant_ratio}:{all_additive_ratio}:{all_recessive_ratio}")
