#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:54:10 2026

@author: mzyp
"""

import pandas as pd
import gzip
import numpy as np
from sys import argv

script, hXX, simID, x = argv

def read_vcf_gz(file_path):
    with gzip.open(file_path, 'rt') as f:
        skip_count = 0
        for line in f:
            if line.startswith('##'):
                skip_count += 1
            else:
                break
    df = pd.read_csv(
        file_path, 
        compression='gzip', 
        sep='\t', 
        #to skip heads. The number of rows to skip is skip_count
        skiprows=skip_count
    )  
    # Drop the columns I don't need
    df = df.drop(['ID', 'REF', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '#CHROM'], axis=1)    
    return df

pop_list = ['p1', 'p2', 'p3']

results_list = []
muts_df = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.subset.{x}.muts.gz', sep='\s+')
for pop in pop_list:
    # Read vcf information of this pop
    vcf_file_path = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.muts.vcf.gz'
    raw_vcf_data = read_vcf_gz(vcf_file_path)
    # Read ROH information of this pop
    roh_file_path = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.txt'
    roh_df_raw = pd.read_csv(roh_file_path, sep='\s+')
    # Modify datatype for following analysis
    raw_vcf_data['POS'] = raw_vcf_data['POS'].astype(int)
    raw_vcf_data['ALT'] = raw_vcf_data['ALT'].astype(str)
    muts_df['pos'] = muts_df['pos'].astype(int)
    muts_df['mid'] = muts_df['mid'].astype(str)
    # Remove non-biallelic locus
    raw_vcf_data = raw_vcf_data[~raw_vcf_data['ALT'].str.contains(',', na=False)].copy()
    # Only choose the intersection of two sets.
    # The reason I do this is, in the vcf file, some mutations are not chosen in the subset.
    merged_vcf = pd.merge(
        raw_vcf_data, 
        muts_df, 
        left_on=['POS', 'ALT'], 
        right_on=['pos', 'mid'], 
        how='inner'
    )
    merged_vcf = merged_vcf.drop(['type'], axis = 1)
    # Modify datatype for following analysis
    merged_vcf['s'] = merged_vcf['s'].astype(float)
    merged_vcf['h'] = merged_vcf['h'].astype(float)
    # Give vcf file a new order
    merged_vcf = merged_vcf.sort_values('POS').reset_index(drop=True)
    # Get positions of all the mutations for folloing mapping analysis
    pos_array = merged_vcf['POS'].values
    # =============================================================================
    # Preparing for effect size calculation    
    # =============================================================================
    # All the genotype columns of all the individuals
    tsk_cols = [col for col in merged_vcf.columns if col.startswith('tsk_')]
    geno_map = {
        '0|0': 0, 
        '0|1': 1, 
        '1|0': 1, 
        '1|1': 2
    }
    # Transfer genotype to the number of alt alleles
    for col in tsk_cols:
        merged_vcf[col] = merged_vcf[col].map(geno_map).astype(int)
    # Number of alt alleles as a matrix
    geno_matrix = merged_vcf[tsk_cols].values
    # Value of h
    h_vec = merged_vcf['h'].values[:, np.newaxis]
    # =============================================================================
    # Preparing for ROH mapping    
    # =============================================================================
    # The name of roh columns
    roh_cols =[col.replace('tsk_', 'roh_') for col in tsk_cols]
    # Define a new data frame. All the values are Non-ROH. It will be change if the mut
    # belongs to any ROH class
    roh_class_df = pd.DataFrame('Non-ROH', index=merged_vcf.index, columns=roh_cols)
    # ROH information is grouped to each ind
    roh_grouped = roh_df_raw.groupby('ind')
    
    for ind in tsk_cols:
        if ind in roh_grouped.groups:
            # Get ROH infor of specific ind
            ind_roh = roh_grouped.get_group(ind)
            # Define the name of this ROH column
            roh_col_name = ind.replace('tsk_', 'roh_')
            # Get the ROH column index of this ind
            col_idx = roh_class_df.columns.get_loc(roh_col_name)
            # All the ROH start position in this ind
            starts = ind_roh['start'].values
            # All the ROH end position in this ind
            ends = ind_roh['end'].values
            # All the ROH class in this ind
            classes = ind_roh['class'].values
            # where to insert start position
            s_indices = np.searchsorted(pos_array, starts, side='left')
            # where to insert end position
            e_indices = np.searchsorted(pos_array, ends, side='right')
            # Each (start, end) combination has a corresponding ROH class
            for s_idx, e_idx, roh_class in zip(s_indices, e_indices, classes):
                if s_idx < e_idx:
                    # All the muts covered by this region has the corresponding class value
                    roh_class_df.iloc[s_idx:e_idx, col_idx] = roh_class
    # Now the vcf file has its ROH class                  
    merged_vcf = pd.concat([merged_vcf, roh_class_df], axis=1)
    
    # =============================================================================
    # Number of causal alleles in different ROH regions and Non-ROH region
    # =============================================================================
    tsk_cols = [col for col in merged_vcf.columns if col.startswith('tsk_')]
    roh_cols = [f"roh_{col.split('_')[1]}" for col in tsk_cols] 
    df_tsk = merged_vcf[tsk_cols]
    #Calculate the total number of causal allele of each individual
    nonzero_per_col = (df_tsk != 0).sum(axis=0)
    summary_results = {}
    for tsk_col in tsk_cols:
        suffix = tsk_col.split('_')[1]     
        roh_col = f'roh_{suffix}'
        
        if roh_col in merged_vcf.columns:
            non_zero_mask = (merged_vcf[tsk_col] != 0)        
            # Only consider the mutations that really exist for each individual
            target_roh_values = merged_vcf.loc[non_zero_mask, roh_col]   
            # Then we can use counts function confidently.
            # Otherwise, for an ind, a mutation can find its position in their 
            # ROH info. But it doesn't mean it really existed in that ind.
            counts = target_roh_values.value_counts()        
            summary_results[tsk_col] = counts
            
    Just_single_allele_df = pd.DataFrame(summary_results).fillna(0).astype(int)
    Just_single_allele_df.loc['total'] = Just_single_allele_df.sum()
    Just_single_allele_df.loc['ROH'] = Just_single_allele_df.loc[['A', 'B', 'C']].sum()
    Just_single_allele_df.loc['prop_ROH'] = Just_single_allele_df.loc['ROH'] / Just_single_allele_df.loc['total']
    Just_single_allele_df.loc['prop_Non_ROH'] = Just_single_allele_df.loc['Non-ROH'] / Just_single_allele_df.loc['total']
    mean_prop_ROH_single = Just_single_allele_df.loc['prop_ROH'].mean()
    mean_prop_Non_ROH_single = Just_single_allele_df.loc['prop_Non_ROH'].mean()
    mean_total_single_count = Just_single_allele_df.loc['total'].mean()
    
    # =============================================================================
    # Number of causal homozygous in different ROH regions and Non-ROH region
    # =============================================================================
    homozygous_per_col = (df_tsk == 2).sum(axis=0)
    summary_results_homo = {}
    for tsk_col in tsk_cols:
        suffix = tsk_col.split('_')[1]     
        roh_col = f'roh_{suffix}'
        
        if roh_col in merged_vcf.columns:
            homo_mask = (merged_vcf[tsk_col] == 2)        
            target_roh_values = merged_vcf.loc[homo_mask, roh_col]        
            counts = target_roh_values.value_counts()        
            summary_results_homo[tsk_col] = counts
            
    homo_allele_df = pd.DataFrame(summary_results_homo).fillna(0).astype(int)
    homo_allele_df.loc['total'] = homo_allele_df.sum()
    homo_allele_df.loc['ROH'] = homo_allele_df.loc[['A', 'B', 'C']].sum()
    homo_allele_df.loc['prop_ROH'] = homo_allele_df.loc['ROH'] / homo_allele_df.loc['total']
    homo_allele_df.loc['prop_Non_ROH'] = homo_allele_df.loc['Non-ROH'] / homo_allele_df.loc['total']
    mean_prop_ROH_homo = homo_allele_df.loc['prop_ROH'].mean()
    mean_prop_Non_ROH_homo = homo_allele_df.loc['prop_Non_ROH'].mean()
    mean_total_homo_count = homo_allele_df.loc['total'].mean()
    results_list.append({
        'hXX': hXX,
        'simID': simID,
        'pop': pop,
        'subset_value': x,
        'prop_ROH_single': mean_prop_ROH_single,
        'prpp_Non_ROH_single': mean_prop_Non_ROH_single,
        'prop_ROH_homo': mean_prop_ROH_homo,
        'prop_Non_ROH_homo': mean_prop_Non_ROH_homo,
        'mean_count_single': mean_total_single_count,
        'mean_count_homo': mean_total_homo_count,
        })

final_results_df = pd.DataFrame(results_list)
final_results_df.to_csv(f"distribution_of_causal_alleles.{hXX}.{simID}.subset.{x}.single_and_homo_all.csv", index=False)