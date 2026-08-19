#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:29:07 2026

@author: mzyp
"""
import pandas as pd
import gzip
import numpy as np
from sys import argv

script, hXX, simID, x = argv

tau_list = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
pop_list = ['p1', 'p2', 'p3']


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

# In some scenario, not all the del mutations are casual allele
# The default setting is all of them are casual, but I introduced x as the paramter to control the proportion of
# del mutations to be casual 
# Yes, read the subset
muts_df = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.subset.{x}.muts.gz', sep='\s+')
all_results_list =[]
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
# Calculate the total length of each roh class for each ind
# =============================================================================
    length_summary = roh_df_raw.groupby(['ind', 'class'])['length'].sum().unstack(fill_value=0)
    
    for roh_class in['A', 'B', 'C']:
        if roh_class not in length_summary.columns:
            length_summary[roh_class] = 0.0
    
    length_summary['Non-ROH'] = 127.8 - length_summary[['A', 'B', 'C']].sum(axis=1)
    
    length_summary = length_summary.add_suffix('_length')
# =============================================================================
# Calculate effect size of each position
# =============================================================================
    for tau in tau_list:
        # Each tau has an empty list to calculate the total phenotype socre of each class
        pheno_summary_list =[]
        # Calculate the effect of each allele
        merged_vcf['effect'] = merged_vcf['s'].abs()**tau
        # Collect effect size of each single allele
        effect_vec = merged_vcf['effect'].values[:, np.newaxis]
        # In following rules:
        #   if it has two alleles, the score is 1*effect;
        #   if it has one allele, the score is h*effect;
        #   if it has 0 allele, the score is 0.
        #   and this should work for all the scenarios I simulated.
        pheno_matrix = effect_vec * ( (geno_matrix == 1) * h_vec + (geno_matrix == 2) ) 
        pheno_cols = [col.replace('tsk_', 'pheno_') for col in tsk_cols]
        pheno_df = pd.DataFrame(pheno_matrix, columns=pheno_cols, index=merged_vcf.index)
        # A full dataframe for this pop and tau.
        final_df = pd.concat([merged_vcf, pheno_df], axis=1)
        
# =============================================================================
# Get the sum of phenotype score for each ROH class
# =============================================================================
        for ind in tsk_cols:
            roh_col = ind.replace('tsk_', 'roh_')
            pheno_col = ind.replace('tsk_', 'pheno_')
            # Get the sum of phenotype score based on the values of roh_col
            sums = final_df.groupby(roh_col)[pheno_col].sum()            
            res_dict = {'ind': ind}
            for roh_class in ['A', 'B', 'C', 'Non-ROH']:
                res_dict[roh_class] = sums.get(roh_class, 0.0)
                
            pheno_summary_list.append(res_dict)
        
        pheno_summary = pd.DataFrame(pheno_summary_list).set_index('ind')
        pheno_summary = pheno_summary.add_suffix('_pheno_sum')
        # Add corresponding roh length colunns. These columns repeat for all the
        # tau values though.
        current_analysis_df = pd.concat([length_summary, pheno_summary], axis=1)
        current_analysis_df['tau'] = tau
        current_analysis_df['pop'] = pop
        # This list collect all the information for this popxtau combination
        all_results_list.append(current_analysis_df)
        
master_analysis_df = pd.concat(all_results_list, ignore_index=False).reset_index(names='ind')
master_analysis_df['tsk_num'] = master_analysis_df['ind'].str.replace('tsk_', '').astype(int)
master_analysis_df = master_analysis_df.sort_values(by=['tau', 'tsk_num'])
master_analysis_df = master_analysis_df.drop(columns=['tsk_num']).reset_index(drop=True)

# =============================================================================
# Proportion of different kinds of inds
# =============================================================================
# Only need a fixed tau value.
target_tau = 1.0
categories =['A', 'B', 'C', 'Non-ROH']
stats_list =[]

target_df = master_analysis_df[master_analysis_df['tau'] == target_tau].copy()

for pop_name, group in target_df.groupby('pop'):
    for roh_class in categories:
        len_col = f'{roh_class}_length'
        score_col = f'{roh_class}_pheno_sum'
        
        cond1_count = ((group[len_col] > 0) & (group[score_col] != 0)).sum()
       
        cond2_count = ((group[len_col] > 0) & (group[score_col] == 0)).sum()        
        cond3_count = ((group[len_col] == 0) & (group[score_col] == 0)).sum()
        
        stats_list.append({
            'pop': pop_name,
            'category': roh_class,
            'len>0_score!=0': cond1_count,
            'len>0_score==0': cond2_count,
            'len==0_score==0': cond3_count,
            'hXX': hXX,
            'simID': simID
        })
stats_df = pd.DataFrame(stats_list)
stats_output_csv = f"different_kind_of_inds.{hXX}.{simID}.subset.{x}.csv"
stats_df.to_csv(stats_output_csv, index=False)

# =============================================================================
# per-cM calculation
# =============================================================================
master_analysis_df['total_abs_pheno'] = (
    master_analysis_df['A_pheno_sum'].abs() +
    master_analysis_df['B_pheno_sum'].abs() +
    master_analysis_df['C_pheno_sum'].abs() +
    master_analysis_df['Non-ROH_pheno_sum'].abs()
)

sd_per_tau = master_analysis_df.groupby('tau')['total_abs_pheno'].std()
cols_to_modify =['A_pheno_sum', 'B_pheno_sum', 'C_pheno_sum', 'Non-ROH_pheno_sum']

for col in cols_to_modify:
    norm_col_name = col.replace('_sum', '_norm')
    master_analysis_df[norm_col_name] = master_analysis_df[col] / master_analysis_df['tau'].map(sd_per_tau)

master_analysis_df = master_analysis_df.drop(columns=['total_abs_pheno'])
master_analysis_df['A_per_cM'] = master_analysis_df['A_pheno_norm']/master_analysis_df['A_length']
master_analysis_df['B_per_cM'] = master_analysis_df['B_pheno_norm']/master_analysis_df['B_length']
master_analysis_df['C_per_cM'] = master_analysis_df['C_pheno_norm']/master_analysis_df['C_length']
master_analysis_df['Non-ROH_per_cM'] = master_analysis_df['Non-ROH_pheno_norm']/master_analysis_df['Non-ROH_length']
# =============================================================================
# Proportion
# =============================================================================
master_analysis_df['total_pheno_norm'] = (
    master_analysis_df['A_pheno_norm'] +
    master_analysis_df['B_pheno_norm'] +
    master_analysis_df['C_pheno_norm'] +
    master_analysis_df['Non-ROH_pheno_norm']
)
master_analysis_df['A_proportion'] = master_analysis_df['A_pheno_norm']/master_analysis_df['total_pheno_norm']
master_analysis_df['B_proportion'] = master_analysis_df['B_pheno_norm']/master_analysis_df['total_pheno_norm']
master_analysis_df['C_proportion'] = master_analysis_df['C_pheno_norm']/master_analysis_df['total_pheno_norm']
master_analysis_df['Non-ROH_proportion'] = master_analysis_df['Non-ROH_pheno_norm']/master_analysis_df['total_pheno_norm']
# =============================================================================
# Aggregation
# =============================================================================
cols_per_cM =['A_per_cM', 'B_per_cM', 'C_per_cM', 'Non-ROH_per_cM']
master_analysis_df[cols_per_cM] = master_analysis_df[cols_per_cM].replace([np.inf, -np.inf], np.nan)

agg_list =[]

for (pop_val, tau_val), group in master_analysis_df.groupby(['pop', 'tau']):
    
    res_dict = {
        'pop': pop_val,
        'tau': tau_val
    }
    
    for cls in['A', 'B', 'C', 'Non-ROH']:
        per_cM_col = f'{cls}_per_cM'
        prop_col = f'{cls}_proportion'        
        res_dict[f'{cls}_per_cM_mean_all'] = group[per_cM_col].mean()
        res_dict[f'{cls}_per_cM_mean_nonzero'] = group[per_cM_col].replace(0.0, np.nan).mean()       
        res_dict[f'{cls}_proportion_mean'] = group[prop_col].mean()
        
    agg_list.append(res_dict)

aggregated_df = pd.DataFrame(agg_list)
aggregated_df['hXX'] = hXX
aggregated_df['simID'] = simID

aggregated_df.to_csv(f"roh_pheno_revision_summary.{hXX}.{simID}.subset.{x}.csv", index=False)