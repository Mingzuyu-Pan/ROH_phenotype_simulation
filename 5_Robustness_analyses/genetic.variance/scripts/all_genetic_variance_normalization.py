#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:44:03 2026

@author: mzyp
"""
import pandas as pd
import gzip
import numpy as np
from sys import argv

#script, hXX, simID, x = argv

tau_list = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
pop_list = ['p1', 'p2', 'p3']

hXX = 'h00'
simID = '103568543989312214'
x = '1.0'

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
all_results_list_uncentered =[]
all_results_list_centered =[]
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
# Adding normalization
# =============================================================================
        tau_value = str(int(float(tau)*100))
        p1_original_pheno_path = f'./{hXX}/{simID}/pheno/gravel.del.only.{hXX}.{simID}.rho.100.tau.{tau_value}.p1.pheno.txt'
        p2_original_pheno_path = f'./{hXX}/{simID}/pheno/gravel.del.only.{hXX}.{simID}.rho.100.tau.{tau_value}.p2.pheno.txt'
        p3_original_pheno_path = f'./{hXX}/{simID}/pheno/gravel.del.only.{hXX}.{simID}.rho.100.tau.{tau_value}.p3.pheno.txt'
        p1_original_pheno = pd.read_csv(p1_original_pheno_path, sep='\s+')
        p1_pheno_vec = p1_original_pheno[['phenoApos', 'phenoBpos', 'phenoCpos', 'phenoNONEpos']].sum(axis=1).values
        p2_original_pheno = pd.read_csv(p2_original_pheno_path, sep='\s+')
        p2_pheno_vec = p2_original_pheno[['phenoApos', 'phenoBpos', 'phenoCpos', 'phenoNONEpos']].sum(axis=1).values
        p3_original_pheno = pd.read_csv(p3_original_pheno_path, sep='\s+')
        p3_pheno_vec = p3_original_pheno[['phenoApos', 'phenoBpos', 'phenoCpos', 'phenoNONEpos']].sum(axis=1).values     
        all_pheno_this_tau = np.concatenate([p1_pheno_vec, p2_pheno_vec, p3_pheno_vec])
        sd_now = np.std(all_pheno_this_tau, ddof=1)
        pheno_cols = [col for col in final_df.columns if col.startswith('pheno_')]
        final_df[pheno_cols] = final_df[pheno_cols] / sd_now
        
# =============================================================================
# Now start the calculation of VA and VD. (normalization and not centered)
# =============================================================================
        tsk_cols = [col for col in final_df.columns if col.startswith('tsk_')]
        pheno_cols = [col for col in final_df.columns if col.startswith('pheno_')]
        roh_cols = [col.replace('tsk_', 'roh_') for col in tsk_cols]
        
        G_matrix = final_df[pheno_cols].values.astype(float)
        ROH_matrix = final_df[roh_cols].values.astype(str)
        
        def get_A_and_D(row):
            # Genotype            
            x = row[tsk_cols].values.astype(float)
            # Genetic value of each locus
            y = row[pheno_cols].values.astype(float)
    
            if x.min() == x.max():
                slope = 0.0
                intercept = np.mean(y)
            else:
            # Linear regression
                slope, intercept = np.polyfit(x, y, 1)
    
            #breeding value (A)
            A_row = intercept + slope * x
            # D = G - A
            D_row = y - A_row    
            
            return pd.Series(np.concatenate([A_row, D_row]))
        
        ad_df = final_df.apply(get_A_and_D, axis=1)
        A_matrix = ad_df.iloc[:, :500].values
        D_matrix = ad_df.iloc[:, 500:].values

        total_A = np.sum(A_matrix, axis=0)
        total_D = np.sum(D_matrix, axis=0)
        total_G = np.sum(G_matrix, axis=0)
        V_A = np.var(total_A, ddof=1)
        V_D = np.var(total_D, ddof=1)        
        V_G = np.var(total_G, ddof=1)
        
        Cov_AD = np.cov(total_A, total_D)[0, 1]
        Sum_Check = V_A + V_D + 2 * Cov_AD - V_G  
# ==================================================================================
# Variance expalined by ROH and non-ROH regions. (normalization and not centered)
# ==================================================================================
        # Variance explained by ROH region
        mask_ROH = (ROH_matrix != 'Non-ROH')    
        ROH_G_sum = np.sum(G_matrix * mask_ROH, axis=0)        
        ROH_V_G = np.var(ROH_G_sum, ddof=1)
        # Variance explained by non-ROH region
        mask_Non_ROH = (ROH_matrix == 'Non-ROH')        
        Non_ROH_G_sum = np.sum(G_matrix * mask_Non_ROH, axis=0)        
        Non_ROH_V_G = np.var(Non_ROH_G_sum, ddof=1)
        cov_ROH = np.cov(ROH_G_sum, Non_ROH_G_sum)[0, 1]
        ROH_check = ROH_V_G + Non_ROH_V_G + 2*cov_ROH - V_G
        all_results_list_uncentered.append({
            'hXX': hXX,
            'simID': simID,
            'pop': pop,
            'tau': tau,
            'subset_value': x,
            'Total_V_G': V_G,
            'Total_V_A': V_A,
            'Total_V_D': V_D,
            'Total_cov': 2 * Cov_AD,
            'Total_check': Sum_Check,
            'ROH_V_G': ROH_V_G,
            'Non_ROH_V_G': Non_ROH_V_G,
            'ROH_cov': cov_ROH,
            'ROH_check': ROH_check
        })                

        
# =============================================================================
# Now start the calculation of VA and VD. (normalization and centered)
# =============================================================================
        tsk_cols = [col for col in final_df.columns if col.startswith('tsk_')]
        pheno_cols = [col for col in final_df.columns if col.startswith('pheno_')]
        roh_cols = [col.replace('tsk_', 'roh_') for col in tsk_cols]

        G_matrix = final_df[pheno_cols].values.astype(float)
        ROH_matrix = final_df[roh_cols].values.astype(str)
        G_mean = np.mean(G_matrix, axis=1, keepdims=True)
        G_matrix_centered = G_matrix - G_mean
        final_df[pheno_cols] = G_matrix_centered
        G_matrix = G_matrix_centered
        
        def get_A_and_D(row):
            # Genotype            
            x = row[tsk_cols].values.astype(float)
            # Genetic value of each locus
            y = row[pheno_cols].values.astype(float)
    
            if x.min() == x.max():
                slope = 0.0
                intercept = np.mean(y)
            else:
            # Linear regression
                slope, intercept = np.polyfit(x, y, 1)
    
            #breeding value (A)
            A_row = intercept + slope * x
            # D = G - A
            D_row = y - A_row    
            
            return pd.Series(np.concatenate([A_row, D_row]))
        
        ad_df = final_df.apply(get_A_and_D, axis=1)
        A_matrix = ad_df.iloc[:, :500].values
        D_matrix = ad_df.iloc[:, 500:].values

        total_A = np.sum(A_matrix, axis=0)
        total_D = np.sum(D_matrix, axis=0)
        total_G = np.sum(G_matrix, axis=0)
        V_A = np.var(total_A, ddof=1)
        V_D = np.var(total_D, ddof=1)        
        V_G = np.var(total_G, ddof=1)
        
        Cov_AD = np.cov(total_A, total_D)[0, 1]
        Sum_Check = V_A + V_D + 2 * Cov_AD - V_G  
# ==================================================================================
# Variance expalined by ROH and non-ROH regions. (normalization and centered)
# ==================================================================================
        # Variance explained by ROH region
        mask_ROH = (ROH_matrix != 'Non-ROH')        
        ROH_G_sum = np.sum(G_matrix * mask_ROH, axis=0)        
        ROH_V_G = np.var(ROH_G_sum, ddof=1)
        # Variance explained by non-ROH region
        mask_Non_ROH = (ROH_matrix == 'Non-ROH')        
        Non_ROH_G_sum = np.sum(G_matrix * mask_Non_ROH, axis=0)        
        Non_ROH_V_G = np.var(Non_ROH_G_sum, ddof=1)
        cov_ROH = np.cov(ROH_G_sum, Non_ROH_G_sum)[0, 1]
        ROH_check = ROH_V_G + Non_ROH_V_G + 2*cov_ROH - V_G
        all_results_list_centered.append({
            'hXX': hXX,
            'simID': simID,
            'pop': pop,
            'tau': tau,
            'subset_value': x,
            'Total_V_G': V_G,
            'Total_V_A': V_A,
            'Total_V_D': V_D,
            'Total_cov': 2 * Cov_AD,
            'Total_check': Sum_Check,
            'ROH_V_G': ROH_V_G,
            'Non_ROH_V_G': Non_ROH_V_G,
            'ROH_cov': cov_ROH,
            'ROH_check': ROH_check
        })       
all_results_df_uncentered = pd.DataFrame(all_results_list_uncentered)                
all_results_df_uncentered.to_csv(f"variance_partition_summary.{hXX}.{simID}.subset.{x}.uncentered.normalization.csv", index=False) 
all_results_df_centered = pd.DataFrame(all_results_list_centered)                
all_results_df_centered.to_csv(f"variance_partition_summary.{hXX}.{simID}.subset.{x}.centered.normalization.csv", index=False)               