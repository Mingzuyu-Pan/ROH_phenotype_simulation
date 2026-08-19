#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:03:32 2026

@author: mzyp
"""
import pandas as pd
import gzip
import numpy as np

#define a function to only read non-head lines into dataframe
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

hXX_list = ['h00','h05','h10','h15']
pop_list = ['p1', 'p2', 'p3']
x_list = [0.3, 1.0]
for x in x_list:
    # Need an empty list for each subset proportion
    results_list = []
    for hXX in hXX_list:
        sim_list = []
        sim_list_path = f'{hXX}_500.list'
        # Import all the simulation id into a list
        with open(sim_list_path, "r", encoding="utf-8") as file:
            for line in file:
                sim_list.append(line.strip())
        for simID in sim_list:
            muts_df = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.subset.{x}.muts.gz', sep='\s+')
            for pop in pop_list:
                vcf_file_path = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.muts.vcf.gz'
                raw_vcf_data = read_vcf_gz(vcf_file_path)
                # Transfer the datatype for following comparison
                raw_vcf_data['POS'] = raw_vcf_data['POS'].astype(str)
                raw_vcf_data['ALT'] = raw_vcf_data['ALT'].astype(str)
                muts_df['pos'] = muts_df['pos'].astype(str)
                muts_df['mid'] = muts_df['mid'].astype(str)
                # Remove non-biallelic locus
                raw_vcf_data = raw_vcf_data[~raw_vcf_data['ALT'].str.contains(',', na=False)].copy()
                # Merge genotype info and variant info
                merged_vcf = pd.merge(
                    raw_vcf_data, 
                    muts_df, 
                    # What I need is intersection.
                    how='inner', 
                    left_on=['POS', 'ALT'], 
                    right_on=['pos', 'mid']
                )
                merged_vcf = merged_vcf.drop(['type'], axis = 1)
                merged_vcf['s'] = merged_vcf['s'].astype(float)
                merged_vcf['h'] = merged_vcf['h'].astype(float)
                
                tau_list = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    
                for tau in tau_list:
                    # Calculate the effect of each allele
                    merged_vcf['effect'] = merged_vcf['s'].abs()**tau
                    # All the genotype columns of all the individuals
                    tsk_cols = [col for col in merged_vcf.columns if col.startswith('tsk_')]
                    geno_map = {
                        '0|0': 0, 
                        '0|1': 1, 
                        '1|0': 1, 
                        '1|1': 2
                    }
                    # Replace genotype by the number of alleles
                    merged_vcf[tsk_cols] = merged_vcf[tsk_cols].replace(geno_map).astype(int)
                    # Collect all the number of alleles
                    geno_matrix = merged_vcf[tsk_cols].values
                    # Collect effect size of each single allele
                    effect_vec = merged_vcf['effect'].values[:, np.newaxis]
                    # Collect dominance coefficient of each single allele
                    h_vec = merged_vcf['h'].values[:, np.newaxis]
                    # In following rules:
                    #   if it has two alleles, the score is 1*effect;
                    #   if it has one allele, the score is h*effect;
                    #   if it has 0 allele, the score is 0.
                    #   and this should work for all the scenarios I simulated.
                    pheno_matrix = effect_vec * ( (geno_matrix == 1) * h_vec + (geno_matrix == 2) ) 
                    pheno_cols = [col.replace('tsk_', 'pheno_') for col in tsk_cols]
                    pheno_df = pd.DataFrame(pheno_matrix, columns=pheno_cols, index=merged_vcf.index)
                    final_df = pd.concat([merged_vcf, pheno_df], axis=1)
                    
                    if pop == 'p1':
                        final_df['p1_freq'] = final_df['p1_freq'].astype(float)
                        vector1 = final_df[pheno_cols].sum()
                        # Use the allele frequency in the whole specific population (not joint allele frequency) to decide if the allele is rare or not
                        vector2 = final_df[final_df['p1_freq'] <= 0.001][pheno_cols].sum()
                        vector3 = final_df[final_df['p1_freq'] > 0.01][pheno_cols].sum()
                        vector4 = final_df[(final_df['p1_freq'] > 0.001) & (final_df['p1_freq'] <= 0.01)][pheno_cols].sum()
                        pop_label = 'African'
                    elif pop == 'p2':
                        final_df['p2_freq'] = final_df['p2_freq'].astype(float)
                        vector1 = final_df[pheno_cols].sum()
                        vector2 = final_df[final_df['p2_freq'] <= 0.001][pheno_cols].sum()
                        vector3 = final_df[final_df['p2_freq'] > 0.01][pheno_cols].sum()
                        vector4 = final_df[(final_df['p2_freq'] > 0.001) & (final_df['p2_freq'] <= 0.01)][pheno_cols].sum()
                        pop_label = 'European'
                    elif pop == 'p3':
                        final_df['p3_freq'] = final_df['p3_freq'].astype(float)
                        vector1 = final_df[pheno_cols].sum()
                        vector2 = final_df[final_df['p3_freq'] <= 0.001][pheno_cols].sum()
                        vector3 = final_df[final_df['p3_freq'] > 0.01][pheno_cols].sum()
                        vector4 = final_df[(final_df['p3_freq'] > 0.001) & (final_df['p3_freq'] <= 0.01)][pheno_cols].sum()
                        pop_label = 'East Asian'
                        
                    var_total = vector1.var()
                    var_less_than_0_001 = vector2.var()
                    var_between_0_001_and_0_01 = vector4.var()
                    var_more_than_0_01 = vector3.var()
                    cov_rare_common = vector2.cov(vector3)+vector2.cov(vector4)+vector3.cov(vector4)
                    results_list.append({
                        'hXX': hXX,
                        'simID': simID,
                        'pop': pop_label,
                        'tau': tau,
                        'Total_V_G': var_total,
                        'less_than_0_001_V_G': var_less_than_0_001,
                        'between_0_001_and_0_01_V_G': var_between_0_001_and_0_01,
                        'more_than_0_01_V_G': var_more_than_0_01,
                        'Covariance': cov_rare_common,
                        'Sum_Check (%)': 100*(var_less_than_0_001 + var_between_0_001_and_0_01 + var_more_than_0_01 + 2 * cov_rare_common)/var_total 
                    })
    tau_results_df = pd.DataFrame(results_list)
    tau_results_df.to_csv(f"all_simulations_variance_partition_common_versus_rare.{x}.csv", index=False)