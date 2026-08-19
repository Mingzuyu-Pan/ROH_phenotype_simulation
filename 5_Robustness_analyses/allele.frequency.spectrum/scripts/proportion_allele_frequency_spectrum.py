#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 13:06:06 2026

@author: mzyp
"""
import pandas as pd
import numpy as np

hXX_list = ['h00','h05','h10','h15']
x_list = ['0.3', '1.0']
for x in x_list:
    for hXX in hXX_list:
        # Each hXX has a new result df
        results_list = []
        sim_list = []
        sim_list_path = f'{hXX}_500.list'
        # Import all the simulation id into a list
        with open(sim_list_path, "r", encoding="utf-8") as file:
            for line in file:
                sim_list.append(line.strip())
        for simID in sim_list:
            # Only need to consider mutation information.
            muts_df = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.subset.{x}.muts.gz', sep='\s+')
            # From 10-6 to 1, 50 bins, so 51 values.
            log_bins = np.logspace(-6, 0, num=51)
            # Add a very small value to 1 to deal with the data close to boundary
            log_bins[-1] = log_bins[-1] + 1e-6 
            # Only consider the mutations really exist in the pop
            muts_df_p1 = muts_df[muts_df['p1_freq'] > 0].copy()
            # The interval is [a,b).
            muts_df_p1['freq_bin'] = pd.cut(muts_df_p1['p1_freq'], bins=log_bins, right=False)
            # Calculate number of allele in each bin. 
            # Show all bins, though some of them may not have any observed values.
            result_p1 = muts_df_p1.groupby('freq_bin', observed=False)[['p1_freq']].count()
            total_p1 = result_p1['p1_freq'].sum()
            result_p1['p1_freq'] = result_p1['p1_freq']/total_p1
            # Change column names of result dataframe
            result_p1 = result_p1.add_suffix(f'_p1_sim{simID}') 
            
            muts_df_p2 = muts_df[muts_df['p2_freq'] > 0].copy()
            muts_df_p2['freq_bin'] = pd.cut(muts_df_p2['p2_freq'], bins=log_bins, right=False)
            result_p2 = muts_df_p2.groupby('freq_bin', observed=False)[['p2_freq']].count()
            total_p2 = result_p2['p2_freq'].sum()
            result_p2['p2_freq'] = result_p2['p2_freq']/total_p2
            result_p2 = result_p2.add_suffix(f'_p2_sim{simID}') 
            
            muts_df_p3 = muts_df[muts_df['p3_freq'] > 0].copy()
            muts_df_p3['freq_bin'] = pd.cut(muts_df_p3['p3_freq'], bins=log_bins, right=False)
            result_p3 = muts_df_p3.groupby('freq_bin', observed=False)[['p3_freq']].count()
            total_p3 = result_p3['p3_freq'].sum()
            result_p3['p3_freq'] = result_p3['p3_freq']/total_p3
            result_p3 = result_p3.add_suffix(f'_p3_sim{simID}') 
            
            #Combine results of three populations
            current_sim_combined = pd.concat([result_p1, result_p2, result_p3], axis=1)
            # Keep the result of each simID into the bigger list
            results_list.append(current_sim_combined)
            
        final_df = pd.concat(results_list, axis=1)
        final_df = final_df.reset_index()
        #Each hXX has its own final_summay_df
        final_summary_df = pd.DataFrame()
        final_summary_df['freq_bin'] = final_df['freq_bin']
        
        p1_cols = [col for col in final_df.columns if col.startswith('p1')]
        final_summary_df['p1_mean'] = final_df[p1_cols].mean(axis=1)
        final_summary_df['p1_sem']  = final_df[p1_cols].sem(axis=1)
        p2_cols = [col for col in final_df.columns if col.startswith('p2')]
        final_summary_df['p2_mean'] = final_df[p2_cols].mean(axis=1)
        final_summary_df['p2_sem']  = final_df[p2_cols].sem(axis=1)        
        p3_cols = [col for col in final_df.columns if col.startswith('p3')]
        final_summary_df['p3_mean'] = final_df[p3_cols].mean(axis=1)
        final_summary_df['p3_sem']  = final_df[p3_cols].sem(axis=1)        
                
        final_summary_df.to_csv(f"proportion_allele_frequency_spectrum_casual_effect_{hXX}.subset.{x}.csv", index=False)