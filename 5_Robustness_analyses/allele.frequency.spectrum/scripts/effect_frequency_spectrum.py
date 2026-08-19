#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 13:06:06 2026

@author: mzyp
"""
import pandas as pd
import numpy as np

hXX_list = ['h00']
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
            muts_df['effect_0.7'] = muts_df['s'].abs()**0.7
            muts_df['effect_1.0'] = muts_df['s'].abs()**1.0
            muts_df['effect_1.3'] = muts_df['s'].abs()**1.3
            # From 10-6 to 1, 50 bins, so 51 values.
            log_bins = np.logspace(-6, 0, num=51)
            # Add a very small value to 1 to deal with the data close to boundary
            log_bins[-1] = log_bins[-1] + 1e-6 
            # Just collect the column names
            effect_cols =['effect_0.7', 'effect_1.0', 'effect_1.3']
            # Only consider the mutations really exist in the pop
            muts_df_p1 = muts_df[muts_df['p1_freq'] > 0].copy()
            # The interval is (a,b].
            muts_df_p1['freq_bin'] = pd.cut(muts_df_p1['p1_freq'], bins=log_bins, right=False)
            # Calculate mean effect in each bin. 
            # Show all bins, though some of them may not have any observed values.
            result_p1 = muts_df_p1.groupby('freq_bin', observed=False)[effect_cols].mean()
            # Change column names of result dataframe
            result_p1 = result_p1.add_suffix(f'_p1_sim{simID}') 
            
            muts_df_p2 = muts_df[muts_df['p2_freq'] > 0].copy()
            muts_df_p2['freq_bin'] = pd.cut(muts_df_p2['p2_freq'], bins=log_bins, right=False)
            result_p2 = muts_df_p2.groupby('freq_bin', observed=False)[effect_cols].mean()
            result_p2 = result_p2.add_suffix(f'_p2_sim{simID}') 
            
            muts_df_p3 = muts_df[muts_df['p3_freq'] > 0].copy()
            muts_df_p3['freq_bin'] = pd.cut(muts_df_p3['p3_freq'], bins=log_bins, right=False)
            result_p3 = muts_df_p3.groupby('freq_bin', observed=False)[effect_cols].mean()
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
        
        effects =['effect_0.7', 'effect_1.0', 'effect_1.3']
        pops = ['p1', 'p2', 'p3']
        
        for eff in effects:
            for p in pops:
                target_prefix = f'{eff}_{p}'
                # Looking for all the comlumns has same pop and tau value
                target_cols = final_df.filter(like=target_prefix).columns
                # Working on this temp df
                target_data = final_df[target_cols]
                # Calculate mean and se of each row, i.e., each bin.
                final_summary_df[f'{target_prefix}_mean'] = target_data.mean(axis=1)
                final_summary_df[f'{target_prefix}_sem']  = target_data.sem(axis=1)
                
        final_summary_df.to_csv(f"effect_allele_frequency_spectrum_casual_effect_{hXX}.subset.{x}.csv", index=False)