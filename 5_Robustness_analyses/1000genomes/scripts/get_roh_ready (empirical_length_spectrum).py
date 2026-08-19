#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:52:29 2025

@author: mzyp
"""
import re
import pandas as pd

all_summary_dfs =[] 
chr1_summary_dfs =[] 
chr1_100M_summary_dfs =[] 

pop_mapping = {
'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
'CEU':'EUR','TSI':'EUR','GBR':'EUR','FIN':'EUR','IBS':'EUR',
'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR','ASW':'AFR','ACB':'AFR'}

with open('pop_list.txt', 'r', encoding='utf-8') as file:
    pop_list =[line.strip() for line in file]

for pop in pop_list:
    file_path_roh = f'./WholeGenome_chr1-22_{pop}_ROH.roh.bed'
    data_rows = []        
    track_pattern = re.compile(r"Ind:\s*([A-Za-z0-9_]+)")
    current_track_id = None
    
    # Save each .bed to a dataframe. Then decide how to process it.    
    with open(file_path_roh, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith("track"):
                match = track_pattern.search(line)
                if match:
                    current_track_id = match.group(1)
                continue
                
            if line:
                parts = line.split('\t') 
                parts.append(current_track_id) 
                data_rows.append(parts)
    # Basic data clean  
    df_python = pd.DataFrame(data_rows, columns= ['chr', 'start', 'end', 'class', 'cM_length','drop1','drop2','drop3','drop4','ind'])
    df_python = df_python.drop(columns=['drop1','drop2','drop3','drop4'])
    df_python['pop'] = pop
    df_python['superpop'] = df_python['pop'].map(pop_mapping)
    new_order = ['chr','ind','pop','superpop','start','end','class','cM_length']
    df_python = df_python[new_order]
    df_python[['start','end','cM_length']] = df_python[['start','end','cM_length']].astype(float)
    df_python['bp_length'] = df_python['end']-df_python['start']+1
    # I want raw roh data, so no further process needed.
    df_chr1 = df_python[df_python['chr'].isin(['chr1'])]
    chr1_summary_dfs.append(df_chr1)
    
    df_chr1_100M = df_chr1[df_chr1['end'] <= 100000000] 
    chr1_100M_summary_dfs.append(df_chr1_100M)
    
    
    all_summary_dfs.append(df_python)
    
final_master_chr1_distr = pd.concat(chr1_summary_dfs, ignore_index=True)        
final_master_chr1_distr.to_csv('All_Populations_chr1_ROH_Summary_distr.csv.gz', index=False, compression='gzip')
final_master_chr1_100M_distr = pd.concat(chr1_100M_summary_dfs, ignore_index=True)        
final_master_chr1_100M_distr.to_csv('All_Populations_chr1_100MB_ROH_Summary_distr.csv.gz', index=False, compression='gzip')
final_master_df_distr = pd.concat(all_summary_dfs, ignore_index=True)        
final_master_df_distr.to_csv('All_Populations_ROH_Summary_distr.csv.gz', index=False, compression='gzip')

