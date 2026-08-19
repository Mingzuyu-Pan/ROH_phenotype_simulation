#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:52:29 2025

@author: mzyp
"""
import re
import gzip
import pandas as pd

hXX_list = ['h00','h05','h10']
track_pattern = re.compile(r"tsk_[0-9]+")

for hXX in hXX_list:
    collection = []
    with open(f'{hXX}_500.list', 'r', encoding='utf-8') as file:
        sim_list =[line.strip() for line in file]
        
    pop_list = ['p1','p2','p3']
    for simID in sim_list:
        for pop in pop_list:
            simID = str(simID)
            file_path_roh = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.bed.gz'
            data_rows = []
            current_track_id = None            
            with gzip.open(file_path_roh, 'rt', encoding='utf-8') as f:
                for line in f:
                    if line.startswith("track"):
                        match = track_pattern.search(line)
                        if match:
                            current_track_id = match.group(0)
                            continue
                    if line:
                        parts = line.split('\t') 
                        parts.append(current_track_id) 
                        data_rows.append(parts)
            
            df_python = pd.DataFrame(data_rows, columns= ['chr', 'start', 'end', 'class', 'cM_length','drop1','drop2','drop3','drop4','ind'])
            df_python = df_python.drop(columns=['chr','drop1','drop2','drop3','drop4'])
            df_python['sim'] = simID
            df_python['pop'] = pop
            new_order = ['sim','ind','pop','start','end','class','cM_length']
            df_python = df_python[new_order]
            df_python[['start','end','cM_length']] = df_python[['start','end','cM_length']].astype(float)
            df_python['bp_length'] = df_python['end']-df_python['start']+1

            summary_df = df_python.groupby(['sim','ind', 'pop','class']).agg(
                nROH=('ind', 'size'),           
                sROH_cM=('cM_length', 'sum'), 
                sROH_bp=('bp_length', 'sum')  
            ).reset_index()
            collection.append(summary_df)
        
    final_summary = pd.concat(collection, ignore_index=True)  
    
    final_summary.to_csv(f'{hXX}_simulated_ROH_Summary.csv', index=False)