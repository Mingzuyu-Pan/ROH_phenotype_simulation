#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:42:39 2025

@author: mzp5919
"""

import pandas as pd
import numpy as np

# To read the hXX list the population list
hXX_list = ['h00','h05','h10','h15']
pop_list = ['p1','p2','p3']

# Set up the cols for summary results
summary_cols = ['hXX','replicate', 'pop', 'A_0', 'A_1', 'A_2','A_nan','B_0', 'B_1', 'B_2', 'B_nan', 'C_0', 'C_1', 'C_2', 'C_nan', 'NONE_0','NONE_1','NONE_2', 'NONE_nan']
final_summary_df = pd.DataFrame(columns=summary_cols)

# The length of chromosome in cM
lenth_of_region = 127.8

for n in range(len(hXX_list)):
    # Set up the folder path of files
    file_path_roh = './'+str(hXX_list[n])+'/' 
    file_path_list = './'+str(hXX_list[n])+'_500.list'   
    # To read the simulation/replicate ID from hXX_500.list
    simID = []
    with open(file_path_list, "r", encoding="utf-8") as file:
        for line in file:
            simID.append(line.strip())
    
    # Start calculation from each replicate
    for s in range(len(simID)):
        # Then consider three populations in each replicate
        for p in range(len(pop_list)):
            # Read the file containing length of different classes of ROH for all the individuals in each replicate
            ROH_path = file_path_roh + str(simID[s]) + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+ '.trees.'+str(pop_list[p])+'.roh.bed.frac.gz'
            # Read the file containing phenotype score contributed by different ROH classes for all the individuals in each repplicate
            Pheno_path = file_path_roh + str(simID[s]) + '/pheno/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+ '.rho.100.tau.100.'+str(pop_list[p])+'.pheno.normalized.txt'
            
            # To store the current simID, popID, and hXX ID
            current_replicate = simID[s]
            current_pop = pop_list[p]
            current_hXX = hXX_list[n]
            
            # Read these two files
            ROH = pd.read_csv(ROH_path, sep=" ")
            Pheno = pd.read_csv(Pheno_path, sep=" ")           
            
            # Drop pop column and calculate length of Non-ROH regions
            del ROH['pop']
            ROH['NONE'] = lenth_of_region - ROH['TOTAL']
            
            # Drop 'TOTAL' ROH length column and the columns containing negative phenotype scores
            del ROH['TOTAL']
            columns_to_drop = Pheno.columns[Pheno.columns.str.endswith('stab')]
            Pheno = Pheno.drop(columns=columns_to_drop)
            
            # Start from dealing with ROH file
            # For each ROH class, if an individual has any non-zero ROH, its status will be set to 'Valid', otherwise, it will be set to '0'.
            # The type of '0' we set here is string instead of float/int
            ROH['ROH_A_status'] = np.where(ROH['A'] == 0,0,'Valid')
            ROH['ROH_B_status'] = np.where(ROH['B'] == 0,0,'Valid')
            ROH['ROH_C_status'] = np.where(ROH['C'] == 0,0,'Valid')
            ROH['ROH_NONE_status'] = np.where(ROH['NONE'] == 0,0,'Valid')
            
            # Same protocol for phenotype score
            Pheno['Pheno_A_status'] = np.where(Pheno['phenoApos'] == 0,0,'Valid')
            Pheno['Pheno_B_status'] = np.where(Pheno['phenoBpos'] == 0,0,'Valid')
            Pheno['Pheno_C_status'] = np.where(Pheno['phenoCpos'] == 0,0,'Valid')
            Pheno['Pheno_NONE_status'] = np.where(Pheno['phenoNONEpos'] == 0,0,'Valid')
            
            # We merged these two files together and dropped the raw data
            merged = pd.concat([ROH, Pheno], axis=1)
            merged_filtered = merged.drop(columns=['A','B','C','NONE','phenoApos','phenoBpos','phenoCpos','phenoNONEpos'])
            
            # Define a new function
            # If an individual doesn't have either a specific class of ROH nor phenotype score contributed by this ROH class, 
            # it will be set to "0"
            # If it has a specific class of ROH but no phenotype score contributed by this ROH class,
            # it will be set to "1"
            # If it has both ROH class and the phenotype score contributed by the ROH class,
            # It will be set to "2"
            # Otherwise, it will be set to NaN. For example, it doesn't have a specific class of ROH, but that ROH class contributed to phenotype score. 
            # This function can cover multiple ROH classes and return value.
            def apply_selection_logic(df, suffix, choices, default=np.nan):
                col_roh = f'ROH_{suffix}_status'
                col_pheno = f'Pheno_{suffix}_status'
                col_result = f'{suffix}_result'
                
                conditions = [
                    (df[col_roh] == '0') & (df[col_pheno] == '0'),
                    (df[col_roh] == 'Valid') & (df[col_pheno] == '0'),
                    (df[col_roh] == 'Valid') & (df[col_pheno] == 'Valid')
                ]
                
                df[col_result] = np.select(conditions, choices, default=default)
            
            # Define the ROH classes we want to consider
            # Define the return values
            suffixes = ['A', 'B', 'C', 'NONE']
            choices = [0, 1, 2]
            
            for q in suffixes:
                apply_selection_logic(merged_filtered, q, choices)
                
            return_values = [0, 1, 2, np.nan]
            
            # To count the number of each return values
            count_A = merged_filtered['A_result'].value_counts().reindex(return_values, fill_value=0)
            count_B = merged_filtered['B_result'].value_counts().reindex(return_values, fill_value=0)
            count_C = merged_filtered['C_result'].value_counts().reindex(return_values, fill_value=0)
            count_NONE = merged_filtered['NONE_result'].value_counts().reindex(return_values, fill_value=0)
            
            new_row_data = {
                'hXX': [current_hXX],
                'replicate': [current_replicate],
                'pop': [current_pop],
                # Note the index is the return value
                'A_0': [count_A.loc[0]],
                'A_1': [count_A.loc[1]],
                'A_2': [count_A.loc[2]],
                'A_nan': [count_A.loc[np.nan]],
                'B_0': [count_B.loc[0]],
                'B_1': [count_B.loc[1]],
                'B_2': [count_B.loc[2]],
                'B_nan': [count_B.loc[np.nan]],
                'C_0': [count_C.loc[0]],
                'C_1': [count_C.loc[1]],
                'C_2': [count_C.loc[2]],
                'C_nan': [count_C.loc[np.nan]],
                'NONE_0': [count_NONE.loc[0]],
                'NONE_1': [count_NONE.loc[1]],
                'NONE_2': [count_NONE.loc[2]],
                'NONE_nan': [count_NONE.loc[np.nan]],
            
            }
            # Add this new row to final dataframe
            new_row_df = pd.DataFrame(new_row_data)
            final_summary_df = pd.concat([final_summary_df, new_row_df], ignore_index=True)
            
# Loop is finished. Output the dataframe to csv file.
final_summary_df.to_csv(
    'calculate_0_proportion_deleterious.csv', 
    index=False,   
    sep=',',       
    encoding='utf-8' 
)