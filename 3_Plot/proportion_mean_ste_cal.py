#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 09:11:29 2025

@author: mzyp
"""

import pandas as pd
import numpy as np

hXX_list = ['h00','h05','h10','h15']

for n in range(len(hXX_list)):
    del_path = './'+str(hXX_list[n])+'_summary_proportion_neutral_normalized.xlsx'              
                
    df = pd.read_excel(del_path)
    df = df.loc[:, ~df.columns.str.startswith('standard_')]
    mean_cols = [col for col in df.columns if col.startswith('mean_')]

    result = df.groupby(['pop','tau', 'rho'])[mean_cols].agg(['mean', 'std', 'count']).reset_index()

    result.columns = ['{}_{}'.format(col[0], col[1]) for col in result.columns]

    result = result.reset_index()
    result = result.drop(columns=['index'])

    df_mean = result.set_index(['pop_', 'tau_'])
    df_mean = df_mean.reset_index()

    df_mean = df_mean.rename(columns={
        'tau_': 'Tau',
        'rho_': 'Rho',
        'pop_': 'Pop',
        'mean_A_mean': 'Short ROH-mean',
        'mean_A_std': 'Short ROH-standard variation',
        'mean_A_count': 'Short ROH-count',
        'mean_B_mean': 'Medium ROH-mean',
        'mean_B_std': 'Medium ROH-standard variation',
        'mean_B_count': 'Medium ROH-count',
        'mean_C_mean': 'Long ROH-mean',
        'mean_C_std': 'Long ROH-standard variation',
        'mean_C_count': 'Long ROH-count',
        'mean_NONE_mean': 'Non-ROH-mean',
        'mean_NONE_std': 'Non-ROH-standard variation',
        'mean_NONE_count': 'Non-ROH-count',
        })

    df_mean['Pop'] = df_mean['Pop'].replace({
        'p1': 'African',
        'p2': 'European',
        'p3': 'East Asian'
        })
    
    df_mean['Short ROH-standard error'] = df_mean['Short ROH-standard variation'] / np.sqrt(df_mean['Short ROH-count'])
    df_mean['Medium ROH-standard error'] = df_mean['Medium ROH-standard variation'] / np.sqrt(df_mean['Medium ROH-count'])
    df_mean['Long ROH-standard error'] = df_mean['Long ROH-standard variation'] / np.sqrt(df_mean['Long ROH-count'])
    df_mean['Non-ROH-standard error'] = df_mean['Non-ROH-standard variation'] / np.sqrt(df_mean['Non-ROH-count'])

    df_mean.to_csv('proportion_'+str(hXX_list[n])+'_plot_neutral.csv', index=False, float_format='%.7f', encoding='utf-8')
    
    
hXX_list = ['h00','h05','h10','h15']

for n in range(len(hXX_list)):
    del_path = './'+str(hXX_list[n])+'_summary_proportion_delterious_normalized.xlsx'              
                
    df = pd.read_excel(del_path)
    df = df.loc[:, ~df.columns.str.startswith('standard_')]
    mean_cols = [col for col in df.columns if col.startswith('mean_')]

    result = df.groupby(['pop','tau', 'rho'])[mean_cols].agg(['mean', 'std', 'count']).reset_index()

    result.columns = ['{}_{}'.format(col[0], col[1]) for col in result.columns]

    result = result.reset_index()
    result = result.drop(columns=['index'])

    df_mean = result.set_index(['pop_', 'tau_'])
    df_mean = df_mean.reset_index()

    df_mean = df_mean.rename(columns={
        'tau_': 'Tau',
        'rho_': 'Rho',
        'pop_': 'Pop',
        'mean_A_mean': 'Short ROH-mean',
        'mean_A_std': 'Short ROH-standard variation',
        'mean_A_count': 'Short ROH-count',
        'mean_B_mean': 'Medium ROH-mean',
        'mean_B_std': 'Medium ROH-standard variation',
        'mean_B_count': 'Medium ROH-count',
        'mean_C_mean': 'Long ROH-mean',
        'mean_C_std': 'Long ROH-standard variation',
        'mean_C_count': 'Long ROH-count',
        'mean_NONE_mean': 'Non-ROH-mean',
        'mean_NONE_std': 'Non-ROH-standard variation',
        'mean_NONE_count': 'Non-ROH-count',
        })

    df_mean['Pop'] = df_mean['Pop'].replace({
        'p1': 'African',
        'p2': 'European',
        'p3': 'East Asian'
        })
    
    df_mean['Short ROH-standard error'] = df_mean['Short ROH-standard variation'] / np.sqrt(df_mean['Short ROH-count'])
    df_mean['Medium ROH-standard error'] = df_mean['Medium ROH-standard variation'] / np.sqrt(df_mean['Medium ROH-count'])
    df_mean['Long ROH-standard error'] = df_mean['Long ROH-standard variation'] / np.sqrt(df_mean['Long ROH-count'])
    df_mean['Non-ROH-standard error'] = df_mean['Non-ROH-standard variation'] / np.sqrt(df_mean['Non-ROH-count'])

    df_mean.to_csv('proportion_'+str(hXX_list[n])+'_plot_deleterious.csv', index=False, float_format='%.7f', encoding='utf-8')