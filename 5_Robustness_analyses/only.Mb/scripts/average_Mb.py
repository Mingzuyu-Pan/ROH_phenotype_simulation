#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 13:46:08 2026

@author: mzyp
"""

import pandas as pd
import numpy as np
from sys import argv

script, hXX, simID = argv

#hXX = 'h00'
#simID = '1021574717937228960'

pop_list = ['p1','p2','p3']
df_all = pd.DataFrame()
for pop in pop_list:    
    roh_cM_raw = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.txt', sep='\s+')
    roh_cM_raw['length_Mb'] = (roh_cM_raw['end'] - roh_cM_raw['start'] + 1) / 1e6
    temp = roh_cM_raw.groupby(['class', 'ind'])[['length_Mb','length']].mean().groupby('class').mean()
    df_result = temp.reset_index()
    df_result['pop'] = pop
    df_result['hXX'] = hXX
    df_result['simID'] = simID
    df_result = df_result[['hXX', 'simID', 'pop', 'class', 'length_Mb','length']]
    df_result = df_result.rename(columns={'length': 'length_cM'})
    df_all = pd.concat([df_all, df_result], ignore_index=True)
df_all.to_csv(f'./gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.length.Mb.in.cM.class.txt', sep='\t', index=False)