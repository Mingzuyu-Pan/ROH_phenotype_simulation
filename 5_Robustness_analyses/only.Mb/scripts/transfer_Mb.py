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


pop_list = ['p1','p2','p3']
# Read genetic distance map for interpolation
for pop in pop_list:    
    # Read raw result from GERMLINE
    roh_cM_raw = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.txt', sep='\s+')
    roh_cM_raw['length_Mb'] = (roh_cM_raw['end'] - roh_cM_raw['start'] + 1) / 1e6

    roh_cM_raw['class_Mb'] = pd.cut(
        roh_cM_raw['length_Mb'],
        bins=[-np.inf, 0.25, 1, np.inf],
        labels=['A', 'B', 'C']
    )
    roh_Mb = roh_cM_raw.drop(columns=['class', 'length'])
    roh_Mb = roh_Mb.rename(columns={'length_Mb': 'length', 'class_Mb': 'class'})
    roh_Mb = roh_Mb[['ind', 'start', 'end', 'class', 'length']]
    roh_Mb.to_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.Mb.txt', sep='\t', index=False)

