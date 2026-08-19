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
gmap = pd.read_csv('./plink.chr1.GRCh37.map', sep='\t', header=None, names=['chr', 'id', 'cM', 'bp'])
for pop in pop_list:    
    # Read raw result from GERMLINE
    germline_raw = pd.read_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.fixed.germline.roh.match', sep='\s+',
                                  header=None, 
                                  names=['familyid1','ind',
                                         'familyid2','indid2',
                                         'chr','start','end',
                                         'dump1','dump2','dump3',
                                         'Mb_length','dump4','dump5','dump6','dump7'])
    # Only choose the column I need
    germline_analysis = germline_raw[['ind','chr','start','end','Mb_length']].copy()
    # Do interpolation to get its cM length
    germline_analysis['start_cM'] = np.interp(germline_analysis['start'], gmap['bp'], gmap['cM'])
    germline_analysis['end_cM'] = np.interp(germline_analysis['end'], gmap['bp'], gmap['cM'])
    germline_analysis['length_cM'] = germline_analysis['end_cM'] - germline_analysis['start_cM']
    germline_analysis['class'] = pd.cut(
        germline_analysis['length_cM'],
        bins=[-np.inf, 0.25, 1, np.inf],
        labels=['A', 'B', 'C'])
    germline_analysis = germline_analysis.drop(columns=['chr', 'Mb_length', 'start_cM', 'end_cM'])
    germline_analysis = germline_analysis.rename(columns={'length_cM': 'length'})
    germline_analysis = germline_analysis[['ind', 'start', 'end', 'class', 'length']]
    germline_analysis['sort_key'] = germline_analysis['ind'].str.split('_').str[1].astype(int)
    germline_analysis = germline_analysis.sort_values(['sort_key', 'start']).drop(columns='sort_key').reset_index(drop=True)
    germline_analysis.to_csv(f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.trees.{pop}.roh.ready.germline.txt', sep='\t', index=False)

