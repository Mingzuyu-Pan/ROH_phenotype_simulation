#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:04:22 2024

@author: mzyp
"""

import pandas as pd
import numpy as np
from sys import argv

script, hXX = argv

sim_id = ['sim_id']
total = ['total']
p1 = ['p1']
p2 = ['p2']
p3 = ['p3']


pd.set_option("display.max_rows", None)
pd.reset_option("display.max_columns")

sim_list = './'+hXX+'_500.list'    
simID = []
with open(sim_list, "r", encoding="utf-8") as file:
    for line in file:
        simID.append(line.strip())

for i in range(len(simID)): 
   path = './'+hXX+'/'+str(simID[i])+'/gravel.del.only.'+hXX+'.'+str(simID[i])+'.muts.gz'
   df = pd.read_csv(path, compression='gzip', sep='\s+')
   total_non_zero = np.count_nonzero(df['freq'])
   p1_non_zero = np.count_nonzero(df['p1_freq'])
   p2_non_zero = np.count_nonzero(df['p2_freq'])
   p3_non_zero = np.count_nonzero(df['p3_freq'])
   
   sim_id.append(str(simID[i]))
   total.append(total_non_zero)
   p1.append(p1_non_zero)
   p2.append(p2_non_zero)
   p3.append(p3_non_zero)
   
summary = pd.DataFrame(list(zip(sim_id,total,p1,p2,p3)))
filename = f"{hXX}_del_count.xlsx"
summary.to_excel(filename, index=False, header=False)  