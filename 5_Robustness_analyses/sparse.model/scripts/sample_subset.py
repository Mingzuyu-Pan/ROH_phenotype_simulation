#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:46:15 2026

@author: mzyp
"""
import random
import gzip


x_list = [1.0, 0.3]
pop_list = ['p1', 'p2', 'p3']
hXX_list = ['h15']

for hXX in hXX_list:
    sim_list = []
    sim_list_path = f'{hXX}_500.list'
    # Import all the simulation id into a list
    with open(sim_list_path, "r", encoding="utf-8") as file:
        for line in file:
            sim_list.append(line.strip())
            
    for simID in sim_list:
        for x in x_list:

            input_filename = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.muts.gz'
            
            output_filename = f'./{hXX}/{simID}/gravel.del.only.{hXX}.{simID}.subset.{x}.muts.gz'
            
            
            with gzip.open(input_filename, 'rt') as f_in, gzip.open(output_filename, 'wt') as f_out:
                # just read a single line, so we could get the header of original file.
                header = f_in.readline()
                # Keep header to the outfile
                f_out.write(header)
                for line in f_in:
                    if random.random() < x:
                        f_out.write(line)