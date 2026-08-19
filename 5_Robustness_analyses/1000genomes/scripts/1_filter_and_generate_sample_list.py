#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 13:09:28 2026

@author: mzyp
"""
import pandas as pd
# Read original sample list
original_list_unfiltered = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep= '\s+')
# Help count and check
number_of_original_inds_unfiltered = len(original_list_unfiltered)
# Read original related ind sample list
related_list = pd.read_csv('20140625_related_individuals.txt', sep= '\t')
# Still, help count and check
related_list.columns = related_list.columns.str.strip()
number_of_related_inds = len(related_list)
# Sample column is the name of related individuals
samples_to_remove = related_list['Sample']
# Remove all the related samples
original_list_filtered = original_list_unfiltered[~original_list_unfiltered['sample'].isin(samples_to_remove)]
# Though looks like the original sample list is already good
number_of_original_inds_filtered = len(original_list_filtered)
# Generating sample list of each subgroup based on super pop list.
for(pop_val, super_pop_val), group_temp in original_list_filtered.groupby(['pop', 'super_pop']):
    file_name = f"{pop_val}_{super_pop_val}.txt"
    #Only interested in population belonging to AFR, EAS, and EUR
    # original_list_filtered['super_pop'].unique()
    if ((super_pop_val != 'SAS') & (super_pop_val != 'AMR')):
        group_temp.to_csv(file_name, sep='\t', index=False)
    
# I'm interested in OoA model
target_super_pops = ['AFR', 'EAS', 'EUR']
# Save sample list of untargeted super pop
not_target_super_pop = original_list_filtered[~original_list_filtered['super_pop'].isin(target_super_pops)]
not_target_super_pop.to_csv('SAS_and_AMR.txt', sep='\t', index=False)
# Save sample list of targeted super pops
target_super_pop = original_list_filtered[original_list_filtered['super_pop'].isin(target_super_pops)]
target_super_pop.to_csv('OoA.txt', sep='\t', index=False)
# Save another version for bcftools
target_super_pop_simple = target_super_pop.drop(['gender','super_pop'],axis =1)
target_super_pop_simple.insert(1, 'new_name', '-')
target_super_pop_simple.to_csv('pop_map_bcftools.txt', sep='\t', index=False, header=False)