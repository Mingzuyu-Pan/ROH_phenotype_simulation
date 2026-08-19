#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:56:03 2026

@author: mzyp
"""

import pandas as pd
import glob
import os
import shutil
from sys import argv

script, x = argv

archive_folder = f"rho_pheno_vcf_revision_raw.subset.{x}"

if not os.path.exists(archive_folder):
    os.makedirs(archive_folder)

result_files = glob.glob(f"roh_pheno_revision_summary.*.subset.{x}.csv")
if result_files:
    final_df = pd.concat([pd.read_csv(f) for f in result_files], ignore_index=True)
    final_output = f"all_simulations_roh_pheno_vcf_revision.subset.{x}.csv"
    final_df.to_csv(final_output, index=False)
    
    for f in result_files:
        shutil.move(f, os.path.join(archive_folder, f))

stats_files = glob.glob(f"different_kind_of_inds.*.subset.{x}.csv")
if stats_files:
    final_stats_df = pd.concat([pd.read_csv(f) for f in stats_files], ignore_index=True)
    stats_output = f"summary_different_kind_of_inds.subset.{x}.csv"
    final_stats_df.to_csv(stats_output, index=False)
    
    for f in stats_files:
        shutil.move(f, os.path.join(archive_folder, f))
