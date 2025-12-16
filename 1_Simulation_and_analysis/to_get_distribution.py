#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 12:37:03 2025

@author: mzyp
"""

import sys
import os
import pandas as pd
#from scipy.stats import norm, anderson
from scipy.stats import norm, shapiro


if len(sys.argv) < 2:
    print("Please offer the value of hXX")
    sys.exit(1)

hXX_list = sys.argv[1:]

output_file = "distribution_fit_results.txt"


with open(output_file, "w") as f:
    f.write("h p_value mu sigma\n")


for hXX in hXX_list:
    file_path = f"./{hXX}_del_count.xlsx"
    df = pd.read_excel(file_path, engine="openpyxl")
    data = df["total"].dropna()

    mu, sigma = norm.fit(data)
    #test_result = anderson(data, dist="norm")
    #p_value = test_result.critical_values[2]
    _, p_value = shapiro(data)


    with open(output_file, "a") as f:
        f.write(f"{hXX} {p_value:.6f} {mu:.6f} {sigma:.6f}\n")

print(f"The result is saved in {output_file}")