#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 16:03:03 2025

@author: mzyp
"""
import bisect
import gzip
import random
import numpy as np
import pandas as pd
from sys import argv



def load_parameters(file_path):
    df = pd.read_csv(file_path, sep="\s+", engine="python")  
    param_dict = df.set_index("h")[["mu", "sigma"]].T.to_dict("list")
    return param_dict

    
def is_in_intervals(value, intervals, starts, ends):
    idx = bisect.bisect_right(starts, value) - 1
    if idx >= 0 and starts[idx] <= value <= ends[idx]:
        return True
    return False

def count_lines(file_path):
    with open(file_path, "r") as f:
        return sum(1 for _ in f)




script, hXX, slimID = argv

count_file = 0

exon_file = "./ccdsGene.chr1.100Mbps.exons.uniq.noverlap.txt"
output_txt_summary = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.exon.in.summary.txt"
output_txt_sampled_position = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.exon.in.sampled.postion.txt"

parameters_file = "distribution_fit_results.txt"
parameters = load_parameters(parameters_file)

if hXX in parameters:
    mu, sigma = parameters[hXX]
    print(f"{hXX}: mu = {mu}, sigma = {sigma}")
else:
    print(f"No {hXX} parameters")



with open(output_txt_summary, "w") as outfile_summary:     
    for w in range(3):    
        input_vcf = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.p"+str(w+1)+".vcf"
        output_vcf_in = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.p"+str(w+1)+".exon.in.vcf"
        output_vcf_out = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.p"+str(w+1)+".exon.out.vcf"
        
        intervals = []
        count_special = 0

        with open(exon_file, "r") as file:
            for line in file:
                line = line.strip()      
                start, end = map(int, line.split("\t"))
                intervals.append((start, end))

        starts = [interval[0] for interval in intervals]
        ends = [interval[1] for interval in intervals]


        with open(input_vcf, "r") as infile, \
            open(output_vcf_in, "w") as outfile_in, \
            open(output_vcf_out, "w") as outfile_out:
        
            for line in infile:
                if line.startswith("#"):
                    outfile_in.write(line)
                    outfile_out.write(line)
                    count_special = count_special+1

                else:
                    fields = line.strip().split("\t")
                    position = int(fields[1])            
                    if is_in_intervals(position, intervals, starts, ends):
                        outfile_in.write(line) 
                        outfile_summary.write(str(position)+"\n")
                    else:
                        outfile_out.write(line) 
                    
        count_file = count_lines(output_vcf_in)-count_special + count_file
        print(count_lines(output_vcf_in)-count_special)
    
print(count_lines(output_txt_summary)-count_file)

with open(output_txt_summary, "r") as file:
    numbers = [int(line.strip()) for line in file]


value_sample = np.random.normal(mu,sigma)
value_sample_round = round(value_sample)
sampled_numbers = random.sample(numbers, value_sample_round)

with open(output_txt_sampled_position, "w") as file:
    for number in sampled_numbers:
        file.write(f"{number}\n")

