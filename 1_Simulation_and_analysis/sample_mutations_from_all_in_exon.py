#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:21:32 2025

@author: mzyp
"""

import random
from collections import Counter, defaultdict
from sys import argv

script, hXX, slimID = argv

# hXX = 'h00'
# seed = '1021574717937228960'
positions_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.exon.in.sampled.postion.txt"


with open(positions_file, "r") as f:
    position_counts = Counter(int(line.strip()) for line in f)
print(f"loaded {sum(position_counts.values())}  position")

for w in range(3):
    vcf_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only."+str(hXX)+"."+str(slimID)+".trees.p"+str(w+1)+".exon.in.vcf"
    output_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p" + str(w+1) + ".exon.in.smapled.vcf"

    vcf_records = defaultdict(list)
    with open(vcf_file, "r") as infile:
        header_lines = [] 
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line) 
                continue
            
            fields = line.strip().split("\t")
            position = int(fields[1])
            vcf_records[position].append(line) 
    

    with open(output_file, "w") as outfile:
        outfile.writelines(header_lines)

        for position, count in position_counts.items():
            if position in vcf_records:
                available_records = vcf_records[position]
                if len(available_records) >= count:
                    selected_records = random.sample(available_records, count)
                else:
                    selected_records = available_records

                
                outfile.writelines(selected_records)

    print(f"The result is saved in {output_file}")