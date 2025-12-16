#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:11:32 2024

@author: mzp5919
"""
import os
from sys import argv

script,  hXX = argv

def list_folders(directory):
    items = os.listdir(directory)
    folders = [item for item in items if os.path.isdir(os.path.join(directory, item))]
    return folders

directory_path = "./"+str(hXX)
folders = list_folders(directory_path)

output_file = f"{hXX}_500.list"

with open(output_file, "w") as f:
    for folder in folders:
        print(folder)  
        f.write(folder + "\n")  

print(f"The data is saved: {output_file}")