#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:02:47 2025

@author: mzp5919
"""
import numpy as np
from sys import argv

script, hXX, slimID, rho, tau = argv


position_list = ['pos']
h_list = ['h']
z_list = ['z']
s = 0 
rho_value = float(int(rho)/100)
tau_value = float(int(tau)/100)
alpha = 0.2
mean = -0.03


def getZ(s, rho, tau, alpha, mean):
	z = 0
	if np.random.ranf() < rho:
		z = abs(s)**tau
	else:
		z = (np.random.gamma(abs(alpha),scale=abs(mean)/abs(alpha)))**tau

	if np.random.ranf() < 0.5:
		z = -1.0*z

	return z

for w in range(3):
    vcf_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p" + str(w+1) + ".exon.in.smapled.vcf"
    output_pheno_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p" + str(w+1)+".exon.in.smapled.z.score.rho."+str(rho_value)+".tau."+str(tau_value)+".txt"


    with open(vcf_file, 'r') as infile: 
        header_lines = []  
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line)  
                continue
            fields = line.strip().split("\t")
        
            position = int(fields[1])
            z_temp = getZ(s, rho_value, tau_value, alpha, mean)
            
            if hXX == 'h00':
               hXX_value = 0.0
                
            if hXX == 'h05':
               hXX_value = 0.5 
                
            if hXX == 'h10':
               hXX_value = 1.0 
               
            if hXX == 'h15':
                h_possible = [1.0, 0.5, 0.0]
                h_probabilities = [0.0454897924605714, 0.8839168539343177, 0.07059335360511085]
                hXX_value = np.random.choice(h_possible, p=h_probabilities)
            
            if hXX == 'h25':
                h_possible = [1.0, 0.5, 0.0]
                h_probabilities = [0.031447682119205296, 0.9084344370860927, 0.06011788079470198]
                hXX_value = np.random.choice(h_possible, p=h_probabilities)
                
        
            position_list.append(position)
            h_list.append(hXX_value)
            z_list.append(z_temp)
        



    with open(output_pheno_file, "w") as file:
        for row in zip(position_list, h_list, z_list):
            file.write(" ".join(map(str, row)) + "\n")         
        
    print(f"The data is written into {output_pheno_file}")
