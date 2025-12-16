#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:22:40 2024

@author: mzp5919
"""
import gzip
from sys import argv
import os

def load_effects(feature_file):
    effect = {}
    h = {}
    with open(feature_file, 'rt') as f:
        for line in f:
            if line.startswith("pos")!= True:
                pos, h_value, z = line.strip().split()
                pos = int(pos)
                effect[pos] = float(z)
                h[pos] = float(h_value)
    return effect, h

def load_roh(roh_file):
    roh = {}
    with gzip.open(roh_file, 'rt') as f:
        ind = None
        for line in f:
            if line.startswith("track"):
                part = line.split()
                ind = part[2]
            else:
                chr_name, start, end, class_name, size, junk_1, junk_2, junk_3, junk_4 = line.strip().split()
                start, end = int(start), int(end) - 1 
                roh.setdefault(ind, {}).setdefault(chr_name, []).append((start, end, class_name))
    return roh

def hits_interval(intervals, pos):
    for start, end, class_name in intervals:
        if start <= pos <= end:
            return class_name
    return '0'



def process_vcf(vcf_file, roh, effect, h):
    stabilizing = {}  
    positive = {}   
    individuals = None

    with open(vcf_file, 'rt') as f:
        for line in f:

            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):

                individuals = line.strip().split()[9:]
                for ind in individuals:
                    stabilizing[ind] = {'A': 0, 'B': 0, 'C': 0, 'NONE': 0}
                    positive[ind] = {'A': 0, 'B': 0, 'C': 0, 'NONE': 0}
            else:

                columns = line.strip().split()
                chr_name = columns[0]
                chr_name = "chr"+str(chr_name)
                pos = int(columns[1])
                genotypes = columns[9:]
                allele_list = {0: 0, 1:pos}

                         
                if pos in effect:
                    for i, genotype in enumerate(genotypes):
                        ind = individuals[i]
                        a1, a2 = genotype.split('|')               
                        if a1 == "." or a2 == ".":
                            continue
                        else:

                            intervals = roh.get(ind, {}).get(chr_name, [])
                            class_name = hits_interval(intervals, pos)

                            if a1 == a2:  # homozygous
                                if a1 == '0':
                                    effect_value = 0
                                else:    
                                    effect_value = effect[pos]
                                    
                                if class_name == '0':  
                                    stabilizing[ind]['NONE'] += effect_value
                                    positive[ind]['NONE'] += abs(effect_value)
                                else:  
                                    stabilizing[ind][class_name] += effect_value
                                    positive[ind][class_name] += abs(effect_value)
                                    
                            else:  # heterozygous
                                mutation_id_a1 = allele_list[int(a1)]
                                effect_value_a1 = effect[mutation_id_a1]*h[mutation_id_a1]
                                mutation_id_a2 = allele_list[int(a2)]
                                effect_value_a2 = effect[mutation_id_a2]*h[mutation_id_a2]
                                
                                if class_name == '0':
                                    stabilizing[ind]['NONE'] += effect_value_a1 + effect_value_a2
                                    positive[ind]['NONE'] += abs(effect_value_a1) + abs(effect_value_a2)
                                else: 
                                    stabilizing[ind][class_name] += effect_value_a1 + effect_value_a2
                                    positive[ind][class_name] += abs(effect_value_a1) + abs(effect_value_a2)       
                else:
                        print(f"{pos} not found in effect dictionary")
                                    

    return stabilizing, positive


script,  hXX, slimID, rho, tau = argv
rho_value = float(int(rho)/100)
tau_value = float(int(tau)/100)

output_folder = "./"+str(hXX)+"/"+str(slimID)+"/pheno_neutral_recalculation"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for w in range(3):
    feature_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p" + str(w+1)+".exon.in.smapled.z.score.rho."+str(rho_value)+".tau."+str(tau_value)+".txt"
    roh_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p"+str(w+1)+".roh.bed.gz"
    vcf_file = "./"+str(hXX)+"/"+str(slimID)+"/gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p"+str(w+1)+".exon.in.smapled.vcf"
    output_file = os.path.join(output_folder, "gravel.del.only." + str(hXX) + "." + str(slimID) + ".trees.p" + str(w+1)+".exon.in.smapled.neutral.mutation.rho."+str(rho_value)+".tau."+str(tau_value)+".pheno.txt")
    effect, h = load_effects(feature_file)
    effect[0]=0
    h[0]= 0
    roh = load_roh(roh_file)
    stabilizing, positive = process_vcf(vcf_file, roh,effect,h)
    with open(output_file, "w") as outfile:
        outfile.write("phenoAstab phenoBstab phenoCstab phenoNONEstab phenoApos phenoBpos phenoCpos phenoNONEpos\n")
    
        for ind in stabilizing:
            outfile.write(f"{ind} "
                          f"{stabilizing[ind]['A']:.15f} "
                          f"{stabilizing[ind]['B']:.15f} "
                          f"{stabilizing[ind]['C']:.15f} "
                          f"{stabilizing[ind]['NONE']:.15f} "
                          f"{positive[ind]['A']:.15f} "
                          f"{positive[ind]['B']:.15f} "
                          f"{positive[ind]['C']:.15f} "
                          f"{positive[ind]['NONE']:.15f}\n")

        print(f"Output has been saved to {output_file}")       
        
        
        
        
        
        
        
        
