#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:39:38 2024

@author: mzp5919
"""

import pandas as pd
import numpy as np

pd.set_option("display.max_rows", None)
pd.reset_option("display.max_columns")


#hXX_list = ['h00','h05','h10']
hXX_list = ['h00','h05','h10','h15']
rhostr = ["100"]
taustr = ["70","80","90","100","110","120","130"]
for n in range(len(hXX_list)):
    file_path_roh = './'+str(hXX_list[n])+'/' 
    file_path_list = './'+str(hXX_list[n])+'_500.list'    
    simID = []
    with open(file_path_list, "r", encoding="utf-8") as file:
        for line in file:
            simID.append(line.strip())
    
    summary_ID = ['simID']
    summary_pop = ['pop']
    summary_rho = ['rho']
    summary_tau = ['tau']
    
    summary_mean_A_frac = ['mean_A_frac']
    # summary_std_dev_A = ['standard_deviation_A']
    # summary_std_error_A = ['standard_error_A']
    
    summary_mean_B_frac = ['mean_B_frac']
    # summary_std_dev_B = ['standard_deviation_B']
    # summary_std_error_B = ['standard_error_B']
    
    summary_mean_C_frac = ['mean_C_frac']
    # summary_std_dev_C = ['standard_deviation_C']
    # summary_std_error_C = ['standard_error_C']
    
    summary_mean_NONE_frac = ['mean_NONE_frac']
    # summary_std_dev_NONE = ['standard_deviation_NONE']
    # summary_std_error_NONE = ['standard_error_NONE']
            
    for s in range(len(simID)):
        file_path_data = './'+str(hXX_list[n])+'/' + str(simID[s]) + '/pheno'
        
        p1_roh_path = file_path_roh + str(simID[s]) + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+ '.trees.p1.roh.bed.frac.gz'
        p2_roh_path = file_path_roh + str(simID[s]) + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+ '.trees.p2.roh.bed.frac.gz'
        p3_roh_path = file_path_roh + str(simID[s]) + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+ '.trees.p3.roh.bed.frac.gz' 
        
        df_p1_roh = pd.read_csv(p1_roh_path, compression='gzip', sep=' ')
        df_p1_roh["NONE"] = 127.8 - df_p1_roh["TOTAL"]
        df_p1_roh = df_p1_roh.drop(columns=["pop"])
        df_p1_roh = df_p1_roh.drop(columns=["TOTAL"])
        
        df_p2_roh = pd.read_csv(p2_roh_path, compression='gzip', sep=' ')
        df_p2_roh["NONE"] = 127.8 - df_p2_roh["TOTAL"] 
        df_p2_roh = df_p2_roh.drop(columns=["pop"])
        df_p2_roh = df_p2_roh.drop(columns=["TOTAL"])
        
        df_p3_roh = pd.read_csv(p3_roh_path, compression='gzip', sep=' ')
        df_p3_roh["NONE"] = 127.8 - df_p3_roh["TOTAL"]
        df_p3_roh = df_p3_roh.drop(columns=["pop"])
        df_p3_roh = df_p3_roh.drop(columns=["TOTAL"])

        for j in range(len(rhostr)):
            for k in range(len(taustr)):
                p1_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p1.pheno.normalized.txt'
                p2_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p2.pheno.normalized.txt'
                p3_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p3.pheno.normalized.txt'
                
                
                df_p1 = pd.read_csv(p1_path, sep=" ")
                df_p1 = df_p1.drop(columns=["phenoAstab", "phenoBstab","phenoCstab","phenoNONEstab"])
                
                df_p2 = pd.read_csv(p2_path, sep=" ")
                df_p2 = df_p2.drop(columns=["phenoAstab", "phenoBstab","phenoCstab","phenoNONEstab"])

                df_p3 = pd.read_csv(p3_path, sep=" ")
                df_p3 = df_p3.drop(columns=["phenoAstab", "phenoBstab","phenoCstab","phenoNONEstab"])
                
                summary_p1 = pd.DataFrame(index=df_p1.index)
                summary_p1['A_frac'] = df_p1['phenoApos'] / df_p1_roh['A']
                summary_p1['B_frac'] = df_p1['phenoBpos'] / df_p1_roh['B']
                summary_p1['C_frac'] = df_p1['phenoCpos'] / df_p1_roh['C']
                summary_p1['NONE_frac'] = df_p1['phenoNONEpos'] / df_p1_roh['NONE']
                
                summary_p2 = pd.DataFrame(index=df_p2.index)
                summary_p2['A_frac'] = df_p2['phenoApos'] / df_p2_roh['A']
                summary_p2['B_frac'] = df_p2['phenoBpos'] / df_p2_roh['B']
                summary_p2['C_frac'] = df_p2['phenoCpos'] / df_p2_roh['C']
                summary_p2['NONE_frac'] = df_p2['phenoNONEpos'] / df_p2_roh['NONE']
                
                summary_p3 = pd.DataFrame(index=df_p3.index)
                summary_p3['A_frac'] = df_p3['phenoApos'] / df_p3_roh['A']
                summary_p3['B_frac'] = df_p3['phenoBpos'] / df_p3_roh['B']
                summary_p3['C_frac'] = df_p3['phenoCpos'] / df_p3_roh['C']
                summary_p3['NONE_frac'] = df_p3['phenoNONEpos'] / df_p3_roh['NONE']
                
                mean_values_p1 = summary_p1.where((summary_p1 != 0) & summary_p1.notna()).mean()
                mean_values_p2 = summary_p2.where((summary_p2 != 0) & summary_p2.notna()).mean()
                mean_values_p3 = summary_p3.where((summary_p3 != 0) & summary_p3.notna()).mean()
                
                
                summary_ID.append(simID[s])
                summary_pop.append("p1")
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                summary_mean_A_frac.append(mean_values_p1['A_frac'])
                summary_mean_B_frac.append(mean_values_p1['B_frac'])
                summary_mean_C_frac.append(mean_values_p1['C_frac'])
                summary_mean_NONE_frac.append(mean_values_p1['NONE_frac'])
                
                summary_ID.append(simID[s])
                summary_pop.append("p2")
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                summary_mean_A_frac.append(mean_values_p2['A_frac'])
                summary_mean_B_frac.append(mean_values_p2['B_frac'])
                summary_mean_C_frac.append(mean_values_p2['C_frac'])
                summary_mean_NONE_frac.append(mean_values_p2['NONE_frac'])
                
                summary_ID.append(simID[s])
                summary_pop.append("p3")
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                summary_mean_A_frac.append(mean_values_p3['A_frac'])
                summary_mean_B_frac.append(mean_values_p3['B_frac'])
                summary_mean_C_frac.append(mean_values_p3['C_frac'])
                summary_mean_NONE_frac.append(mean_values_p3['NONE_frac'])
    df = pd.DataFrame(list(zip(summary_ID,summary_pop,summary_tau,summary_rho,summary_mean_A_frac,summary_mean_B_frac,summary_mean_C_frac,summary_mean_NONE_frac)))  
    df.columns = df.iloc[0]
    df = df.drop(index=0).reset_index(drop=True)  
    filename = f"{hXX_list[n]}_summary_frac_deleterious_normalized.xlsx"
    df.to_excel(filename, index=True)

                         

                