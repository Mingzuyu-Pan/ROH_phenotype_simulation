#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:44:31 2024

@author: mzp5919
"""
import pandas as pd
import numpy as np


#hXX_list = ['h00','h05','h10']
hXX_list = ['h00','h05','h10','h15']
rhostr = ["100"]
taustr = ["70","80","90","100","110","120","130"]
for n in range(len(hXX_list)):
    file_path_list = './'+str(hXX_list[n])+'_500.list'    
    simID = []
    with open(file_path_list, "r", encoding="utf-8") as file:
        for line in file:
            simID.append(line.strip())
    
    summary_ID = ['simID']
    summary_pop = ['pop']
    summary_rho = ['rho']
    summary_tau = ['tau']
    
    summary_mean_A = ['mean_A']
    summary_std_dev_A = ['standard_deviation_A']
    summary_std_error_A = ['standard_error_A']
    
    summary_mean_B = ['mean_B']
    summary_std_dev_B = ['standard_deviation_B']
    summary_std_error_B = ['standard_error_B']
    
    summary_mean_C = ['mean_C']
    summary_std_dev_C = ['standard_deviation_C']
    summary_std_error_C = ['standard_error_C']
    
    summary_mean_NONE = ['mean_NONE']
    summary_std_dev_NONE = ['standard_deviation_NONE']
    summary_std_error_NONE = ['standard_error_NONE']
            
    for s in range(len(simID)):        
        for j in range(len(rhostr)):
            for k in range(len(taustr)):
                file_path_data = './'+str(hXX_list[n])+'/' + str(simID[s]) + '/pheno'
                p1_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p1.pheno.normalized.txt'
                p2_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p2.pheno.normalized.txt'
                p3_path = file_path_data + '/gravel.del.only.'+str(hXX_list[n])+'.'+str(simID[s])+'.rho.'+str(rhostr[j])+'.tau.'+str(taustr[k])+'.p3.pheno.normalized.txt'                
                
                df_p1 = pd.read_csv(p1_path, sep=" ")
                df_p1['phenoSUMpos'] = df_p1['phenoApos']+df_p1['phenoBpos']+df_p1['phenoCpos']+df_p1['phenoNONEpos']
                df_p1['phenoAportion_pos'] = df_p1['phenoApos']/df_p1['phenoSUMpos']
                df_p1['phenoBportion_pos'] = df_p1['phenoBpos']/df_p1['phenoSUMpos']
                df_p1['phenoCportion_pos'] = df_p1['phenoCpos']/df_p1['phenoSUMpos']
                df_p1['phenoNONEportion_pos'] = df_p1['phenoNONEpos']/df_p1['phenoSUMpos']
                
                p1_A_mean = df_p1['phenoAportion_pos'].mean()
                p1_A_std_dev = df_p1['phenoAportion_pos'].std()
                p1_A_std_error = p1_A_std_dev/(500**0.5)
                
                p1_B_mean = df_p1['phenoBportion_pos'].mean()
                p1_B_std_dev = df_p1['phenoBportion_pos'].std()
                p1_B_std_error = p1_B_std_dev/(500**0.5)
                
                p1_C_mean = df_p1['phenoCportion_pos'].mean()
                p1_C_std_dev = df_p1['phenoCportion_pos'].std()
                p1_C_std_error = p1_C_std_dev/(500**0.5)
                
                p1_NONE_mean = df_p1['phenoNONEportion_pos'].mean()
                p1_NONE_std_dev = df_p1['phenoNONEportion_pos'].std()
                p1_NONE_std_error = p1_NONE_std_dev/(500**0.5)
                
                summary_ID.append(simID[s])
                summary_pop.append('p1')
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                
                summary_mean_A.append(p1_A_mean)
                summary_std_dev_A.append(p1_A_std_dev)
                summary_std_error_A.append(p1_A_std_error)
                
                summary_mean_B.append(p1_B_mean)
                summary_std_dev_B.append(p1_B_std_dev)
                summary_std_error_B.append(p1_B_std_error)
                
                summary_mean_C.append(p1_C_mean)
                summary_std_dev_C.append(p1_C_std_dev)
                summary_std_error_C.append(p1_C_std_error)
                
                summary_mean_NONE.append(p1_NONE_mean)
                summary_std_dev_NONE.append(p1_NONE_std_dev)
                summary_std_error_NONE.append(p1_NONE_std_error)

                
                df_p2 = pd.read_csv(p2_path, sep=" ")
                df_p2['phenoSUMpos'] = df_p2['phenoApos']+df_p2['phenoBpos']+df_p2['phenoCpos']+df_p2['phenoNONEpos']
                df_p2['phenoAportion_pos'] = df_p2['phenoApos']/df_p2['phenoSUMpos']
                df_p2['phenoBportion_pos'] = df_p2['phenoBpos']/df_p2['phenoSUMpos']
                df_p2['phenoCportion_pos'] = df_p2['phenoCpos']/df_p2['phenoSUMpos']
                df_p2['phenoNONEportion_pos'] = df_p2['phenoNONEpos']/df_p2['phenoSUMpos']
                
                p2_A_mean = df_p2['phenoAportion_pos'].mean()
                p2_A_std_dev = df_p2['phenoAportion_pos'].std()
                p2_A_std_error = p2_A_std_dev/(500**0.5)
                
                p2_B_mean = df_p2['phenoBportion_pos'].mean()
                p2_B_std_dev = df_p2['phenoBportion_pos'].std()
                p2_B_std_error = p2_B_std_dev/(500**0.5)
                
                p2_C_mean = df_p2['phenoCportion_pos'].mean()
                p2_C_std_dev = df_p2['phenoCportion_pos'].std()
                p2_C_std_error = p2_C_std_dev/(500**0.5)
                
                p2_NONE_mean = df_p2['phenoNONEportion_pos'].mean()
                p2_NONE_std_dev = df_p2['phenoNONEportion_pos'].std()
                p2_NONE_std_error = p2_NONE_std_dev/(500**0.5)
                
                summary_ID.append(simID[s])
                summary_pop.append('p2')
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                
                summary_mean_A.append(p2_A_mean)
                summary_std_dev_A.append(p2_A_std_dev)
                summary_std_error_A.append(p2_A_std_error)
                
                summary_mean_B.append(p2_B_mean)
                summary_std_dev_B.append(p2_B_std_dev)
                summary_std_error_B.append(p2_B_std_error)
                
                summary_mean_C.append(p2_C_mean)
                summary_std_dev_C.append(p2_C_std_dev)
                summary_std_error_C.append(p2_C_std_error)
                
                summary_mean_NONE.append(p2_NONE_mean)
                summary_std_dev_NONE.append(p2_NONE_std_dev)
                summary_std_error_NONE.append(p2_NONE_std_error)

                
                df_p3 = pd.read_csv(p3_path, sep=" ")
                df_p3['phenoSUMpos'] = df_p3['phenoApos']+df_p3['phenoBpos']+df_p3['phenoCpos']+df_p3['phenoNONEpos']
                df_p3['phenoAportion_pos'] = df_p3['phenoApos']/df_p3['phenoSUMpos']
                df_p3['phenoBportion_pos'] = df_p3['phenoBpos']/df_p3['phenoSUMpos']
                df_p3['phenoCportion_pos'] = df_p3['phenoCpos']/df_p3['phenoSUMpos']
                df_p3['phenoNONEportion_pos'] = df_p3['phenoNONEpos']/df_p3['phenoSUMpos']
                
                p3_A_mean = df_p3['phenoAportion_pos'].mean()
                p3_A_std_dev = df_p3['phenoAportion_pos'].std()
                p3_A_std_error = p3_A_std_dev/(500**0.5)
                
                p3_B_mean = df_p3['phenoBportion_pos'].mean()
                p3_B_std_dev = df_p3['phenoBportion_pos'].std()
                p3_B_std_error = p3_B_std_dev/(500**0.5)
                
                p3_C_mean = df_p3['phenoCportion_pos'].mean()
                p3_C_std_dev = df_p3['phenoCportion_pos'].std()
                p3_C_std_error = p3_C_std_dev/(500**0.5)
                
                p3_NONE_mean = df_p3['phenoNONEportion_pos'].mean()
                p3_NONE_std_dev = df_p3['phenoNONEportion_pos'].std()
                p3_NONE_std_error = p3_NONE_std_dev/(500**0.5)
                
                summary_ID.append(simID[s])
                summary_pop.append('p3')
                summary_rho.append(rhostr[j])
                summary_tau.append(taustr[k])
                
                summary_mean_A.append(p3_A_mean)
                summary_std_dev_A.append(p3_A_std_dev)
                summary_std_error_A.append(p3_A_std_error)
                
                summary_mean_B.append(p3_B_mean)
                summary_std_dev_B.append(p3_B_std_dev)
                summary_std_error_B.append(p3_B_std_error)
                
                summary_mean_C.append(p3_C_mean)
                summary_std_dev_C.append(p3_C_std_dev)
                summary_std_error_C.append(p3_C_std_error)
                
                summary_mean_NONE.append(p3_NONE_mean)
                summary_std_dev_NONE.append(p3_NONE_std_dev)
                summary_std_error_NONE.append(p3_NONE_std_error)

    df = pd.DataFrame(list(zip(summary_ID,summary_pop,summary_tau,summary_rho,summary_mean_A,summary_std_dev_A,summary_std_error_A,summary_mean_B,summary_std_dev_B,summary_std_error_B,summary_mean_C,summary_std_dev_C,summary_std_error_C,summary_mean_NONE,summary_std_dev_NONE,summary_std_error_NONE)))                
    df.columns = df.iloc[0]
    df = df.drop(index=0).reset_index(drop=True)     
    filename = f"{hXX_list[n]}_summary_proportion_delterious_normalized.xlsx"
    df.to_excel(filename, index=True)
                
                