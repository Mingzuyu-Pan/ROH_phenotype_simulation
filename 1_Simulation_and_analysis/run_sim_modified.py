#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:18:33 2024

@author: mzp5919
"""

#!/usr/bin/env python3

import subprocess
from subprocess import PIPE
from sys import argv

script, slimSim = argv
 
#gravel_chr1_100Mbp_treeseq_del_only_hXX.slim

hXX = slimSim[36:39]
base = "gravel.del.only."+hXX

print("Simulating...\n")
response = subprocess.run("cat "+slimSim+" | /storage/group/zps5164/default/bin/slim", shell=True, stdout=PIPE, stderr=PIPE)#, capture_output=True)
output_lines = response.stdout.decode('UTF-8').splitlines()

for line in output_lines:
    if line.startswith("##simID:"):
        slimID = line[len("##simID:"):]
        print(f"Found line: {line}")
        print(f"Remaining content after '##simID:': {slimID}")
        break  
    
base = base + "." + slimID
basedir = hXX + "/" + slimID

treefile = base + ".trees"
mutfile = base + ".muts"

print("Sampling and generating vcf files...\n")
subprocess.run("python generate_neutral_muts_and_vcf_by_pop2.py " + treefile, shell=True)
subprocess.run("gzip " + treefile, shell=True)

p1base = treefile + ".p1"
p2base = treefile + ".p2"
p3base = treefile + ".p3"

print("Generating tped files...\n")
subprocess.run("/storage/group/zps5164/default/bin/plink --vcf " + p1base + ".vcf.gz --recode transpose --const-fid --maf 0.05 --out " + p1base, shell = True)
subprocess.run("gzip " + p1base + ".tped",shell=True)
subprocess.run("/storage/group/zps5164/default/bin/plink --vcf " + p2base + ".vcf.gz --recode transpose --const-fid --maf 0.05 --out " + p2base, shell = True)
subprocess.run("gzip " + p2base + ".tped",shell=True)
subprocess.run("/storage/group/zps5164/default/bin/plink --vcf " + p3base + ".vcf.gz --recode transpose --const-fid --maf 0.05 --out " + p3base, shell = True)
subprocess.run("gzip " + p3base + ".tped",shell=True)

print("Calling ROH...\n")
subprocess.run("/storage/group/zps5164/default/bin/garlic --tped "+p1base+".tped.gz --tfam "+p1base+".tfam --cm --map plink.chr1.GRCh37.map --size-bounds 0.25 1 --centromere dummy.centromeres.txt --winsize 10 --auto-winsize --auto-overlap-frac --out "+p1base+" --error 0.000001 --lod-cutoff 0 --tped-missing 9", shell=True)
subprocess.run("/storage/group/zps5164/default/bin/garlic --tped "+p2base+".tped.gz --tfam "+p2base+".tfam --cm --map plink.chr1.GRCh37.map --size-bounds 0.25 1 --centromere dummy.centromeres.txt --winsize 10 --auto-winsize --auto-overlap-frac --out "+p2base+" --error 0.000001 --lod-cutoff 0 --tped-missing 9", shell=True)
subprocess.run("/storage/group/zps5164/default/bin/garlic --tped "+p3base+".tped.gz --tfam "+p3base+".tfam --cm --map plink.chr1.GRCh37.map --size-bounds 0.25 1 --centromere dummy.centromeres.txt --winsize 10 --auto-winsize --auto-overlap-frac --out "+p3base+" --error 0.000001 --lod-cutoff 0 --tped-missing 9", shell=True)

print("Cleanup PLINK files...\n")
subprocess.run("rm "+p1base+".tped.gz", shell=True)
subprocess.run("rm "+p1base+".tfam", shell=True)
subprocess.run("rm "+p1base+".nosex", shell=True)
subprocess.run("rm "+p1base+".error", shell=True)
subprocess.run("rm "+p1base+".log", shell=True)
subprocess.run("rm "+p2base+".tped.gz", shell=True)
subprocess.run("rm "+p2base+".tfam", shell=True)
subprocess.run("rm "+p2base+".nosex", shell=True)
subprocess.run("rm "+p2base+".error", shell=True)
subprocess.run("rm "+p2base+".log", shell=True)
subprocess.run("rm "+p3base+".tped.gz", shell=True)
subprocess.run("rm "+p3base+".tfam", shell=True)
subprocess.run("rm "+p3base+".nosex", shell=True)
subprocess.run("rm "+p3base+".error", shell=True)
subprocess.run("rm "+p3base+".log", shell=True)

print("Calculating ROH coverage fractions...\n")
subprocess.run("perl calculate_ROH_fractions.pl " + p1base + ".roh.bed > " + p1base + ".roh.bed.frac",shell=True)
subprocess.run("perl calculate_ROH_fractions.pl " + p2base + ".roh.bed > " + p2base + ".roh.bed.frac",shell=True)
subprocess.run("perl calculate_ROH_fractions.pl " + p3base + ".roh.bed > " + p3base + ".roh.bed.frac",shell=True)

# print("Counting damaging homozygotes...\n")
# subprocess.run("perl count_damaging_genotypes_in_roh_sims.pl " + mutfile + " " + p1base + ".roh.bed " + p1base + ".muts.vcf.gz " + p1base + ".roh.dam.count" ,shell=True)
# subprocess.run("perl count_damaging_genotypes_in_roh_sims.pl " + mutfile + " " + p2base + ".roh.bed " + p2base + ".muts.vcf.gz " + p2base + ".roh.dam.count" ,shell=True)
# subprocess.run("perl count_damaging_genotypes_in_roh_sims.pl " + mutfile + " " + p3base + ".roh.bed " + p3base + ".muts.vcf.gz " + p3base + ".roh.dam.count" ,shell=True)

print("Zipping files...\n")
subprocess.run("gzip " + mutfile, shell=True)
subprocess.run("gzip " + p1base + ".*.kde", shell=True)
subprocess.run("gzip " + p2base + ".*.kde", shell=True)
subprocess.run("gzip " + p3base + ".*.kde", shell=True)
subprocess.run("gzip " + p1base + ".roh.bed", shell=True)
subprocess.run("gzip " + p2base + ".roh.bed", shell=True)
subprocess.run("gzip " + p3base + ".roh.bed", shell=True)
subprocess.run("gzip " + p1base + ".roh.bed.frac", shell=True)
subprocess.run("gzip " + p2base + ".roh.bed.frac", shell=True)
subprocess.run("gzip " + p3base + ".roh.bed.frac", shell=True)
# subprocess.run("gzip " + p1base + ".roh.dam.count", shell=True)
# subprocess.run("gzip " + p2base + ".roh.dam.count", shell=True)
# subprocess.run("gzip " + p3base + ".roh.dam.count", shell=True)

print("Final cleanup...\n")
subprocess.run("mkdir -p " + basedir, shell=True)
subprocess.run("mv " + base + ".* " + basedir + "/",shell=True)
