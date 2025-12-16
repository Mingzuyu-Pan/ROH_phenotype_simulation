#!/usr/bin/env python3

import subprocess
from subprocess import PIPE
from sys import argv

script, hXX, slimID = argv

basedir = hXX + "/" + slimID

basename = "gravel.del.only." + hXX + "." + slimID
mutfile = basedir + "/" + basename + ".muts.gz"

subprocess.run("mkdir -p " + basedir + "/pheno", shell=True)

rhostr = ["100"]
taustr = ["70","80","90","100","110","120","130"]

for r in rhostr:
	for t in taustr:
		rho = float(r)/100
		tau = float(t)/100
		phenoMutfile = basedir + "/pheno/" + basename + ".rho." + r + ".tau." + t + ".muts.gz"
		subprocess.run("python generate_phenotype_scores.py " + mutfile + " " + str(rho) + " " + str(tau) + " -0.03 0.2 | gzip -c > " + phenoMutfile, shell=True)
		for p in ["1","2","3"]:
			rohfile = basedir + "/" + basename + ".trees.p" + p + ".roh.bed.gz"
			vcffile = basedir + "/" + basename + ".trees.p" + p + ".muts.vcf.gz"
			phenofile = basedir + "/pheno/" + basename + ".rho." + r + ".tau." + t + ".p" + p + ".pheno.txt"
			subprocess.run("perl compute_phenotype_in_roh_sims.pl " + phenoMutfile + " " + rohfile + " " + vcffile + " " + phenofile, shell=True)



