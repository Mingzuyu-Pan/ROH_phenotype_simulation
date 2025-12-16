from sys import argv
import msprime, pyslim, gzip,tskit
import numpy as np 

script, infile = argv
#infile = "gravel.del.only.h00.2281253538061.trees"
outfilep1 = infile + ".p1.vcf.gz"
outfilep2 = infile + ".p2.vcf.gz"
outfilep3 = infile + ".p3.vcf.gz"
mutfilep1 = infile + ".p1.muts.vcf.gz"
mutfilep2 = infile + ".p2.muts.vcf.gz"
mutfilep3 = infile + ".p3.muts.vcf.gz"
#mutfile = infile + ".muts.all.vcf.gz"

ts = tskit.load(infile)

#ensure coalescence of all lineages back in time
recap = pyslim.recapitate(ts,recombination_rate=1e-8, ancestral_Ne=7310, random_seed=np.random.randint(1,100000000))

#build lists for all the ids in each population
p1 = []
p2 = []
p3 = []

for pop,ind in zip(recap.individual_populations,pyslim.individuals_alive_at(recap, 0)):
	if pop == 1:
		p1.append(ind)
	elif pop == 2:
		p2.append(ind)
	elif pop == 3:
		p3.append(ind)

p1 = np.array(p1)
p2 = np.array(p2)
p3 = np.array(p3)

#sample 500 random individuals
p1subsample = np.random.choice(p1, size=500, replace=False) 
p2subsample = np.random.choice(p2, size=500, replace=False) 
p3subsample = np.random.choice(p3, size=500, replace=False)

#build lists that contain the nodes of each individual
p1subsamplenodes = []
p2subsamplenodes = []
p3subsamplenodes = []

for i in np.sort(p1subsample):
	p1subsamplenodes.append(recap.individual(i).nodes[0])
	p1subsamplenodes.append(recap.individual(i).nodes[1])

for i in np.sort(p2subsample):
	p2subsamplenodes.append(recap.individual(i).nodes[0])
	p2subsamplenodes.append(recap.individual(i).nodes[1])

for i in np.sort(p3subsample):
	p3subsamplenodes.append(recap.individual(i).nodes[0])
	p3subsamplenodes.append(recap.individual(i).nodes[1])

p1subsamplenodes = np.array(p1subsamplenodes) 
p2subsamplenodes = np.array(p2subsamplenodes) 
p3subsamplenodes = np.array(p3subsamplenodes) 


#treeseq for each population sample
ts_p1 = recap.simplify(p1subsamplenodes)
ts_p2 = recap.simplify(p2subsamplenodes)
ts_p3 = recap.simplify(p3subsamplenodes)


#output to VCF, this could contain only the deleterious mutations
with gzip.open(mutfilep1, "wt") as vcf_file:
	ts_p1.write_vcf(vcf_file, allow_position_zero=True)
	#ts.write_vcf(output=vcf_file,individuals=p1subsample)

with gzip.open(mutfilep2, "wt") as vcf_file:
	ts_p2.write_vcf(vcf_file, allow_position_zero=True)
	#ts.write_vcf(output=vcf_file,individuals=p2subsample)

with gzip.open(mutfilep3, "wt") as vcf_file:
	ts_p3.write_vcf(vcf_file, allow_position_zero=True)
	#ts.write_vcf(output=vcf_file,individuals=p3subsample)

#add neutral mutations keep=Flase so deleterious mutations won't show up here
mutated = msprime.mutate(recap, rate=2.124e-8, random_seed=np.random.randint(1,100000000), keep=False)

#treeseq for each population sample
ts_p1 = mutated.simplify(p1subsamplenodes)
ts_p2 = mutated.simplify(p2subsamplenodes)
ts_p3 = mutated.simplify(p3subsamplenodes)

#output to VCF, contains only neutral mutations
with gzip.open(outfilep1, "wt") as vcf_file:
	ts_p1.write_vcf(output=vcf_file, allow_position_zero=True)

with gzip.open(outfilep2, "wt") as vcf_file:
	ts_p2.write_vcf(output=vcf_file, allow_position_zero=True)

with gzip.open(outfilep3, "wt") as vcf_file:
	ts_p3.write_vcf(output=vcf_file, allow_position_zero=True)

