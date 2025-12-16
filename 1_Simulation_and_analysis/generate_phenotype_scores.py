#!/usr/bin/env python3

import gzip
from sys import argv
import numpy as np


def getZ(s, rho, tau, alpha, mean):

	z = 0

	if np.random.ranf() < rho:
		z = abs(s)**tau
	else:
		z = (np.random.gamma(abs(alpha),scale=abs(mean)/abs(alpha)))**tau

	if np.random.ranf() < 0.5:
		z = -1.0*z

	return z



script, mutsfile, rho, tau, mean, alpha = argv
rho = float(rho)
tau = float(tau)
alpha = float(alpha)
mean = float(mean)


with gzip.open(mutsfile,'rt') as fin:
	for line in fin:
		if line[0] == 'p':
			header = line.strip().split()
			header.insert(4,"z")
			print(" ".join(header))
		else:
			data = line.strip().split()
			z = getZ(float(data[4]),rho,tau,alpha,mean)
			data.insert(4,str(z))
			print(" ".join(data))



