import numpy as np
import os

# read GF
inf = open('Swl.dat_old','r')
inf2 = open('../g_bath.dat','r')
outf =open('Gwl.dat','w')

for line in inf:
	dataS  = [float(x) for x in line.split()]
	dataG0 = [float(x) for x in inf2.readline().split()]

	w = dataS[0]
	norb = (len(dataS)-1)/4

	outf.write( str( w ) + '\t')

	for s in range(2):
		for m1 in range(norb):
			for m2 in range(norb):
				if (m1==m2):
					sig = dataS[1+m1*4+2*s+0] + dataS[1+m1*4+2*s+1]*1.0j
					g0 = dataG0[1] + dataG0[2]*1.0j
					g = 1.0/( 1.0/g0 - sig)

					outf.write( str(g.real) + '\t' + str(g.imag) + '\t' )
				else:
					outf.write( str(0.0) + '\t' + str(0.0) + '\t' )

	outf.write('\n')
outf.close()

#os.system("mv Gwl_mat.dat Gwl.dat")
