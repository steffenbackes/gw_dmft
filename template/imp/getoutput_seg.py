import numpy as np
import os

# read GF
inf = open('Swl.dat','r')
inf2 = open('../g_bath.dat','r')
outf =open('Gwl.dat','w')

ii = 0
for line in inf:
	dataS  = [float(x) for x in line.split()]
	dataG0 = [float(x) for x in inf2.readline().split()]

	w = dataS[0]
	norb = int( (len(dataS)-1)/4 )
	norbDmft = int( np.sqrt( (len(dataG0)-1)/4 ) )

	if (ii==0):
		print( "norbImp = ",norb )
		print( "norbDMFT = ",norbDmft )
	ii+=1
   
	outf.write( str( w ) + '\t')

	for s in range(2):
		g0 = np.zeros((norb,norb),dtype=complex)
		sig = np.zeros((norb,norb),dtype=complex)
		for m1 in range(norb):
			sig[m1,m1] = dataS[1+m1*4+2*s+0] + dataS[1+m1*4+2*s+1]*1.0j
			for m2 in range(norb):
				g0[m1,m2] = dataG0[1+s*(norbDmft**2)*2+m2*norbDmft*2+m1*2+0] + dataG0[1+s*(norbDmft**2)*2+m2*norbDmft*2+m1*2+1]*1.0j

		g = np.linalg.inv( np.linalg.inv(g0) - sig )

		for m1 in range(norb):
			for m2 in range(norb):
				outf.write( str(g[m1,m2].real) + '\t' + str(g[m1,m2].imag) + '\t' )
	outf.write('\n')
outf.close()

#os.system("mv Gwl_mat.dat Gwl.dat")
