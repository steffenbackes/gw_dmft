import numpy as np

#
inf = open('input.ini','r')
for line in inf:
	if ( 'model.beta' in line ):
		beta = float( line.replace('=',' ').split()[1] )
	if ( 'model.U' in line ):
		U= float( line.replace('=',' ').split()[1] )
inf.close()

# read GF
inf = open('sgf.dat','r')
for i in range(3):
	inf.readline()

data = inf.readline().replace(',',' ').split()
nw = int(data[4])
norb = int(data[5])

print( 'Getoutput:' )
print( 'Found beta=',beta )
print( 'Found U=',U )
print( 'Found nw=',nw )
print( 'Found norb=',norb )
inf.readline()

# read gbath
gbath = np.zeros((norb,norb,2,nw),dtype=complex)
inf2 = open('g0_solver.dat','r')
for n in range(nw):
   data = [float(x) for x in inf2.readline().split()]
   norb_gbath = int( np.rint( np.sqrt( (len(data)-1)/4 ) ) )

   # here we assume that the dmft atom in Gbath is the first one

   for s in range(2):
      for m1 in range(norb):
         for m2 in range(norb):
            gbath[m1,m2,s,n] = data[1 + s*2*norb_gbath**2 + m1*2*norb_gbath + m2*2 + 0] + data[1 + s*2*norb_gbath**2 + m1*2*norb_gbath + m2*2 + 1]*1.0j
inf2.close()
###########

# we have to shift the local mu level for ctint
for n in range(nw):
   for s in range(2):
      gbath[:,:,s,n] = np.linalg.inv( np.linalg.inv( gbath[:,:,s,n] ) - 0.5*U*np.identity(norb) )

outf = open('Gwl.dat','w')
for n in range(nw):
	outf.write( "{:.7e}".format( (2*n+1)*np.pi/beta ) + '\t')

	SigG = np.zeros((norb,norb,2),dtype=complex)

	for m1 in range(norb):
		for m2 in range(norb):
			for s in range(2):
				data = [float(x) for x in inf.readline().replace(',',' ').split()[5:] ] 
				SigG[m1,m2,s] = data[0] + data[1]*1.0j

	for s in range(2):
		gimp = gbath[:,:,s,n] - np.dot( gbath[:,:,s,n] , SigG[:,:,s] )	

		for m1 in range(norb):
			for m2 in range(norb):
				outf.write("{:.7e}".format(gimp[m1,m2].real) + '\t' + "{:.7e}".format(gimp[m1,m2].imag) + '\t') 

	outf.write('\n')
outf.close()
inf.close()

# legendre coeffs
inf = open('gf_legendre.dat','r')
for i in range(3):
	inf.readline()

data = inf.readline().replace(',',' ').split()
nlegendre = int(data[4])
print( 'Found nlegendre=',nlegendre )

legendre = np.zeros((nlegendre,norb,norb,2))
inf.readline()

for l in range(nlegendre):
	for m1 in range(norb):
		for m2 in range(norb):
			for s in range(2):
				data = [float(x) for x in inf.readline().replace(',',' ').split()[5:] ] 
				legendre[l,m1,m2,s] = data[0]
inf.close()

outf = open('legendre_coeffs.dat','w')
for l in range(nlegendre):
	outf.write( str(l) + '\t' )
	for m1 in range(norb):
		for m2 in range(norb):
			for s1 in range(2):		
				outf.write( "{:.7e}".format(legendre[l,m1,m2,s1]) + '\t' )
	outf.write( '\n' )

outf.close()
