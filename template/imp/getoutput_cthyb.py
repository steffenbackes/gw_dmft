import numpy as np

#
inf = open('input.ini','r')
for line in inf:
	if ( 'model.beta' in line ):
		beta = float( line.replace('=',' ').split()[1] )
inf.close()

# read GF
inf = open('gf.dat','r')
outf =open('Gwl.dat','w')
for i in range(3):
	inf.readline()

data = inf.readline().replace(',',' ').split()
nw = int(data[4])
norb = int(data[5])/2

print 'Getoutput:'
print 'Found beta=',beta
print 'Found nw=',nw
print 'Found norb=',norb

inf.readline()

for n in range(nw):
	outf.write( "{:.7e}".format( (2*n+1)*np.pi/beta ) + '\t')

	gimp = np.zeros((norb,norb,2),dtype=complex)
	gimps = np.zeros((norb,norb,2,2),dtype=complex)

	for m1 in range(norb):
		for s1 in range(2):
			for m2 in range(norb):
				for s2 in range(2):
					data = [float(x) for x in inf.readline().replace(',',' ').split()[4:] ] 
					gimps[m1,m2,s1,s2] = data[0] + data[1]*1.0j

	for s in range(2):
		for m1 in range(norb):
			for m2 in range(norb):
				outf.write("{:.7e}".format(gimps[m1,m2,s,s].real) + '\t' + "{:.7e}".format(gimps[m1,m2,s,s].imag) + '\t') 

	outf.write('\n')
outf.close()
inf.close()

# legendre coeffs
inf = open('gf_legendre.dat','r')
for i in range(3):
	inf.readline()

data = inf.readline().replace(',',' ').split()
nlegendre = int(data[6])
print 'Found nlegendre=',nlegendre

legendre = np.zeros((nlegendre,norb,norb,2))
inf.readline()

for m1 in range(norb):
	for s1 in range(2):
		for m2 in range(norb):
			for s2 in range(2):
				for l in range(nlegendre):
					data = [float(x) for x in inf.readline().replace(',',' ').split()[4:] ] 

					if (s1==s2):
						legendre[l,m1,m2,s1] = data[0]
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
