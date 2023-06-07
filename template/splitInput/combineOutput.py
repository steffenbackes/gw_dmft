import numpy as np
np.set_printoptions(precision=3, suppress=True, linewidth=300)

norbPerAtom = 5
copyTo = [ [0,3], [1,2,4,5] ]


noIneqAtoms = len(copyTo)
norbTotal = sum([len(e) for e in copyTo]) * norbPerAtom

# read number of frequencies
inf = open('gw_input.dat','r')
inf.readline()
inf.readline()
nw = int( inf.readline().split()[1] )
inf.close()
############################

inf = [ open('imp'+str(A)+'/Gwl.dat','r') for A in range(len(copyTo)) ]

outf = open('Gwl.dat','w')
for n in range(nw):
   Gimp = np.zeros((norbTotal,norbTotal,2),dtype=complex)
   w = 0.0
   for A in range(len(copyTo)):
      data = [float(x) for x in inf[A].readline().split()]
      w = data[0]

      for s in range(2):
         for m1 in range(norbPerAtom):
            for m2 in range(norbPerAtom):

               for cpTo in copyTo[A]:
                  i = norbPerAtom*cpTo + m1
                  j = norbPerAtom*cpTo + m2
                  Gimp[i,j,s] =       data[1+s*norbPerAtom**2*2 + m1*norbPerAtom*2 + m2*2 + 0]  # real
                  Gimp[i,j,s] += 1.0j*data[1+s*norbPerAtom**2*2 + m1*norbPerAtom*2 + m2*2 + 1]  # imag
      
   outf.write( "{:.7e}".format(w) + '\t' )
   for s in range(2):
      for m1 in range(norbTotal):
         for m2 in range(norbTotal):
            outf.write( "{:.7e}".format(Gimp[m1,m2,s].real) + '\t' + "{:.7e}".format(Gimp[m1,m2,s].imag) + '\t' )
   outf.write('\n')


for f in inf:
   f.close()
outf.close()
