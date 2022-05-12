sigma = 0.001
inf = open('g_loc.dat','r')
useorb = 0

outf = open('gf.dat','w')
for line in inf:
   data = [float(x) for x in line.split()]
   outf.write( str(data[0])+ '\t' + str( (data[1+2*useorb]) ) + '\t'  + str(sigma) + '\t' \
                                  + str( (data[2+2*useorb]) ) + '\t'  + str(sigma) + '\n' )

outf.close()
inf.close()
