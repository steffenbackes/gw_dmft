import os

sigma = 0.002
inf = open('../../s_imp.dat','r')
useorb = 6
w2 = 0
w1= 50

outf = open('gf.dat','w')
simp = []
for line in inf:
	simp.append( [float(x) for x in line.split()] )

s0 =  ( simp[-w2][1+2*useorb]*simp[-w2][0]**2 - simp[-w1][1+2*useorb]*simp[-w1][0]**2 )/( simp[-w2][0]**2 - simp[-w1][0]**2 )
s1 = -( simp[-w2][2+2*useorb]*simp[-w2][0]**3 - simp[-w1][2+2*useorb]*simp[-w1][0]**3 )/( simp[-w2][0]**2 - simp[-w1][0]**2 )

print s0
print s1

os.system("sed -i s/s0=.*/s0="+str(s0)+"/ kk.py")
os.system("sed -i s/s1=.*/s1="+str(s1)+"/ kk.py")

for data in simp:
   outf.write( str(data[0])+ '\t' + str( (data[1+2*useorb]-s0)/s1 ) + '\t'  + str(sigma) + '\t' \
                                  + str( (data[2+2*useorb])/s1 ) + '\t'  + str(sigma) + '\n' )

outf.close()
inf.close()
