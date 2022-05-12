import numpy as np
np.set_printoptions(precision=3)

xorb=1
yorb=1
norb=xorb*yorb
U=4.0
t = 0.25 # NN hopping
Nk = 64

nw = 51
wmin = -5.0
wmax = 5.0
dw = (wmax-wmin)/nw

outfh = open('gw_h_dft.dat','w')
outfh.write('#Hamiltonian in the local diagonal basis \n')
ii = 1
spec = np.zeros((norb,nw))
w = np.linspace(wmin,wmax,nw)
for ikx in range(Nk):
	print ikx*100.0/Nk,'% done'

	for iky in range(Nk):

		kx = ikx*2*np.pi/Nk 
		ky = iky*2*np.pi/Nk 

		outfh.write( str(ii) + '\t' + str(kx) + '\t' + str(ky) + '\t' + str(0.0) + '\t'  )
		ii += 1

		Hk = np.zeros((norb,norb),dtype=complex)
#		for x in range(xorb):
#			for y in range(yorb):
#				Hk[x*yorb+y,(x+1)%xorb*yorb+y] += -t*( 1*(x+1<xorb) + np.exp(-1.0j*kx)*(x+1>=xorb) )
#				Hk[x*yorb+y,(x-1)%xorb*yorb+y] += -t*( 1*(x-1>=0) + np.exp(+1.0j*kx)*(x-1<0) )
#				Hk[x*yorb+y,x%xorb*yorb+(y+1)%xorb] += -t*( 1*(y+1<yorb) + np.exp(-1.0j*ky)*(y+1>=yorb) )
#				Hk[x*yorb+y,x%xorb*yorb+(y-1)%xorb] += -t*( 1*(y-1>=0) + np.exp(+1.0j*ky)*(y-1<0) )

		Hk[0,0] = -2*t*( np.cos(kx) + np.cos(ky) )


		# 2x2
		#Hk = -t*np.array([[ 0.0, 1.0 + np.exp(-1.0j*kx) , 1.0 + np.exp(+1.0j*ky) , 0.0 ],
      #                  [ 0.0, 0.0 , 0.0 ,                    1.0 + np.exp(+1.0j*ky) ],
      #                  [ 0.0, 0.0 , 0.0 ,                    1.0 + np.exp(-1.0j*kx) ],
      #                  [ 0.0, 0.0 , 0.0 ,                    0.0                    ] ])


		for i in range(norb):
			for j in range(norb):
				outfh.write( str(Hk[j,i].real/27.211386) + '\t' + str(Hk[j,i].imag/27.211386) + '\t'  )
		outfh.write( '\n'  )

		for n in range(nw):
			spec[:,n] += np.diag(  -( np.linalg.inv( (w[n] + 0.02j)*np.identity(norb) - Hk ) ).imag/np.pi ) / Nk

outfh.close()

outfs = open('spec_h_dft.dat','w')
for n in range(nw):
	outfs.write( str(w[n]) +'\t' )
	for j in range(norb):
		outfs.write( str(spec[j,n]) + '\t' )
	outfs.write('\n')
outfs.close()

outf = open("umatrix_in.dat",'w')
for i in range(norb):
	for j in range(norb):
		Uval=0.0
		if (i==j):
			Uval = U
		outf.write(str(Uval)+' ')
	outf.write('\n')
outf.write('\n')
for i in range(norb):
	for j in range(norb):
		outf.write(str(0.0)+' ')
	outf.write('\n')
outf.write('\n')
outf.close()
