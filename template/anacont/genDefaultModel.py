import numpy as np

nfreq  = 2000
wstart = -41.0
wend   = 41.0

gauss_mu  = [-30.0, -15.0, 1.0, 15.0, 30.0]
gauss_sig = [  2.0,   2.0, 1.0,  2.0,  2.0]
gauss_weig= [  1.0,   2.0, 10.0,  2.0,  1.0]

norm = sum(gauss_weig)

def gauss(x,mu,sig):
	return np.exp(-0.5*(x-mu)**2/sig**2) / (sig*np.sqrt(2*np.pi) )

outf = open('defaultmodel.dat','w')
dw = (wend-wstart)/nfreq
norm2 = 0.0
for n in range(nfreq+1):
	w = wstart + n*dw
	val = 0.0
	for g in range(len(gauss_mu)):
		val += gauss(w,gauss_mu[g],gauss_sig[g]) * gauss_weig[g]

	#outf.write(str(val/norm) + '\n')
	outf.write(str(w) + '\t' + str(val/norm) + '\n')
	norm2 += dw*val/norm
outf.close()
print 'Final norm: ',norm2
