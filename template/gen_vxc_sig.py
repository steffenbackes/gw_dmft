nkx = 160
nky = 1
nkz = 1
norb = 1
nw = 500

outf = open('gw_sigma.dat','w')
outf2 = open('gw_v_xc.dat','w')
outf.write('# Sxc_wann interpolated to Matsubara freq. \n')
outf2.write('# Vxc in the Wannier representation \n')
ik = 1
for kx in range(nkx):
   for ky in range(nky):
      for kz in range(nkz):


         outf.write( '# ' + str(ik) + '  ' + str(kx*1.0/nkx) + '  ' + str(ky*1.0/nky) + '  ' + str(kz*1.0/nkz) + '   ' + str(0.00000) + '   ' + str(0.0000) + '\n'   )

         for n in range(nw):
            outf.write( str(n) + '  ' )
            for m1 in range(norb):
               for m2 in range(norb):
                  outf.write( str(0.0) + '  ' + str(0.0) + ' ' )
            outf.write( '\n' )

         outf.write( '\n' )
	
			#############

         outf2.write( str(ik) + '  ' + str(kx*1.0/nkx) + '  ' + str(ky*1.0/nky) + '  ' + str(kz*1.0/nkz) + '\t'   )
         for m1 in range(norb):
            for m2 in range(norb):
               outf2.write( str(0.0) + '  ' + str(0.0) + ' ' )
            outf2.write( '\n' )

         ik += 1

outf.close()
