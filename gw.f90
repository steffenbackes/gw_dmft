!!! ============================================================ !!!
!!! = GW module calculates P=GG and Sigma=GW =================== !!!
!!! ============================================================ !!!

module gw
   use constants
   use params
   use mathfunc
   use matsfunc_ops
   implicit none

   contains
!  ============================================================
!  == Constructs the bare interaction matrix V_gw if not read from GW input
!  ============================================================
   subroutine setup_bare_interaction(V)
      complex(kr),intent(inout)    :: V(:,:,:,:,:)      ! norb**2,norb**2,nspin,nkpts,nnu

      integer(ki)            :: m1,m2,a,b, iounit=13, i,j, ikx,iky,ikz, ik,x,y
      real(kr)               :: F2,F4,F0,U0,J1,J2,J3,J4,Jmat(5,5), kx, ky, kz, dist
      real(kr)               :: readUtmp(norb_dmft,norb_dmft),readJtmp(norb_dmft,norb_dmft)

      write(*,'(A)') ' Set up the bare interaction...'
write(*,'(A)') '!! WARNING !!! We use dmft_orbitals in setup_bare_interaction(V) instead of norbPerAtom !!!'
write(*,'(A)') '!! WARNING !!! All the other GW routines work on the full norb !!!'

      V = (0.0_kr,0.0_kr)

!      if usePolK==0 then  ! No PolK: V_gw only has one k-point and frequency, so we either read from file or use Uinput,Jinput
      ! For now always set up the matrix yourself, since reading V_gw from file is not yet implemented

         if (readUmatrix==1) then ! Use data from umatrix_in.dat
            write(*,'(A)') ' Reading interaction matrix from umatrix_in.dat ...'

            open(unit=iounit,file="umatrix_in.dat",status="old")
            ! read Umat
            do m1=1,norb_dmft
               read(iounit,*) readUtmp(m1,:)
            enddo
            ! read Jmat
            do m1=1,norb_dmft
               read(iounit,*) readJtmp(m1,:)
            enddo

            ! copy all U   values to U_iiii 
            ! copy all U' values to U_ijij = U_(ii)(jj) in the two-particle basis matrix representation
            ! copy all J values to U_iijj = U_(ij)(ij) and U_ijji  
            do i=1,norb_dmft
               do j=1,norb_dmft
                  ! U and U' term
                  m1 = dmftOrbsIndex(i)     
                  m2 = dmftOrbsIndex(j)     
                  a = (m1-1)*norb+m1          
                  b = (m2-1)*norb+m2          
                  V(a,b,:,:,:) = readUtmp(i,j)*U0scaleFac
      
                  if (i/=j) then
                     ! first J term iijj
                     a = (m1-1)*norb+m2          
                     b = (m1-1)*norb+m2          
                     V(a,b,:,:,:) = readJtmp(i,j)
                     ! second J term ijji
                     a = (m1-1)*norb+m2          
                     b = (m2-1)*norb+m1          
                     V(a,b,:,:,:) = readJtmp(i,j)
                  endif
      
               enddo
            enddo

            close(iounit)

         else ! we construct the Umatrix ourselves
             ! For one orbital
             if (norb_dmft==1) then
                V(1,1,:,:,:) = Uinput*U0scaleFac
   
             ! For three orbitals use normal Kanamori
             else if (norb_dmft .lt. 5) then
                 do i=1,norb_dmft
                    ! U term
                    m1 = dmftOrbsIndex(i)     
                    a = (m1-1)*norb+m1          
                    b = (m1-1)*norb+m1          
                    V(a,b,:,:,:) = Uinput*U0scaleFac
     
                    do j=1,norb_dmft
                       if (i/=j) then
                          ! U' term
                          m2 = dmftOrbsIndex(j)     
                          a = (m1-1)*norb+m1          
                          b = (m2-1)*norb+m2          
                          V(a,b,:,:,:) = (Uinput-2*Jinput)*U0scaleFac
        
                          ! first J term iijj
                          a = (m1-1)*norb+m2          
                          b = (m1-1)*norb+m2          
                          V(a,b,:,:,:) = Jinput
                          ! second J term ijji
                          a = (m1-1)*norb+m2          
                          b = (m2-1)*norb+m1          
                          V(a,b,:,:,:) = Jinput
                      endif
               enddo
            enddo
             ! For five orbitals use Slater form!
             ! Assume ordering dz2, dx2y2, dxy, dxz, dyz like Wien2K
             else if (norb_dmft==5) then
                F0 = Uinput
                F2 = Jinput*8*14/13.0_kr  ! use F4/F2 = 5/8
                F4 = 14*Jinput - F2       ! follow thesis of Steffen Backes, page 97
   
                U0 =  Uinput            + Jinput*8/7.0_kr
                J1 =  F2*3/49.0_kr      + F4*20/(9*49.0_kr)
                J2 = -Jinput*10/7.0_kr + 3*J1
                J3 =  Jinput*30/7.0_kr - 5*J1
                J4 =  Jinput*20/7.0_kr - 3*J1
   
                ! For Jmat now the orbital order comes into play! See above!
                ! If orbital order is different, this must be modified!
                Jmat = reshape( [ zero,   J2,   J2,   J4,   J4,  &
                              &     J2, zero,   J3,   J1,   J1,  &
                              &     J2,   J3, zero,   J1,   J1,  &
                              &     J4,   J1,   J1, zero,   J1,  &
                              &     J4,   J1,   J1,   J1, zero ] ,[5,5] )
   
                 do i=1,norb_dmft
                    ! U term
                    m1 = dmftOrbsIndex(i)     
                    a = (m1-1)*norb+m1          
                    b = (m1-1)*norb+m1          
                    V(a,b,:,:,:) = U0*U0scaleFac
     
                    do j=1,norb_dmft
                       if ( i/=j) then
                          ! U' term
                          m2 = dmftOrbsIndex(j)     
                          a = (m1-1)*norb+m1          
                          b = (m2-1)*norb+m2          
                          V(a,b,:,:,:) = (U0-2*Jmat(i,j))*U0scaleFac
        
                          ! first J term iijj
                          a = (m1-1)*norb+m2          
                          b = (m1-1)*norb+m2          
                          V(a,b,:,:,:) = Jmat(i,j)
                          ! second J term ijji
                          a = (m1-1)*norb+m2          
                          b = (m2-1)*norb+m1          
                          V(a,b,:,:,:) = Jmat(i,j)
                     endif
                  enddo
               enddo

               write(*,'(A)') 'Slater 5-orb Umat:'
               do i=1,norb_dmft
                  do j=1,norb_dmft
                     m1 = dmftOrbsIndex(i)     
                     m2 = dmftOrbsIndex(j)     
                     a = (m1-1)*norb+m1          
                     b = (m2-1)*norb+m2          
                     write(*,'(F9.5, A)',advance='no') real(V(a,b,1,1,1)), ' '
                  enddo
                  write(*,'(A)') ' '
               enddo
               write(*,'(A)') 'Slater 5-orb Jmat:'
               do i=1,norb_dmft
                  do j=1,norb_dmft
                     m1 = dmftOrbsIndex(i)     
                     m2 = dmftOrbsIndex(j)     
                     a = (m1-1)*norb+m2          
                     b = (m1-1)*norb+m2          
                     write(*,'(F9.5, A)',advance='no') real(V(a,b,1,1,1)), ' '
                  enddo
                  write(*,'(A)') ' '
               enddo
               write(*,'(A)') ' '

             else
                write(*,'(A)') 'ERROR: Interaction matrix of #orbitals != 1,3,5 for each atom is not supported!'
                stop
             endif ! norb_dmft
         endif ! read Umatrix
!      else
         ! We have usePolK==1

      ! Now add a long-range interaction (if needed)
      if ( size(V(1,1,1,:,1))==nkpts ) then
         do x=-dist_interaction,dist_interaction
            do y=-dist_interaction,dist_interaction
               dist = sqrt(1.0*x*x+y*y)
               if ( ( x .ne. 0 .or. y .ne. 0) .and. dist .le. dist_interaction ) then
                  write(*,'(A,F6.3,A,I3,A,I3,A)') 'Adding interaction V=',Unninput/dist,' for site at (',x,',',y,')'

                  do ikx=0,nkx-1
                  do iky=0,nky-1
                  do ikz=0,nkz-1
                     kx = ikx*2*pi/nkx
                     ky = iky*2*pi/nky
                     kz = ikz*2*pi/nkz
                     ik = ikx*nky*nkz + iky*nkz + ikz + 1
            
                     do i=1,norb_dmft
                        do j=1,norb_dmft
                           m1 = dmftOrbsIndex(i)     
                           m2 = dmftOrbsIndex(j)
                           a = (m1-1)*norb+m1
                           b = (m2-1)*norb+m2
         
                           !V(a,b,:,ik,:) = V(a,b,:,ik,:) + 2*Unninput*( cos(kx)+cos(ky)+0*cos(kz) )
                           V(a,b,:,ik,:) = V(a,b,:,ik,:) + Unninput * exp(ci*(kx*x + ky*y)) / dist
                        enddo ! j orb
                     enddo ! i orb
                  enddo ! kz    
                  enddo ! ky
                  enddo ! kx

               endif ! dist < cutoff
            enddo ! y
         enddo ! x
      endif

!      endif ! usePolK
   end subroutine setup_bare_interaction

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculates the polarization from GG !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine polarization_gg(pol,gf_tau,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      complex(kr),intent(inout) :: pol(:,:,:,:,:)      ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(inout) :: gf_tau(:,:,:,:,:)   ! norb,norb,nspin,nkpts,ntau+1
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT

      !!!!!!!
      integer(ki)              :: iqx,iqy,iqz,ipx,ipy,ipz, ip,iq, ipminq, n,m, s, i,j,k,l, a,b,it
      complex(kr)              :: gtmp(norb,norb), ptmp(ntau+1), gtmp2(norb,norb), coeffs(5)
      real(kr)                 :: tau, c2(norb,norb), wf, ident(norb,norb)
      
      if (useSigK==2) then

      pol = (0.0_kr,0.0_kr)
      gf_tau = (0.0_kr,0.0_kr)
      ident = (0.0_kr)
      ptmp = (0.0_kr,0.0_kr)
      do a=1,norb
         ident(a,a) = 1.0_kr
      enddo
      !!!!!!!

      write(*,'(A)') ' Filling the G(tau) array for the GG calculation...'
      do iqx=0,nkx-1
         call progress(iqx,nkx)
      do iqy=0,nky-1
      do iqz=0,nkz-1
         iq = iqx*nky*nkz + iqy*nkz + iqz +1
      
         do s=1,nspin
            ! first get the c2 coeff at the largest matsubara freq
            gtmp = gfk(s,iq,nomega,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
            c2 = -real(gtmp)*wn(nomega-1)**2

            do n=1,nomega
               wf = wn(n-1)
               gtmp = gfk(s,iq,n,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
               do it=1,ntau+1
                  tau = (it-1)*beta/ntau
                  gf_tau(:,:,s,iq,it) = gf_tau(:,:,s,iq,it) + 2*(real(gtmp)*cos(wf*tau) + aimag(gtmp)*sin(wf*tau)   &
                                   &   + ident*sin(wf*tau)/wf + c2*cos(wf*tau)/wf**2 )/beta
               enddo ! it
            enddo ! n
            
            ! now add the analytic terms
            do it=1,ntau+1
               tau = (it-1)*beta/ntau
               gf_tau(:,:,s,iq,it) = gf_tau(:,:,s,iq,it) - ident*0.5 + c2*(2*tau-beta)/4

!write(*,*) real(gf_tau(1,1,s,iq,it))
            enddo
!stop 1

         enddo ! spin s
      enddo ! iqz
      enddo ! iqy
      enddo ! iqx

      ! Now calculate GG
      write(*,'(A)') ' Calculating the lattice polarization from GG...'
      write(*,'(A)') ' WARNING: Currently only for paramagnetic case!'
      s=1
         do i=1,norb
         do j=1,norb
         do k=1,norb
         do l=1,norb

            ! in the end we save it in pol(a,b)
            a = (i-1)*norb + k 
            b = (j-1)*norb + l 

            ! only calculate the upper diagonal and then copy
            if (a .gt. b) then
               continue
            endif

            do iqx=0,nkx-1
               call progress(iqx,nkx)
            do iqy=0,nky-1
            do iqz=0,nkz-1
               !iq =  i*norb**3*nkpts + j*norb**2*nkpts + k*norb*nkpts + l*nkpts & 
               !   & + iqx*nky*nkz + iqy*nkz + iqz +1
               iq = iqx*nky*nkz + iqy*nkz + iqz +1

               ptmp = (0.0_kr,0.0_kr)               

               do ipx=0,nkx-1
               do ipy=0,nky-1
               do ipz=0,nkz-1
                  ip = ipx*nky*nkz + ipy*nkz + ipz +1

                  ! modulo gives us directly the right rewinding for the k coordinates
                  ipminq = modulo(ipx-iqx,nkx)*nky*nkz + modulo(ipy-iqy,nky)*nkz + modulo(ipz-iqz,nkz) +1

                  do it=1,ntau+1
                     ptmp(it) = ptmp(it) + (-gf_tau(l,i,s,ipminq,ntau-it+2)) * gf_tau(k,j,s,ip,it) 
                  enddo


               enddo ! ipz
               enddo ! ipy
               enddo ! ipx
               ptmp = ptmp/nkpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do it=1,ntau+1
!  write(*,*) real(ptmp(it) )
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   do it=1,ntau+1
!      ptmp(it) = ptmp(it) - ptmp(ntau/2)
!   enddo
!   write(*,*) iq,ptmp(ntau/2) ,staticpart

               ! now we have Pol(tau,q), just do fourier trafo
               do n=1,nnu
                  wf = vn(n-1)
!                  do it=1,ntau+1
!                     tau = (it-1)*beta/ntau
!                     pol(:,:,s,iq,n) = pol(:,:,s,iq,n) + ptmp(:,:,it)*exp(ci*wf*tau)
!                  enddo

                   !simpson 
                  pol(a,b,s,iq,n) = pol(a,b,s,iq,n) + ptmp(1) + ptmp(ntau+1)
                  do it=1,ntau/2-1
                     tau = (2*it-0)*beta/ntau
                     pol(a,b,s,iq,n) = pol(a,b,s,iq,n) + 2*ptmp(2*it+1)*exp(ci*wf*tau)
                  enddo
                  do it=1,ntau/2
                     tau = (2*it-1)*beta/ntau
                     pol(a,b,s,iq,n) = pol(a,b,s,iq,n) + 4*ptmp(2*it)*exp(ci*wf*tau)
                  enddo
               
               enddo ! n
               pol(a,b,s,iq,:) = pol(a,b,s,iq,:) * beta/(3*ntau)  ! 3 from Simpson rule

               ! interpolate to zero
!               coeffs = get_lowfreq_coeff_bos(pol(a,b,s,iq,2:nnu),5,6)
!               pol(a,b,s,iq,1) = coeffs(1)

!write(*,*) iq
!do n=1,10
!   write(*,*) real(pol(a,b,s,iq,n))
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! now we have pol at all freqs, do fourier
!write(*,*) "BACKTRAFO of P:"
!            ptmp = (0.0_kr)
!            c2 = -real(pol(:,:,s,iq,nnu))*vn(nnu-1)**2
!            do it=1,ntau+1
!               tau = (it-1)*beta/ntau
!               ! w=0 part
!               ptmp(it) = ptmp(it) + pol(a,b,s,iq,1)/beta
!               do n=2,nnu
!                  wf = vn(n-1)
!                  ptmp(it) = ptmp(it) + 2*(real(pol(a,b,s,iq,n))*cos(wf*tau)  & 
!                                      & + aimag(pol(a,b,s,iq,n))*sin(wf*tau)  &
!                                      & + c2(a,b)*cos(wf*tau)/wf**2 )/beta
!               enddo ! n
!               ! now add the analytic terms
!               ptmp(it) = ptmp(it) - c2(a,b)*(tau**2 - tau*beta + beta*beta/6)/(2*beta)
!   write(*,*) real(ptmp(it))
!            enddo ! it
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

            enddo ! iqz
            enddo ! iqy
            enddo ! iqx

            pol(b,a,s,:,:) = pol(a,b,s,:,:)  ! copy to lower diagonal
         enddo ! l
         enddo ! k
         enddo ! j
         enddo ! i

      pol(:,:,2,:,:) = pol(:,:,1,:,:)
      pol = 2*pol  ! twice for sum over spin

      endif ! useSigK==2

   end subroutine polarization_gg


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculates W from P and V !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_W(W,P,pimp,p_dc,V)
      complex(kr),intent(inout) :: W(:,:,:,:,:)      ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: P(:,:,:,:,:)      ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: pimp(:,:,:,:)     ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: p_dc(:,:,:,:)     ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: V(:,:,:,:,:)      ! norb**2,norb**2,nspin,nkpts,nnu

      integer(ki) :: q,n,s
      complex(kr) :: tmp(norb**2,norb**2), unity(norb**2,norb**2), pol(norb**2,norb**2)

      if (usePolK/=0 .or. useSigK==2) then

         unity=(0.0_kr,0.0_kr)
         do n=1,norb**2
            unity(n,n) = 1.0_kr
         enddo
   
         do s=1,nspin
            do q=1,size(W(1,1,1,:,1))
               do n=1,size(W(1,1,1,1,:))

                  if (usePolK==2) then ! we use Pol impurity
                     pol = P(:,:,s,q,n) + pimp(:,:,s,n) - p_dc(:,:,s,n)
                  else       ! No update of Polarization with Pimp, just use GW
                     pol = P(:,:,s,q,n)
                  endif

                  tmp = unity - matmul( pol, V(:,:,s,q,n) )
                  call get_inverse(tmp, tmp, norb**2)
                  W(:,:,s,q,n) = matmul(V(:,:,s,q,n),tmp)
               enddo
            enddo
         enddo
         write(*,*) "V at M=",V(1,1,1,(nkx/2)*nky*nkz+(nky/2)*nkz+1,1)
         write(*,*) "W at M=",W(1,1,1,(nkx/2)*nky*nkz+(nky/2)*nkz+1,1)

      else
         W = V
      endif

   end subroutine get_W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculates the selfenergy from GW !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine sigma_gw(sigma,W,gf_tau)
      complex(kr),intent(inout) :: sigma(:,:,:,:,:)    ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: W(:,:,:,:,:)        ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: gf_tau(:,:,:,:,:)   ! norb,norb,nspin,nkpts,ntau+1

      !!!!!!!
      integer(ki)              :: iqx,iqy,iqz,ipx,ipy,ipz, ip,iq, ipminq, n,m, s, i,j,k,l, a,b,it
      complex(kr),allocatable  :: w_tau(:,:,:,:,:)         ! norb**2,norb**2,nspin,nkpts,ntau
      complex(kr)              :: stmp(ntau+1), staticpart         
      real(kr)                 :: tau, c2(norb**2,norb**2), wf, coeffs(2)
      real(kr),allocatable     :: Vqbare(:,:,:)  !norb**2,norb**2,nkpts

      if (useSigK==2) then

      allocate( w_tau(norb**2,norb**2,nspin,nkpts,ntau+1) )
      allocate( Vqbare(norb**2,norb**2,nkpts) )

      sigma = (0.0_kr,0.0_kr)
      w_tau = (0.0_kr,0.0_kr)
      stmp = (0.0_kr,0.0_kr)
      !!!!!!!

      write(*,'(A)') ' Filling the W(tau) array for the GW calculation...'
      do iqx=0,nkx-1
         call progress(iqx,nkx)
      do iqy=0,nky-1
      do iqz=0,nkz-1
         iq = iqx*nky*nkz + iqy*nkz + iqz +1
      
         do s=1,nspin
            ! get the high freq coefficients to take out the bare constant part
            do i=1,norb**2
               do j=1,norb**2
                  coeffs = get_highfreq_coeff_bosonic( W(i,j,s,iq,:) )
                  Vqbare(i,j,iq) = coeffs(1)
                  c2(i,j) = coeffs(2)
               enddo
            enddo
!write(*,*) c2

            do it=1,ntau+1
               tau = (it-1)*beta/ntau
               ! w=0 part
               w_tau(:,:,s,iq,it) = w_tau(:,:,s,iq,it) + (W(:,:,s,iq,1)-Vqbare(:,:,iq))/beta
               do n=2,nnu
                  wf = vn(n-1)
                  w_tau(:,:,s,iq,it) = w_tau(:,:,s,iq,it) + 2*(real(W(:,:,s,iq,n)-Vqbare(:,:,iq))*cos(wf*tau)  & 
                                                          & + aimag(W(:,:,s,iq,n)               )*sin(wf*tau)  &
                                                          & + c2*cos(wf*tau)/wf**2 )/beta
               enddo ! n
               ! now add the analytic terms
               w_tau(:,:,s,iq,it) = w_tau(:,:,s,iq,it) - c2*(tau**2 - tau*beta + beta*beta/6)/(2*beta)
!   write(*,*) real(gf_tau(1,1,s,iq,it))
            enddo ! it
!stop 1

         enddo ! spin s
      enddo ! iqz
      enddo ! iqy
      enddo ! iqx

      ! Now calculate GW
      write(*,'(A)') ' Calculating the Selfenergy from GW...'
      write(*,'(A)') ' WARNING: Currently only for paramagnetic case!'
      s=1
         do a=1,norb
         do b=1,norb
            ! only calculate the upper diagonal and then copy
            if (a .gt. b) then
               continue
            endif

            do iqx=0,nkx-1
               call progress(iqx,nkx)
            do iqy=0,nky-1
            do iqz=0,nkz-1
               !iq =  a*norb*nkpts + b*nkpts & 
               !   & + iqx*nky*nkz + iqy*nkz + iqz +1
               iq = iqx*nky*nkz + iqy*nkz + iqz +1

               stmp = (0.0_kr,0.0_kr)               

               do i=1,norb
               do j=1,norb
                  do ipx=0,nkx-1
                  do ipy=0,nky-1
                  do ipz=0,nkz-1
                     ip = ipx*nky*nkz + ipy*nkz + ipz +1
   
                     ! modulo gives us directly the right rewinding for the k coordinates
                     ipminq = modulo(iqx-ipx,nkx)*nky*nkz + modulo(iqy-ipy,nky)*nkz + modulo(iqz-ipz,nkz) +1
   
                     do it=1,ntau+1
                        stmp(it) = stmp(it) - w_tau((a-1)*norb+i,(j-1)*norb+b,s,ip,it) * gf_tau(i,j,s,ipminq,it) 
                     enddo
   
                  enddo ! ipz
                  enddo ! ipy
                  enddo ! ipx
               enddo ! j
               enddo ! i
               stmp = stmp/nkpts

!do it=1,ntau+1
!   write(*,*) real(stmp(it))
!enddo

               ! now we have sigma(tau,q), just do fourier trafo
               do n=1,nomega
                  wf = wn(n-1)
!                  do it=1,ntau+1
!                     tau = (it-1)*beta/ntau
!                     pol(:,:,s,iq,n) = pol(:,:,s,iq,n) + ptmp(:,:,it)*exp(ci*wf*tau)
!                  enddo

                   !simpson 
                  sigma(a,b,s,iq,n) = sigma(a,b,s,iq,n) + stmp(1) + stmp(ntau+1)
                  do it=1,ntau/2-1
                     tau = (2*it-0)*beta/ntau
                     sigma(a,b,s,iq,n) = sigma(a,b,s,iq,n) + 2*stmp(2*it+1)*exp(ci*wf*tau)
                  enddo
                  do it=1,ntau/2
                     tau = (2*it-1)*beta/ntau
                     sigma(a,b,s,iq,n) = sigma(a,b,s,iq,n) + 4*stmp(2*it)*exp(ci*wf*tau)
                  enddo
               
               enddo ! n
               sigma(a,b,s,iq,:) = sigma(a,b,s,iq,:) * beta/(3*ntau)  ! 3 from Simpson rule

!write(*,*) iq
!do n=1,10
!   write(*,*) real(pol(a,b,s,iq,n))
!enddo
               ! Now add the Hartree-Fock part from Ubare
               stmp = (0.0_kr,0.0_kr)
               do i=1,norb
               do j=1,norb
                  do ipx=0,nkx-1
                  do ipy=0,nky-1
                  do ipz=0,nkz-1
                     ip = ipx*nky*nkz + ipy*nkz + ipz +1
                     ipminq = modulo(iqx-ipx,nkx)*nky*nkz + modulo(iqy-ipy,nky)*nkz + modulo(iqz-ipz,nkz) +1
   
                     ! Hartree
                     stmp(1) = stmp(1) - Vqbare((a-1)*norb+b,(i-1)*norb+j,1     ) * gf_tau(i,j,1,ip,ntau+1)
                     stmp(1) = stmp(1) - Vqbare((a-1)*norb+b,(i-1)*norb+j,1     ) * gf_tau(i,j,2,ip,ntau+1)
      
                     ! Fock
                     stmp(1) = stmp(1) + Vqbare((a-1)*norb+i,(j-1)*norb+b,ipminq) * gf_tau(i,j,s,ip,ntau+1)
                  enddo ! ipz
                  enddo ! ipy
                  enddo ! ipx
               enddo ! j
               enddo ! i
               sigma(a,b,s,iq,:) = sigma(a,b,s,iq,:) + stmp(1)/nkpts

            enddo ! iqz
            enddo ! iqy
            enddo ! iqx

            sigma(b,a,s,:,:) = sigma(a,b,s,:,:)  ! copy to lower diagonal
         enddo ! b
         enddo ! a

      sigma(:,:,2,:,:) = sigma(:,:,1,:,:)

      deallocate( w_tau )
      deallocate( Vqbare )

      endif ! (useSigK==2) then
   end subroutine sigma_gw

end module gw
