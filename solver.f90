!  ============================================================
!  == Impurity solver interface
!  ============================================================

module solver
   use constants
   use params
   implicit none
   
 contains

!  ============================================================
!  == Write the input for the CTHYB: umatrix, local mu, K(t),K'(t) functions
!  == We will only write the data for the first equiv. atom
!  ============================================================   
   subroutine write_solver_input(uloc,hybrid,gbath,mu_loc)
      use matsfunc_ops
      use io_ops
      use mathfunc
      complex(kr),intent(in) :: uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in) :: hybrid(:,:,:,:)      ! norb,norb,nspin,nomega
      complex(kr),intent(in) :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr),intent(in)    :: mu_loc(:,:,:)         ! norb,nspin

      integer(ki),parameter    :: umatrix_out = 10, mu_out = 11, k_out=12, mu_matrix_out=13, mu_matrix_solver_out=14
      integer(ki)              :: m1,m2,m3,m4,s1,s2,it,n , a,b,c,d, i,j, atom
      real(kr)                 :: K,Kp,tau, hfcoeff(2), Uval
      complex(kr),allocatable  :: Utensortrafo(:,:,:)         ! norbPerAtom**2,norbPerAtom**2,nspin)
      complex(kr),allocatable  :: UtensortrafoT(:,:,:)        ! norbPerAtom**2,norbPerAtom**2,nspin)
      complex(kr),allocatable  :: mu_loc_trafo(:,:,:)         ! norbPerAtom,norbPerAtom,nspin
      complex(kr),allocatable  :: gbath_trafo(:,:,:,:)        ! norbPerAtom,norbPerAtom,nspin,nomega
      complex(kr),allocatable  :: hybrid_trafo(:,:,:,:)       ! norbPerAtom,norbPerAtom,nspin,nomega
      complex(kr),allocatable  :: Uijkl(:,:)                  ! norbPerAtom**2,norbPerAtom**2                 
      complex(kr),allocatable  :: Uijkl_trafo(:,:,:,:)        ! norbPerAtom**2,norbPerAtom**2 ,nspin,nspin
      complex(kr),allocatable     :: UtrafoA(:,:,:)              ! norbPerAtom,norbPerAtom ,nspin,  extract from Utrafo
      complex(kr),allocatable     :: UtrafoTA(:,:,:)             ! norbPerAtom,norbPerAtom ,nspin,  extract from Utrafo
      complex(kr),allocatable  :: tmpmatrix(:,:)              ! norbPerAtom,norbPerAtom
!      complex(kr)              :: Uavg(noAtoms,nnu)
      real(kr)                 :: coeffs(2)
      character(len=1024)      :: path

      ! We do everything separately for each atom


      ! Read basis transformation matrix, init with identity or read from file
      call read_utrafo()


!      if (diagBath==1) then
!         do s1=1,nspin
!            write(*,'(A,I2)') 'Basis transformation for spin=',s1
!            call print_matrix(Utrafo(:,:,s1),norbPerAtom)
!         enddo
!      endif

      ! ==========================================================
      ! Write Kt.dat for frequency dependent interactions, atoms are handled
      ! within
      call write_Kt_solver(uloc,"Kt.dat")
      ! ==========================================================

      do atom=1,noAtoms    
         write(*,'(A,I2)') 'Write impurity solver input for atom ',atom
         allocate( mu_loc_trafo(norbPerAtom(atom),norbPerAtom(atom),nspin) )                  ; mu_loc_trafo=(0.0_kr,0.0_kr)
         allocate( hybrid_trafo(norbPerAtom(atom),norbPerAtom(atom),nspin,nomega) )           ; hybrid_trafo=(0.0_kr,0.0_kr)
         allocate( Uijkl(norbPerAtom(atom)**2,norbPerAtom(atom)**2) )                         ; Uijkl=(0.0_kr,0.0_kr)
         allocate( Uijkl_trafo(norbPerAtom(atom)**2,norbPerAtom(atom)**2, nspin,nspin) )      ; Uijkl_trafo=(0.0_kr,0.0_kr)
         allocate( Utensortrafo(norbPerAtom(atom)**2,norbPerAtom(atom)**2,nspin) )            ; Utensortrafo=(0.0_kr,0.0_kr)
         allocate( UtensortrafoT(norbPerAtom(atom)**2,norbPerAtom(atom)**2,nspin) )           ; UtensortrafoT=(0.0_kr,0.0_kr)
         allocate( UtrafoA(norbPerAtom(atom),norbPerAtom(atom),nspin) )                       ; UtrafoA=(0.0_kr,0.0_kr)
         allocate( UtrafoTA(norbPerAtom(atom),norbPerAtom(atom),nspin) )                      ; UtrafoTA=(0.0_kr,0.0_kr)
         allocate( tmpmatrix(norbPerAtom(atom),norbPerAtom(atom)) )                           ; tmpmatrix=(0.0_kr,0.0_kr)
  
         ! extract the Utrafo matrix for the given atom
         do s1=1,nspin
            call get_dmftpart_atom(UtrafoA(:,:,s1),Utrafo(:,:,s1),atom)
            call get_dmftpart_atom(UtrafoTA(:,:,s1),UtrafoT(:,:,s1),atom)
         enddo
         if (diagBath==1) then
            do s1=1,nspin
               write(*,'(A,I2,A,I2)') 'Basis transformation for atom= ',atom,', spin=',s1
               call print_matrix(UtrafoA(:,:,s1),norbPerAtom(atom))
            enddo
         endif

         ! =======================================================
         ! Transform interaction matrix ===========================
         ! =======================================================
          do a=0,norbPerAtom(atom)-1
             do b=0,norbPerAtom(atom)-1
   
                do c=0,norbPerAtom(atom)-1
                   do d=0,norbPerAtom(atom)-1
   
                     ! Here we use proper ijkl ordering
                     m1 = norb_dmft(atom,a+1)-1
                     m2 = norb_dmft(atom,b+1)-1
                     m3 = norb_dmft(atom,c+1)-1
                     m4 = norb_dmft(atom,d+1)-1
                     if ( size(uloc(1,1,1,:)) > nfit ) then
                        coeffs = get_highfreq_coeff_bosonic( uloc(m1*norb+m3+1,m2*norb+m4+1,1,:) )
                        !Uijkl(a*norbPerAtom(atom)+b+1,c*norbPerAtom(atom)+d+1) = uloc(m1*norb+m3+1,m2*norb+m4+1,1,size(uloc(1,1,1,:)))
                        Uijkl(a*norbPerAtom(atom)+b+1,c*norbPerAtom(atom)+d+1) = coeffs(1)
                     else
                        Uijkl(a*norbPerAtom(atom)+b+1,c*norbPerAtom(atom)+d+1) = uloc(m1*norb+m3+1,m2*norb+m4+1,1,1)
                     endif
                  enddo
               enddo
   
   !             ! U and U' values
   !             m1 = a*norbPerAtom(atom)+b +1
   !             m2 = a*norbPerAtom(atom)+b +1
   !             Uijkl(m1,m2) = umatrix(a+1,1,b+1,2)
   !
   !             ! J-abba values
   !             m1 = a*norbPerAtom(atom)+b +1 
   !             m2 = b*norbPerAtom(atom)+a +1
   !             Uijkl(m1,m2) = umatrix(a+1,1,b+1,2) - umatrix(a+1,1,b+1,1)
   !
   !             ! J-aabb values
   !             m1 = a*norbPerAtom(atom)+a +1
   !             m2 = b*norbPerAtom(atom)+b +1
   !             Uijkl(m1,m2) = umatrix(a+1,1,b+1,2) - umatrix(a+1,1,b+1,1)
   !
             enddo
          enddo
   
         if (diagBath==1) then
            write(*,'(A)') 'Uijkl before trafo:'
            call print_matrix(Uijkl,norbPerAtom(atom)**2)
         endif
   
         ! create the tensor transformation matrix
         do s1=1,nspin
            do a=0,norbPerAtom(atom)-1
               do b=0,norbPerAtom(atom)-1
                  do c=0,norbPerAtom(atom)-1
                     do d=0,norbPerAtom(atom)-1
               
                        m1 = a*norbPerAtom(atom)+b
                        m2 = c*norbPerAtom(atom)+d
   
                        Utensortrafo(m1+1,m2+1,s1)  = UtrafoA(a+1,c+1,s1) * UtrafoA(b+1,d+1,s1) 
                        UtensortrafoT(m1+1,m2+1,s1) = UtrafoTA(a+1,c+1,s1) * UtrafoTA(b+1,d+1,s1)  
                        
                     enddo ! l
                  enddo ! k
               enddo ! j
            enddo ! i
         enddo ! s1
        
         ! Now transform the Uijkl tensor    
         do s1=1,nspin
            do s2=1,nspin
    
               Uijkl_trafo(:,:,s1,s2) = matmul( UtensortrafoT(:,:,s1), matmul( Uijkl , Utensortrafo(:,:,s2) ) )
   
               if (diagBath==1) then
                  write(*,'(A,I2,A,I2)') 'Uijkl after trafo for s1=',s1,', s2=',s2
                  call print_matrix(Uijkl_trafo(:,:,s1,s2),norbPerAtom(atom)**2)
               endif
   
            enddo
         enddo
   
         ! =======================================================
         ! Transform mu_loc matrix ===============================
         ! =======================================================
         do s1=1,nspin
            call get_dmftpart_atom_cr(tmpmatrix, mu_loc(:,:,s1) ,atom)
            mu_loc_trafo(:,:,s1) = matmul( UtrafoTA(:,:,s1), matmul(tmpmatrix, UtrafoA(:,:,s1) ) )
   
            ! Make it hermitian to account for rounding errors
            mu_loc_trafo(:,:,s1) = 0.5_kr*( mu_loc_trafo(:,:,s1) + TRANSPOSE(CONJG(mu_loc_trafo(:,:,s1))) )
         enddo
       
         
   
         ! =======================================================
         ! Transform hybridization and gbath matrix ========================
         ! =======================================================
         do s1=1,nspin
            do n=1,nomega
               call get_dmftpart_atom(tmpmatrix, hybrid(:,:,s1,n) ,atom)
               hybrid_trafo(:,:,s1,n) = matmul( UtrafoTA(:,:,s1), matmul(tmpmatrix, UtrafoA(:,:,s1) ) )
            enddo
         enddo
   
         ! =======================================================
         ! Neglect remaining offdiagonal elements? ===============
         ! =======================================================
   !      if (neglectOffdiagBath==1) then
   !         write(*,'(A)') 'We neglect any (remaining) offdiagonal components in Hyb(iw) and mu_loc for the solver!'
   !         
   !         do s1=1,nspin
   !            do m1=1,norbPerAtom(atom)
   !               do m2=1,norbPerAtom(atom)
   !
   !                  if (m1/=m2) then
   !                     mu_loc_trafo(m1,m2,s1) = (0.0_kr,0.0_kr)
   !                     hybrid_trafo(m1,m2,s1,:) = (0.0_kr,0.0_kr)
   !                  endif
   !
   !               enddo ! m2
   !            enddo ! m1
   !         enddo ! s1
   !      endif
       
         ! =======================================================
         ! Write everything into files ===========================
         ! =======================================================
         ! FROM HERE ON:
         ! For the solver output we have a different ordering !!
         ! spin runs first, then orbitals !
         ! But we only write it out that way
         ! =============================================
   
         ! ==========================================================
         ! First write the hybridization function on imaginary time
         write(*,'(A)') ' Calculate imaginary time hybrid function...'
         if (atom>9)  write(path, "(A,I2,A)") "imp",atom,"/delta0"
         if (atom<10) write(path, "(A,I1,A)") "imp",atom,"/delta0"
         call write_imag_time_hybrid(trim(path),hybrid_trafo)
         ! ==========================================================
         ! And the bath Green's function on imaginary time, submit delta because we
         ! modify the muloc by subtracting U/2 for the CT-INT solver
         if (atom>9)  write(path, "(A,I2,A)") "imp",atom,"/gbath_tau"
         if (atom<10) write(path, "(A,I1,A)") "imp",atom,"/gbath_tau"
         write(*,'(A)') " Calculate imaginary time bath Green's function..."
         call write_imag_time_bath(trim(path),hybrid_trafo,mu_loc_trafo)
   
         ! ==========================================================
         ! Write Uijkl matrix
         if (atom>9)  write(path, "(A,I2,A)") "imp",atom,"/Uijkl.dat"
         if (atom<10) write(path, "(A,I1,A)") "imp",atom,"/Uijkl.dat"
         call write_Uijkl_matrix_solver(Uijkl_trafo,trim(path))
         ! ==========================================================


         ! ==========================================================
         ! Write density-density interaction matrix
         if (atom>9)  write(path, "(A,I2,A)") "imp",atom,"/umatrix_hybseg.dat"
         if (atom<10) write(path, "(A,I1,A)") "imp",atom,"/umatrix_hybseg.dat"
         call write_umatrix_segment(Uijkl_trafo,trim(path))
         ! ==========================================================

         ! ==========================================================
         ! Write mu_loc
         if (atom>9)  write(path, "(A,I2,A)") "imp",atom,"/mu_loc"
         if (atom<10) write(path, "(A,I1,A)") "imp",atom,"/mu_loc"
         call write_muloc_solver(mu_loc_trafo,trim(path))
         ! ==========================================================

   
         ! ==========================================================
!         ! Write mu_loc matrix and density-density interaction matrix
!         open(unit=mu_out,file="mu_vector.dat",status="unknown")
!         open(unit=mu_matrix_out,file="mu_matrix.dat",status="unknown")
!         open(unit=mu_matrix_solver_out,file="mu_matrix_solver.dat",status="unknown")
!         open(unit=umatrix_out,file="umatrix.dat",status="unknown")
!   
!         do m1=1,norbPerAtom(atom)
!            do s1=1,nspin
!               write(mu_out,'((F11.6),(4X))',advance='no') real(mu_loc_trafo(m1,m1,s1))
!               do m2=1,norbPerAtom(atom)
!                  do s2=1,nspin
!                     a = dmftOrbsIndex(m1)-1
!                     b = dmftOrbsIndex(m2)-1
!   
!                     !coeffs = get_highfreq_coeff_bosonic( uloc(a*norb+a+1,b*norb+b+1,1,:) )
!   
!                     !Uval = real(uloc(a*norb+a+1,b*norb+b+1,1,size(uloc(1,1,1,:))))
!                     !Uval = coeffs(1)
!                     Uval = real( Uijkl_trafo((m1-1)*norbPerAtom(atom)+m2,(m1-1)*norbPerAtom(atom)+m2,1,1) )   ! U, U'
!                     if (m1/=m2 .and. s1==s2) then
!                        Uval = Uval - real( Uijkl_trafo((m1-1)*norbPerAtom(atom)+m2,(m2-1)*norbPerAtom(atom)+m1,1,1) ) ! U'-J
!                     endif
!   
!                     write(umatrix_out,'((F10.5),(4X))',advance='no') Uval
!   
!                     if (s1==s2) then
!                        write(mu_matrix_solver_out,'(I3,4X,I3,4X,(F10.5),(4X),(F10.5))')          &
!                              & ((m1-1)*2+s1-1),((m2-1)*2+s2-1),  -real( mu_loc_trafo(m1,m2,s1)), &
!                              &                                   -aimag(mu_loc_trafo(m1,m2,s1))
!                     else  
!                        write(mu_matrix_solver_out,'(I3,4X,I3,4X,(F10.5),(4X),(F10.5))') &
!                                      & ((m1-1)*2+s1-1),((m2-1)*2+s2-1), 0.0, 0.0
!                     endif
!   
!                  enddo ! s2
!               enddo ! m2
!               write(umatrix_out,'(A1)') ' '
!            enddo ! s1
!         enddo ! m1
!   
!         ! general matrix format for plotting
!         do s1=1,nspin
!            do m1=1,norbPerAtom(atom)
!               do m2=1,norbPerAtom(atom)
!                write(mu_matrix_out,'(2(F11.6,4X))',advance='no') real(mu_loc_trafo( m1,m2,s1) ),&
!                                     &                            aimag(mu_loc_trafo(m1,m2,s1) )
!               enddo ! m2
!               write(mu_matrix_out,'(A1)') ' '
!            enddo ! m1
!            write(mu_matrix_out,'(A1)') ' '
!         enddo ! s1
!   
!         close(umatrix_out)
!         close(mu_out)
!         close(mu_matrix_out)
!         close(mu_matrix_solver_out)
!   
!   
!         ! ==========================================================
!         ! Now write K(tau) if needed
!         ! ==========================================================
!         if ( useUw/=0 ) then
!             ! now create the K(tau),K'(tau) functions, maybe put this somewhere else?
!             ! We assume Uavg = F0 and only use this
!             call get_Uavg(Uavg,uloc)
!             hfcoeff = get_highfreq_coeff_bosonic( Uavg ) ! get constant term
!   
!             write(*,'(A,F9.5)') 'Using a fitted F0 when generating K(tau) of : ' &
!                       & ,hfcoeff(1)
!   
!             open(unit=k_out,file="Kt.dat",status="unknown")
!             do it=0,ntau
!                K = 0.0_kr
!                Kp = 0.0_kr
!                tau = it*beta*1.0_kr/ntau
!                do n=1,nnu-1
!                   !K = K   - ( Uavg(n+1) - Ubare )*( cos(vn(n)*tau)-1 )*2.0_kr/(beta*vn(n)**2)
!                   !Kp = Kp + ( Uavg(n+1) - Ubare )*  sin(vn(n)*tau)    *2.0_kr/(beta*vn(n)   )
!                   K = K   - ( real(Uavg(n+1)) - hfcoeff(1) )*( cos(vn(n)*tau)-1 )*2.0_kr/(beta*vn(n)**2)
!                   Kp = Kp + ( real(Uavg(n+1)) - hfcoeff(1) )*  sin(vn(n)*tau)    *2.0_kr/(beta*vn(n)   )
!                enddo
!                ! n=0 case separately
!                !K = K   + ( Ubare - Uavg(1) )*(beta-tau)*tau/(2*beta)
!                !Kp = Kp + ( Ubare - Uavg(1) )*(beta-2*tau)/(2*beta)
!                K = K   + ( hfcoeff(1) - real(Uavg(1)) )*(beta-tau)*tau/(2*beta)
!                Kp = Kp + ( hfcoeff(1) - real(Uavg(1)) )*(beta-2*tau)/(2*beta)
!   
!                if (K<0.0_kr) K=0.0_kr
!   
!                write(k_out,'((ES14.7,3X,ES14.7,3X,ES14.7))') tau,K,Kp
!             enddo
!             close(k_out)
!         endif
         
   
         deallocate( mu_loc_trafo ) 
         deallocate( hybrid_trafo ) 
         deallocate( Uijkl ) 
         deallocate( Uijkl_trafo ) 
         deallocate( Utensortrafo ) 
         deallocate( UtensortrafoT ) 
         deallocate( UtrafoA ) 
         deallocate( UtrafoTA ) 
         deallocate( tmpmatrix ) 
      enddo !  atoms
      write(*,'(A)') ' '
   
   end subroutine write_solver_input
   
!  ============================================================
!  == Calls the solver: Adjust here to any solver, modify MPI threads!
!  ============================================================
   subroutine run_solver()

         call system("./runsolver.sh")

   end subroutine run_solver

!  ============================================================
!  == A simple Hartree-Fock solver for k-dependent interactions
!  == For now do a simple 2D lattice with NN interaction
!  ============================================================

   function V_q(qx, qy, qz)
      real(kr)            :: V_q
      real(kr),intent(in) :: qx,qy,qz

      V_q = 2* 1.0_kr * ( cos( qx ) + cos( qy )  )!*cos( qz/4 ) )

   end function V_q

!  =============================================================
! assume paramgnetic hamiltonian

!   subroutine hf_solver(s_gw,h_dft,umatrix,s_gw_dc,simp,dft_hartree,dft_exchange)
!       use matsfunc_ops, only : wn, filling
!       use mathfunc
!       complex(kr),intent(inout) :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
!       complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
!       real(kr),intent(in)       :: umatrix(:,:,:,:)    ! norbPerAtom,nspin,norbPerAtom,nspin
!       complex(kr),intent(inout)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
!       complex(kr),intent(inout)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
!       real(kr),intent(inout)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
!       real(kr),intent(inout)       :: dft_exchange(:,:,:) ! norb,norb,nspin
!
!       complex(8),allocatable :: gloc_tmp(:,:,:) ! norb,norb,nomega
!       real(kr),allocatable   :: unity(:,:)    ! norb,norb
!       real(kr),allocatable   :: nk(:,:,:)    ! norbPerAtom,spin,nkpts
!       integer(ki)            :: w,k,q,s,i,m1,m2,s1,s2, kx,ky,kz,qx,qy,qz, iq, ik, a
!       real(kr)               :: qxr,qyr,qzr
!
!    !   Lapack stuff
!       complex(8), dimension(norb) :: work  ! work array for LAPACK
!       integer, dimension(norb)    :: ipiv  ! pivot indices
!       integer :: info_lpck
!       external ZGETRF
!       external ZGETRI
!
!       if ( useSigK ==1 ) then
!          stop "ERROR: useSigK != 1 but HF solver is called! gw_sigma is not properly allocated !!!"
!       endif
!
!       s_gw_dc = (0.0_kr,0.0_kr)
!       simp = (0.0_kr,0.0_kr)
!       dft_hartree = (0.0_kr)
!       dft_exchange = (0.0_kr)
!
!       allocate( gloc_tmp(norb,norb,nomega) )
!       allocate( unity(norb,norb) )
!       allocate( nk(norbPerAtom,nspin,nkpts) )
!
!    !  ====================================================
!    !  == Create identity matrix for use with iw_n + mu ==
!       unity = (0.0_kr)
!       do i=1,norb
!          unity(i,i) = 1.0_kr
!       enddo
!    !  ====================================================
!
!!  ====================================================
!!  == Calculate the filling for all k ==
!   write(*,'(A)') 'Calculate the filling for all k...'
!
!   nk = (0.0_kr)
!
!   do k=1,nkpts
!      do s=1,nspin
!         do w=1,nomega
!            ! First create the inverse Green's function
!            gloc_tmp(:,:,w) = (ci*wn(w-1)+mu)*unity - h_dft(:,:,s,k) - s_gw(:,:,s,k,w)
!
!            ! Then invert with LAPACK
!            call ZGETRF(norb, norb, gloc_tmp(:,:,w), norb, ipiv, info_lpck)
!            if (info_lpck /= 0) then
!               write(*,'(A,I3)') 'ERROR: Greensfunction matrix is numerically singular! Return value', &
!                         & info_lpck
!               stop
!            end if
!            call ZGETRI(norb, gloc_tmp(:,:,w), norb, ipiv, work, norb, info_lpck)
!            if (info_lpck /= 0) then
!               stop 'Matrix inversion failed!'
!            end if
!        enddo
!
!        ! Then calculate filling
!        !do m1=1,norbPerAtom
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        do m1=1,6
!           nk(m1,s,k) = filling( gloc_tmp(10+m1,10+m1,:) )
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           !nk(m1,s,k) = filling( gloc_tmp(dmftOrbsIndex(m1),dmftOrbsIndex(m1),:) )
!           !nk(m1,s,k) = 0.95_kr
!        enddo
!
!      enddo ! s
!   enddo ! k
!
!       write(*,*) 'Filling at Gamma: UP'
!       do m1=1,norbPerAtom
!          write(*,'((ES10.3),(2X))',advance='no') nk(m1,1,1)
!       enddo
!       write(*,*) ' '
!       write(*,*) 'Filling at Gamma: DN'
!       do m1=1,norbPerAtom
!          write(*,'((ES10.3),(2X))',advance='no') nk(m1,2,1)
!       enddo
!       write(*,*) ' '
!
!       write(*,*) 'Filling at M: UP'
!       do m1=1,norbPerAtom
!          write(*,'((ES10.3),(2X))',advance='no') nk(m1,1, nky*nkz*nkx/2 + nkz*nky/2 +1 )
!       enddo
!       write(*,*) ' '
!       write(*,*) 'Filling at M: DN'
!       do m1=1,norbPerAtom
!          write(*,'((ES10.3),(2X))',advance='no') nk(m1,2, nky*nkz*nkx/2 + nkz*nky/2 +1 )
!       enddo
!       write(*,*) ' '
!
!!  ====================================================
!
!   deallocate( unity )
!   deallocate( gloc_tmp )
!
!       ! Then update the Selfenergy
!       write(*,*) 'Calculate Hartree-Fock Selfenergy...'
!       s_gw = (0.0_kr,0.0_kr)
!
!       do kx=0,nkx-1
!       do ky=0,nky-1
!       do kz=0,nkz-1
!          ik = nkz*nky*kx + nkz*ky + kz +1
!
!         call progress(ik,nkpts)
!
!          do m1=1,norbPerAtom
!             do s1=1,nspin
!                do qx=0,nkx-1
!                do qy=0,nky-1
!                do qz=0,nkz-1
!                   iq = nkz*nky*qx + nkz*qy + qz +1
!
!                   do m2=1,norbPerAtom
!              !        do s2=1,nspin
!              !           ! local term
!              !           s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1) = &
!              !      &    s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1)  + 0*umatrix(m1,s1,m2,s2) * nk(m2,s2,iq)
!              !        enddo
!
!                      ! different spin
!                      s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1) = &
!                 &    s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1) + V_q(0.0_kr,0.0_kr,0.0_kr) * nk(m2,3-s1,iq)
!
!                      ! same spin
!                      qxr = (kx-qx)*2.0*pi/nkx
!                      qyr = (ky-qy)*2.0*pi/nky
!                      qzr = (kz-qz)*2.0*pi/nkz
!
!!                      if (m1==m2) then
!!                         s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1) = &
!!                    &    s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1)   &
!!                    & + ( V_q(0.0_kr,0.0_kr,0.0_kr) - V_q(qxr,qyr,qzr) )* nk(m2,s1,iq)
!!                      else
!                          s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1) = &
!                     &    s_gw(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s1,ik,1)   &
!                     & +  V_q(0.0_kr,0.0_kr,0.0_kr) * nk(m2,s1,iq)
!!                      endif
!                    
!                   enddo
!                enddo
!                enddo
!                enddo
!            enddo
!          enddo
!       enddo
!       enddo
!       enddo
!       write(*,*) 'Done HF'
!
!       deallocate( nk )
!
!       ! now copy to the other equivalent atoms
!       do w=1,nomega
!          do ik=1,nkpts
!             do a=1,noEquivAtoms-1
!                do m1=1,norbPerAtom
!                   s_gw( dmftOrbsIndex(a*norbPerAtom+m1), dmftOrbsIndex(a*norbPerAtom+m1) ,:,ik,w)   &
!                           &  = s_gw( dmftOrbsIndex(m1),dmftOrbsIndex(m1) ,:,ik,w)
!                enddo
!             enddo
!          enddo
!       enddo
! 
!       write(*,*) 'Copy done'
!
!       s_gw(:,:,:,:,1) = s_gw(:,:,:,:,1)/nkpts
!       do w=2,nomega
!          do ik=1,nkpts
!             s_gw(:,:,:,ik,w) = s_gw(:,:,:,ik,1)
!          enddo
!       enddo
!
!
!       write(*,*) 'HF Selfenergy at Gamma:'  ! norb,norb,nspin,nkpts,nomega
!       do m1=1,norb
!          write(*,'((ES10.3),(2X))',advance='no') real(s_gw(m1,m1,1,1,1))
!       enddo
!       write(*,*) ' '
!       do m1=1,norb
!          write(*,'((ES10.3),(2X))',advance='no') real(s_gw(m1,m1,2,1,1))
!       enddo
!       write(*,*) ' '
!
!       write(*,*) 'HF Selfenergy at M:'
!       do m1=1,norb
!          write(*,'((ES10.3),(2X))',advance='no') real(s_gw(m1,m1,1,nky*nkz*nkx/2 + nkz*nky/2 +1,1))
!       enddo
!       write(*,*) ' '
!       do m1=1,norb
!          write(*,'((ES10.3),(2X))',advance='no') real(s_gw(m1,m1,2,nky*nkz*nkx/2 + nkz*nky/2 +1,1))
!       enddo
!       write(*,*) ' '
!
!!       do w=1,nomega
!!          do ik=1,nkpts
!!             s_gw(:,:,1,ik,w) = ( s_gw(:,:,1,ik,w) + s_gw(:,:,2,ik,w) )/2.0
!!             s_gw(:,:,2,ik,w) = s_gw(:,:,1,ik,w)
!!          enddo
!!       enddo
!!       stop 'AFTER HF SOLVER'
!
!   end subroutine hf_solver


!  ============================================================
!  == Solve using DCA
!  ============================================================
   subroutine solver_DCA(gimp,simp_new,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use io_ops
      use matsfunc_ops
      use fourier_interpolation
      complex(kr),intent(inout) :: gimp(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(inout) :: simp_new(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(inout) :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(inout) :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

      integer(ki)             :: x,y,i,j,k, kdim, ip,ikx,iky, s
      real(kr)                :: kpos(npx*npy,2), kx, ky, pos(npx*npy,2)
      complex(kr),allocatable :: g0K(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),allocatable :: g0KUshift(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),allocatable :: g0R(:,:,:,:)          ! npx*npy,npx*npy,nspin,nomega
      complex(kr),allocatable :: g0RUshift(:,:,:,:)          ! npx*npy,npx*npy,nspin,nomega
      complex(kr),allocatable :: gimpR(:,:,:,:)          ! npx*npy,npx*npy,nspin,nomega
      complex(kr),allocatable :: gimpsym(:,:,:)          ! npx*npy,npx*npy,nomega
      complex(kr),allocatable :: gimpK(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),allocatable :: simpK(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),allocatable :: s_gw_coarse(:,:,:,:,:)     ! norb,norb,nspin,npx*npy,nomega

      complex(kr),allocatable :: gbathtmp(:,:), gimptmp(:,:)          ! npx*npy

      allocate(g0K(npx*npy,nspin,nomega))
      allocate(g0KUshift(npx*npy,nspin,nomega))
      allocate(g0R(npx*npy,npx*npy,nspin,nomega))
      allocate(g0RUshift(npx*npy,npx*npy,nspin,nomega))
      allocate(gimpR(npx*npy,npx*npy,nspin,nomega))
      allocate(gimpK(npx*npy,nspin,nomega))
      allocate(simpK(npx*npy,nspin,nomega))
      allocate(gbathtmp(npx*npy,npx*npy))
      allocate(gimptmp(npx*npy,npx*npy))
      allocate( s_gw_coarse(1,1,nspin,npx*npy,nomega ) )

      k=1
      do i=0,npx-1
         do j=0,npy-1
            pos(k,:) = (/ i, j /)
            kpos(k,:) = (/ i*2*pi/npx, j*2*pi/npy /)
            !if (npx .ne. npy) then
            !   pos(k,:) = (/ i+j, j /)/2
            !   kpos(k,:) = (/ i*2*pi/npx + j*2*pi/npy , j*2*pi/npy /)
            !endif
            k=k+1
         enddo
      enddo

      call get_g0K_DCA(g0K,g0KUshift,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      call write_loc_matsfuncvec("g0K_DCA_Ushift.dat",g0KUshift)
      call write_loc_matsfuncvec("g0K_DCA.dat",g0K)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Fourier trafo of g0K
      g0R = (0.0_kr,0.0_kr)
      g0RUshift = (0.0_kr,0.0_kr)
      do i=1,npx*npy
         do j=1,npx*npy
            x = pos(i,1) - pos(j,1)
            y = pos(i,2) - pos(j,2)
   
            do k=1,npx*npy
               kx = kpos(k,1)
               ky = kpos(k,2)
               g0R(i,j,:,:)       = g0R(i,j,:,:)       + g0K(k,:,:)*exp(-ci*( x*kx+y*ky ) )  
               g0RUshift(i,j,:,:) = g0RUshift(i,j,:,:) + g0KUshift(k,:,:)*exp(-ci*( x*kx+y*ky ) )  
            enddo
         enddo ! iy
      enddo ! ix

      g0R = g0R/(npx*npy)
      g0RUshift = g0RUshift/(npx*npy)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call write_loc_matsfunc_matrix("g0R_DCA.dat",g0R)
      call write_loc_matsfunc_matrix("g0RUshift_DCA.dat",g0RUshift)
      call transform_to_imag_time("g0R_tau.dat",g0RUshift)

      call run_solver()
      call read_loc_matsfunc_matrix("Gwl.dat",gimpR)

      ! symmetrize spin components
      if (spinorder==1) then
         write(*,'(A)') ' Symmetrize spin components of the DCA cluster GF (PM)...'
         allocate ( gimpsym(npx*npy,npx*npy,nomega) )
         gimpsym = (0.0_kr,0.0_kr)
         do s=1,nspin
            gimpsym = gimpsym + gimpR(:,:,s,:)
         enddo
         do s=1,nspin
            gimpR(:,:,s,:) = gimpsym/nspin
         enddo
         deallocate(gimpsym)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Fourier trafo from gimpR to gimpK
      gimpK = (0.0_kr,0.0_kr)
      do k=1,npx*npy
         kx = kpos(k,1)
         ky = kpos(k,2)
         do i=1,npx*npy
            do j=1,npx*npy
               x = pos(i,1) - pos(j,1)
               y = pos(i,2) - pos(j,2)
      
               gimpK(k,:,:) = gimpK(k,:,:) + gimpR(i,j,:,:)*exp(+ci*( x*kx+y*ky ) )   
            enddo
         enddo ! iy
      enddo ! ix
 
      gimpK = gimpK/(npx*npy)
      call write_loc_matsfuncvec("g_impK_DCA.dat",gimpK)
 
      ! Get selfenergy
      simp_new = (0.0_kr,0.0_kr)
      gimp = (0.0_kr,0.0_kr)
      do k=1,npx*npy
         simpK(k,:,:) = 1.0/g0K(k,:,:) - 1.0/gimpK(k,:,:)
 
         simp_new(1,1,:,:) = simp_new(1,1,:,:) + simpK(k,:,:)/(npx*npy)
         gimp(1,1,:,:) = gimp(1,1,:,:) + gimpK(k,:,:)/(npx*npy)

         s_gw_coarse(1,1,:,k,:)  = simpK(k,:,:)
      enddo
      call write_loc_matsfuncvec("s_impK_DCA.dat",simpK)

      ! copy selfenergy to s_gw, this is still patchy
      do ikx=0,nkx-1
         do iky=0,nky-1
            k = ikx*nky + iky +1 
            !ip = patchindex(ikx,iky)
            ip = kPatchIndex(k)
            s_gw(1,1,:,k,:)     = simpK(ip,:,:)
            ! norb,norb,nspin,nkpts,nomega

         enddo
      enddo      
      !call write_kdep_static_f("sgw_DCA_patch.dat",s_gw(:,:,:,:,1),0)
      call write_kdep_w0_f("sgw_DCA_patch.dat",s_gw)
      
      ! or interpolate the DCA selfenergy ??
      if (interpDCA==1) then
         if (npx .ne. npy) stop 'ERROR: Interpolation of DCA Selfenergy only possible for npx==npy !'

         do i=1,nomega
            call interpolate_kmatrix_cubic(s_gw(:,:,:,:,i), s_gw_coarse(:,:,:,:,i),nkx,nky,nkz,npx,npy,1)
         enddo
      endif
      !call write_kdep_static_f("sgw_DCA_interpolated.dat",s_gw(:,:,:,:,1),0)
      call write_kdep_w0_f("sgw_DCA_interpolated.dat",s_gw)

      ! set doublecounting to local part of DCA selfenrgy
      s_gw_dc = simp_new
      
!      do k=1,nkpts
!         s_gw(1,1,:,k,:)  = s_gw(1,1,:,k,:) - simp_new(1,1,:,:) 

         !!!!!!!!!!!
!         s_gw(1,1,:,k,1)  = s_gw(1,1,:,k,1) - simp_new(1,1,:,1) 
!         do i=2,nomega
!            s_gw(1,1,:,k,i)  = s_gw(1,1,:,k,1)
!         enddo
!      enddo      

      deallocate(g0K)
      deallocate(g0KUshift)
      deallocate(g0R)
      deallocate(g0RUshift)
      deallocate(gimpK)
      deallocate(gimpR)
      deallocate(simpK)
      deallocate(gbathtmp)
      deallocate(gimptmp)
      deallocate( s_gw_coarse )

   end subroutine solver_DCA
   
end module solver
