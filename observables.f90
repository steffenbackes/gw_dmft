!  ============================================================
!  == Everything related to observables, Pade, Fermisurface, eff.masses etc
!  ============================================================

module observables
   use constants
   implicit none

 contains
 !  ============================================================
 !  == Write Pade anacont of a local Matsubara function into a file
 !  ============================================================
    subroutine write_pade_loc(gf,filename,fermionic)
       use params
       use anacont
       complex(kr),intent(in)       :: gf(:,:,:,:)    ! norb,norb,nspin,nomega
       character(len=*), intent(in) :: filename
       logical,intent(in)           :: fermionic

       integer                 :: iounit=10,n,m,s, realnw = 1000
       real(kr)                :: dx,              a=-20.0_kr, b=10.0_kr
       complex(kr),allocatable :: cont_gf(:,:,:)

       allocate( cont_gf(norb,nspin,realnw) )
       dx = (b-a)/realnw

       ! get continuation of all orbitals, spins, only diagonals
       do m=1,norb
          do s=1,nspin
             call init_pade(gf(m,m,s,:),nomega,fermionic)
             do n=1,realnw
                cont_gf(m,s,n) = get_pade(a + dx*n)
             enddo
          enddo
       enddo

       ! Then write all stuff into a file
       open(unit=iounit,file=filename,status="unknown")
       do n=1,realnw
          write(iounit,'(ES12.5,3X)',advance='no') a + dx*n
          do s=1,nspin
             do m=1,norb
                write(iounit,'(ES12.5,3X,ES12.5,3X)',advance='no') cont_gf(m,s,n)
             enddo
          enddo
          write(iounit,'(A)') ' '
       enddo
       close(iounit)
       deallocate( cont_gf )
    end subroutine write_pade_loc


!  ============================================================
!  == Write the k-resolved pade-continued greensfunction into a file
!  == We only write the components for the first equivalent atom to save space
!  ============================================================
   subroutine write_pade_bandstructure(filename,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use params
      use anacont
      use matsfunc_ops
!      use io_ops, only : write_loc_matsfunc
      character(len=*),intent(in) :: filename
      complex(kr),intent(in)      :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)      :: s_gw_dc(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(kr),intent(in)      :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)         :: dft_hartree(:,:,:) ! norb,norb,nspin
      real(kr),intent(in)         :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)          :: onlyDFT              ! create only the DFT Green's function

      integer                 :: iounit=10,iounitg=11,   realnw = 500
      real(kr)                :: dx,              a=-6.0_kr, b=6.0_kr
      complex(kr),allocatable :: gf_tmp(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(8),allocatable  :: inv_gf(:,:)     ! norb, for inverting
      complex(kr),allocatable :: cont_gf(:,:,:)   ! norb,nspin,realnw
      real(kr),allocatable    :: unity(:,:)    ! norb,norb
      integer(ki)             :: w,k,ik,s,i,n,m

  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gf_tmp(norb,norb,nspin,nomega) )
      allocate( inv_gf(norb,norb) )      
      allocate( unity(norb,norb) )

      ! if there are uncorrelated orbitals, increase the lower bound...
      if (norb /= norb_dmft) a = -10.0_kr

!  ====================================================
!  == Create identity matrix for use with iw_n + mu ==
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
      allocate( cont_gf(norb,nspin,realnw) )
      dx = (b-a)/realnw
!  ====================================================

!  ====================================================
      write(*,'(A)') 'Calculate the bandstructure with Pade...'
!     Here one should parallelize over something

      open(unit=iounit,file=filename,status="unknown")
      open(unit=iounitg,file=trim(filename) // "_imag",status="unknown")

      do ik=0,nkx/2 + nky/2 + nky/2 + nkz/2
          ! This is really bad stuff, change this
          k = ik*nkz*nky +1

          ! X to M
          if (ik .ge. nkx/2) k = (nkx/2)*nkz*nky  + (ik-nkx/2)*nkz +1

          ! M to Gamma
          if (ik .ge. nkx/2 + nky/2) k = (nkx/2)*nkz*nky + (nky/2)*nkz - (ik-nkx/2-nky/2)*( nkz*nky + nkz ) +1

          ! Gamma to Z
          !if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = ik - (nkx/2 + nky/2 + nky/2) +1

          ! Gamma to R
          if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = (ik-(nkx/2+nky/2+nky/2))*(  nkz*nky + nkz  + 1 )    +1

         do w=1,nomega
            do s=1,nspin
                ! First create the inverse DFT Green's function
!                inv_gf = (ci*wn(w-1)+mu + ci*0.01 + hfield(s))*unity - h_dft(:,:,s,k)

                ! Then add all GW, DC and DMFT contributions
!                if (.not. onlyDFT) then
!                   inv_gf = inv_gf                   &
!       &                    + dft_hartree(:,:,s)     & ! Subtract the DFT Hartree term (no Hartree and Exchange left)
!       &                    - simp(:,:,s,w)            ! Add the impurity DMFT Selfenergy including Hartree-Fock
!
!                  if ( useGW == 1 ) then
!                     inv_gf = inv_gf                 &
!       &                    + v_xc(:,:,s,k)          & ! Substract V_XC from DFT
!       &                    - s_gw(:,:,s,k,w)        & ! Add the GW ex-cor. Selfenergy (add Exchange back)
!       &                    + s_gw_dc(:,:,s,w)         ! remove the doublecounting between GW and DMFT (remove local exchange)
!                  else
!                      ! if we do LDA+DMFT without GW, we also want to replace the local exchange by DMFT
!                      inv_gf = inv_gf + dft_exchange(:,:,s)   ! Subtract the DFT Exchange term, its already in simp
!                   endif
!                endif
!
!               ! Then invert with LAPACK
!               call ZGETRF(norb, norb, inv_gf, norb, ipiv, info_lpck)
!               if (info_lpck /= 0) then
!                  write(*,*)' ERROR: Greensfunction matrix is numerically singular! Return value',&
!                            & info_lpck
!                  stop
!               end if
!               call ZGETRI(norb, inv_gf, norb, ipiv, work, norb, info_lpck)
!              if (info_lpck /= 0) then
!                 stop 'Matrix inversion failed!'
!              end if
!              
!              gf_tmp(:,:,s,w) = inv_gf
               gf_tmp(:,:,s,w) = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)

              ! Then we have a proper G(k,omega) orbital matrix for given spin
            enddo ! spin loop
            ! Now we have the full Greensfunction stored in gf_tmp(:,:,s,:) norb,norb,nspin
            ! for a given frequency and k-point
         enddo ! Matsubara frequency loop
         ! Now we have everything for a given k-point, so let's do the analytic continuation

! Write out the Greensfunction at Gamma and X !!!!!!!!!!!!!!!!
!   if (ik .eq. 0) then
!     call write_loc_matsfunc("g_gamma.dat",gf_tmp)
!     call write_pade_loc(gf_tmp,"g_gamma_pade.dat",.true.)
!   end if 
!   if (ik .eq. nkx/2) then
!     call write_loc_matsfunc("g_X.dat",gf_tmp)
!     call write_pade_loc(gf_tmp,"g_X_pade.dat",.true.)
!   end if

         do w=1,nomega
            write(iounitg,'(I5,3X,ES12.5,3X)',advance='no') ik, wn(w-1)
            do s=1,nspin
               do m=1,norb/noEquivAtoms
                  write(iounitg,'(ES12.5,3X,ES12.5,3X)',advance='no') real(gf_tmp(m,m,s,w)),aimag(gf_tmp(m,m,s,w))
               enddo
            enddo
            write(iounitg,'(A)') ' '
         enddo
         write(iounitg,'(A)') ' '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! get continuation of all orbitals, spins, only diagonals
         do m=1,norb/noEquivAtoms
            do s=1,nspin
               call init_pade( gf_tmp(m,m,s,:),nomega,.true.)
               do n=1,realnw
                  cont_gf(m,s,n) = get_pade(a + dx*n )
               enddo
            enddo
         enddo

        ! Then write all stuff into a file
!        open(unit=iounit,file=filename,status="unknown")
        do n=1,realnw
           write(iounit,'(I5,3X,ES12.5,3X)',advance='no') ik, a + dx*n
           do s=1,nspin
              do m=1,norb/noEquivAtoms
                 write(iounit,'(ES12.5,3X)',advance='no') -AIMAG(cont_gf(m,s,n))/pi
              enddo
           enddo
           write(iounit,'(A)') ' '
         enddo
         write(iounit,'(A)') ' '

       enddo ! k-loop
!stop 'After writing bandstructure'
       close(iounit)
       close(iounitg)
!  ====================================================

      deallocate( unity )
      deallocate( gf_tmp )
      deallocate( inv_gf )
      deallocate( cont_gf )
   end subroutine write_pade_bandstructure

!  ============================================================
!  == Write the k-resolved pade-continued greensfunction into a file
!  == We only write the components for the first equivalent atom to save space
!  ============================================================
   subroutine write_fermisurface(filename,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use params
      use anacont
      use matsfunc_ops
!      use io_ops, only : write_loc_matsfunc
      character(len=*),intent(in) :: filename
      complex(kr),intent(in)      :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)      :: s_gw_dc(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(kr),intent(in)      :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)         :: dft_hartree(:,:,:) ! norb,norb,nspin
      real(kr),intent(in)         :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)          :: onlyDFT

      integer                 :: iounit=10
      complex(kr),allocatable :: gf_tmp(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(8),allocatable  :: inv_gf(:,:)     ! norb, for inverting
      real(kr),allocatable    :: unity(:,:)    ! norb,norb
      integer(ki)             :: w,k,ik,s,i,n,m

      ! matrices for least square fit
      real(8) :: Xmat(6,5), Xmatinv(5,5), betaMat(5), Ymat(6)
      real(8), dimension(5) :: ls_work  ! work array for LAPACK
      integer, dimension(5) :: ls_ipiv  ! pivot indices

  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gf_tmp(norb,norb,nspin,nomega) )
      allocate( inv_gf(norb,norb) )      
      allocate( unity(norb,norb) )

!  ====================================================
!  == Create identity matrix for use with iw_n + mu ==
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

!  ====================================================
      write(*,*) 'Calculate the Fermisurface with Pade...'
      write(*,*) ' '
!     Here one should parallelize over something

      open(unit=iounit,file=filename,status="unknown")

      do ik=0,nkx*nky-1
         k = ( mod(ik,nkx) )*nkz*nky + int(ik/nkx)*nkz  +1

         ! For the Fermisurface we don't need more than 6 since we will only fit
         ! the first 6 frequencies
         do w=1,6
            do s=1,nspin
                ! First create the inverse DFT Green's function
!                inv_gf = (ci*wn(w-1)+mu + ci*0.003 + hfield(s) )*unity - h_dft(:,:,s,k)
!
!                ! Then add all GW, DC and DMFT contributions
!                if (.not. onlyDFT) then
!                   inv_gf = inv_gf                   & 
!       &                    + dft_hartree(:,:,s)     & ! Subtract the DFT Hartree term
!       &                    - simp(:,:,s,w)            ! Add the impurity DMFT Selfenergy including Hartree-Fock
!
!                  if ( useGW == 1 ) then
!                   inv_gf = inv_gf                   & 
!       &                    + v_xc(:,:,s,k)          & ! Substract V_XC from DFT
!       &                    - s_gw(:,:,s,k,w)        & ! Add the GW ex-cor. Selfenergy
!       &                    + s_gw_dc(:,:,s,w)         ! remove the doublecounting between GW and DMFT
!                   else
!                      ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                      inv_gf = inv_gf + dft_exchange(:,:,s)   ! Subtract the DFT Hartree term, its already in simp
!                   endif
!                endif ! (.not. onlyDFT)
!
!               ! Then invert with LAPACK
!               call ZGETRF(norb, norb, inv_gf, norb, ipiv, info_lpck)
!               if (info_lpck /= 0) then
!                  write(*,*)' ERROR: Greensfunction matrix is numerically singular! Return value',&
!                            & info_lpck
!                  stop
!               end if
!               call ZGETRI(norb, inv_gf, norb, ipiv, work, norb, info_lpck)
!              if (info_lpck /= 0) then
!                 stop 'Matrix inversion failed!'
!              end if
!              
!              gf_tmp(:,:,s,w) = inv_gf
              gf_tmp(:,:,s,w) = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)

              ! Then we have a proper G(k,omega) orbital matrix for given spin
            enddo ! spin loop
            ! Now we have the full Greensfunction stored in gf_tmp(:,:,s,:) norb,norb,nspin
            ! for a given frequency and k-point
         enddo ! Matsubara frequency loop
         ! Now we have everything for a given k-point, so let's do the analytic continuation

         ! get continuation of all orbitals, spins, only diagonals
         write(iounit,'(I4,3X,I4,3X)',advance='no') mod(ik,nkx), int(ik/nkx)

         do s=1,nspin
            do m=1,norb
               !write(*,*) 'm=',m,' s=',s
               ! Use Pade
               !call init_pade( gf_tmp(m,m,s,:),nomega,.true.)
               !write(iounit,'(ES10.3,3X)',advance='no') -AIMAG( get_pade( 0.0_kr ) )/pi
               
               ! Use linear fit
               !write(iounit,'(ES10.3,3X)',advance='no') -AIMAG( 1.5*gf_tmp(m,m,s,1)-0.5*gf_tmp(m,m,s,2) )/pi
               
               ! Use least square fit with polynomial of 4th order to the 6
               ! lowest frequencies
               do i=0,5 ! frequencies
                  Ymat(i+1) = aimag(gf_tmp(m,m,s,i+1))
                  do n=0,4 ! order of the monome
                     Xmat(i+1,n+1) = wn(i)**n
                  enddo
               enddo
               !write(*,*) 'Ymat= ',Ymat
               !write(*,*) 'Xmat= ',Xmat

               Xmatinv = matmul( transpose(Xmat),Xmat )
               !write(*,*) 'Xmatinv before= ',Xmatinv
               !write(*,*) 'Xmatinv before shape= ',shape(Xmatinv)
               call DGETRF(5, 5, Xmatinv, 5, ls_ipiv, info_lpck)
               call DGETRI(5, Xmatinv, 5, ls_ipiv, ls_work, 5, info_lpck)
               !write(*,*) 'Xmatinv= ',Xmatinv
               betaMat = matmul( Xmatinv, matmul(transpose(Xmat), Ymat) )
               !write(*,*) 'betaMat= ',betaMat
               write(iounit,'(ES10.3,3X)',advance='no') -betaMat(1)/pi
            enddo
         enddo
         write(iounit,'(A)') ' '
         if ( mod(ik+1,nkx)==0 ) write(iounit,'(A)') ' '

       enddo ! k-loop
       close(iounit)
!  ====================================================

      deallocate( unity )
      deallocate( gf_tmp )
      deallocate( inv_gf )
   end subroutine write_fermisurface


!  ============================================================
!  == Write the k-resolved pade-continued greensfunction into a file
!  == We only write the components for the first equivalent atom to save space
!  == INterpolate the Selfenergy instead of the Green's function
!  ============================================================
   subroutine write_fermisurface_sinterp(filename,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use params
      use anacont
      use matsfunc_ops
!      use io_ops, only : write_loc_matsfunc
      character(len=*),intent(in) :: filename
      complex(kr),intent(in)      :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)      :: s_gw_dc(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(kr),intent(in)      :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)         :: dft_hartree(:,:,:) ! norb,norb,nspin
      real(kr),intent(in)         :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)          :: onlyDFT

      integer                 :: iounit=10
      complex(8),allocatable  :: inv_gf(:,:), stmp(:,:)     ! norb, for inverting
      real(kr),allocatable    :: unity(:,:)    ! norb,norb
      integer(ki)             :: w,k,ik,s,i,n,m, m1,m2

      ! matrices for least square fit
      real(8) :: Xmat(6,5), Xmatinv(5,5), betaMat(5), YmatR(6), YmatI(6)
      real(8), dimension(5) :: ls_work  ! work array for LAPACK
      integer, dimension(5) :: ls_ipiv  ! pivot indices

  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( inv_gf(norb,norb) )      
      allocate( stmp(norb,norb) )      
      allocate( unity(norb,norb) )

!  ====================================================
!  == Create identity matrix for use with iw_n + mu ==
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

!  ====================================================
      write(*,*) 'Calculate the Fermisurface with Pade...'
      write(*,*) ' '
!     Here one should parallelize over something

      write(*,*) ' WARNING: Nonlocal Selfenergies in fermisurface_sinterp set to w=1 !...'
      w = 1

      open(unit=iounit,file=filename,status="unknown")

      do ik=0,nkx*nky-1
         k = ( mod(ik,nkx) )*nkz*nky + int(ik/nkx)*nkz  +1
         write(iounit,'(I4,3X,I4,3X)',advance='no') mod(ik,nkx), int(ik/nkx)

         do s=1,nspin

               ! Use least square fit with polynomial of 4th order to the 6
               ! lowest frequencies

               stmp = (0.0_kr,0.0_kr)

               do m1=1,norb
                  do m2=1,norb

                     do i=0,5 ! frequencies
                        YmatR(i+1) = real(simp(m1,m2,s,i+1))
                        YmatI(i+1) = aimag(simp(m1,m2,s,i+1))
                        do n=0,4 ! order of the monome
                           Xmat(i+1,n+1) = wn(i)**n
                        enddo
                     enddo
                     Xmatinv = matmul( transpose(Xmat),Xmat )
                     call DGETRF(5, 5, Xmatinv, 5, ls_ipiv, info_lpck)
                     call DGETRI(5, Xmatinv, 5, ls_ipiv, ls_work, 5, info_lpck)

                     ! Real part
                     betaMat = matmul( Xmatinv, matmul(transpose(Xmat), YmatR) )
                     stmp(m1,m2) = betaMat(1)
                     
                     ! Imag part
                     betaMat = matmul( Xmatinv, matmul(transpose(Xmat), YmatI) )
                     stmp(m1,m2) = stmp(m1,m2) + betaMat(1)*ci

                  enddo ! m2
                  ! sanity check: Diagonal Imaginary part should always be negative
                  if ( aimag(stmp(m1,m1)) > 0.0_kr ) then
                     stmp(m1,m1) = real(stmp(m1,m1)) - 0.003*ci
                  endif
               enddo ! m1

                ! First create the inverse DFT Green's function, we are at w=0
!                inv_gf = ( mu + ci*0.003 + hfield(s) )*unity - h_dft(:,:,s,k)
!
!                ! Then add all GW, DC and DMFT contributions
!                if (.not. onlyDFT) then
!                   inv_gf = inv_gf                   & 
!       &                    + dft_hartree(:,:,s)     & ! Subtract the DFT Hartree term
!       &                    - stmp                     ! Add the impurity DMFT Selfenergy including Hartree-Fock
!
!                  if ( useGW == 1 ) then
!                   inv_gf = inv_gf                   & 
!       &                    + v_xc(:,:,s,k)          & ! Substract V_XC from DFT
!       &                    - s_gw(:,:,s,k,w)        & ! Add the GW ex-cor. Selfenergy
!       &                    + s_gw_dc(:,:,s,w)         ! remove the doublecounting between GW and DMFT
!                   else
!                      ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                      inv_gf = inv_gf + dft_exchange(:,:,s)   ! Subtract the DFT Hartree term, its already in simp
!                   endif
!                endif ! (.not. onlyDFT)
!
!               ! Then invert with LAPACK
!               call ZGETRF(norb, norb, inv_gf, norb, ipiv, info_lpck)
!               if (info_lpck /= 0) then
!                  write(*,*)' ERROR: Greensfunction matrix is numerically singular! Return value',&
!                            & info_lpck
!                  stop
!               end if
!               call ZGETRI(norb, inv_gf, norb, ipiv, work, norb, info_lpck)
!               if (info_lpck /= 0) then
!                  stop 'Matrix inversion failed!'
!               end if
                inv_gf = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
              
               do m=1,norb
                  write(iounit,'(ES10.3,3X)',advance='no') -aimag(inv_gf(m,m))/pi
               enddo

         enddo ! spin loop

         write(iounit,'(A)') ' '
         if ( mod(ik+1,nkx)==0 ) write(iounit,'(A)') ' '

       enddo ! k-loop
       close(iounit)
!  ====================================================

      deallocate( unity )
      deallocate( inv_gf )
      deallocate( stmp )
   end subroutine write_fermisurface_sinterp

!  ============================================================
!  == Calculate and print the total energy calculated from Migdal Formula
!  ============================================================
   subroutine print_Etot(h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use params
      use matsfunc_ops
      complex(kr),intent(in)      :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)      :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)      :: s_gw_dc(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(kr),intent(in)      :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)         :: dft_hartree(:,:,:) ! norb,norb,nspin
      real(kr),intent(in)         :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)          :: onlyDFT

      complex(kr),allocatable     :: SigGw(:,:,:), g_tmp(:,:,:), trSigGw(:)
      integer(ki)     :: w,k,s,m
      complex(kr)     :: density(norb,norb),Hdens(norb,norb),H0(norb,norb)
      complex(kr)     :: unity(norb,norb), sigk(norb,norb), Etot, trHdens, trSigG
      real(kr)        :: coeffs(5)

      write(*,'(A)') " "
      write(*,'(A)') "Total Energy calculation:"

      allocate( SigGw(norb,norb,nomega) )
      allocate( g_tmp(norb,norb,nomega) )
      allocate( trSigGw(nomega) )

      unity = (0.0_kr,0.0_kr)
      do m=1,norb
         unity(m,m) = 1.0_kr
      enddo

      Etot = 0.0_kr
      do s=1,nspin
         SigGw = (0.0_kr,0.0_kr)
         density = (0.0_kr,0.0_kr)
         Hdens = (0.0_kr,0.0_kr)
         do k=1,nkpts   

            ! Create G and Sigma
            do w=1,nomega
               g_tmp(:,:,w) = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)

               sigk = simp(:,:,s,w)
               if ( useSigK/=0 )then
                  sigk = sigk + s_gw(:,:,s,k,w) - s_gw_dc(:,:,s,w) ! add nonlocal selfenergy and subtract Doublecounting
               endif

               SigGw(:,:,w) = SigGw(:,:,w) + matmul( sigk, g_tmp(:,:,w) )/nkpts

            enddo ! w

            ! Create Density Matrix
            call density_matrix(density,g_tmp)

            ! Create H0
            H0 = h_dft(:,:,s,k) - ( mu + hfield(s) )*unity 
            ! Then add all static contributions that consitute H0
            if (.not. onlyDFT) then
               H0 = H0  - dft_hartree(:,:,s) - dft_exchange(:,:,s) ! Subtract the DFT Hartree term (no Hartree and Exchange left)
               if ( useSigK/=0 )then
                  H0 = H0 - v_xc(:,:,s,k)       ! Subtract V_XC from DFT
               endif
           endif ! (.not. onlyDFT)

           Hdens = Hdens + matmul(H0,density)/nkpts

         enddo ! k
      
         ! Now we have Hdens and SigGw, now sum over orbitals so we only have one freq.sum at the end
         trHdens = 0.0_kr
         trSigGw = (0.0_kr,0.0_kr)
         do m=1,norb
            trHdens = trHdens + Hdens(m,m)
            trSigGw = trSigGw + SigGw(m,m,:)
         enddo

         ! now frequency sum of trSigGw
         coeffs = get_highfreq_coeff(trSigGw,666)
         !write(*,*) coeffs

         trSigG = 0.5_kr*coeffs(2) - coeffs(3)*beta/4 + coeffs(5)*(beta**3)/48
         do w=1,nomega
            trSigG = trSigG + ( real(trSigGw(w)) + coeffs(3)/wn(w-1)**2 - coeffs(5)/wn(w-1)**4  )*2.0_kr/beta
         enddo
         ! Done!

         write(*,'(A,I2,A,F9.3,A,F9.3)') "spin=",s,": trHdens=",real(trHdens),",  trSigG=",real(0.5*trSigG)
         Etot = Etot + trHdens + 0.5_kr*trSigG

      enddo ! s
      write(*,'(A,F12.5)') "Etot = ",real(Etot)
      write(*,'(A)') " "

      deallocate( SigGw )
      deallocate( g_tmp )
      deallocate( trSigGw)

   end subroutine print_Etot

end module observables
