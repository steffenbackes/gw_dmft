!  ============================================================
!  == Operations on Matsubara functions
!  ============================================================

module matsfunc_ops
   use constants
   use params
   use mathfunc
   implicit none

 contains

!  ============================================================
!  == Extracts the submatrix of DMFT orbitals on a given atom
!  ============================================================   
   subroutine get_dmftpart_atom(fout,fin,atom)
      complex(kr),intent(inout) :: fout(:,:)             ! norbPerAtom,norbPerAtom
      complex(kr),intent(in)    :: fin(:,:)             ! norb,norb
      integer(ki),intent(in)    :: atom

      integer(ki) :: a,m1,m2,i,j

      do m1=1,norbPerAtom(atom)
         do m2=1,norbPerAtom(atom)
            i = norb_dmft(atom,m1)
            j = norb_dmft(atom,m2)
            fout(m1,m2) = fin(i,j)
         enddo
      enddo
   end subroutine get_dmftpart_atom

!  ============================================================
!  == Extracts the submatrix of DMFT orbitals on a given atom for real functions
!  ============================================================   
   subroutine get_dmftpart_atom_cr(fout,fin,atom)
      complex(kr),intent(inout)   :: fout(:,:)             ! norbPerAtom,norbPerAtom
      real(kr),intent(in)       :: fin(:,:)             ! norb,norb
      integer(ki),intent(in)    :: atom

      integer(ki) :: a,m1,m2,i,j

      do m1=1,norbPerAtom(atom)
         do m2=1,norbPerAtom(atom)
            i = norb_dmft(atom,m1)
            j = norb_dmft(atom,m2)
            fout(m1,m2) = fin(i,j)
         enddo
      enddo
   end subroutine get_dmftpart_atom_cr

!  ============================================================
!  == Calculate the Doublecounting between the GW and DMFT Selfenergy
!  ============================================================
   subroutine get_sigma_gw_doublecounting(s_gw_dc,s_gw,gloc,wloc)
      complex(kr),intent(inout) :: s_gw_dc(:,:,:,:)   ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)    ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: gloc(:,:,:,:)      ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: wloc(:,:,:,:)      ! norb**2,norb**2,nspin,nnu

      integer                 :: w,k,s,m1,m2,a
      real(kr)                :: coeffs(5)
      complex(kr),allocatable :: dc_tmp(:,:,:,:)   ! norb,norb,nspin,nomega

!!!!!!!!!!!!!
             integer(ki) :: m
             integer(ki),parameter :: iounit = 10
             real(kr)              :: readindata(1+2*nspin*norb*norb) ! temparray for readin
             complex(kr)           :: tmp
!!!!!!!!!!!!!

      s_gw_dc = (0.0_kr, 0.0_kr)

      ! None of the implemented DC below have a Hartree part, like the GW Selfenergy
      ! Take care about this when implementing new types of DC's

      if (useSigK/=0) then

         allocate( dc_tmp(norb,norb,nspin,nomega) )
         dc_tmp  = (0.0_kr,0.0_kr)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! which type of DC?
         if (DCtype==1) then
            ! This is the [GW]_loc doublecounting
            ! First get the full local GW Selfenergy but take only the
            ! diagonal elements for the DMFT orbitals, rest is zero

            write(*,'(A)') ' Calculate the doublecounting as [GW]_loc on the DMFT atoms ...'
            write(*,'(A)') ' We assume the DC to be orbital-diagonal, since s_imp is ...'

            call get_locpart_mats(dc_tmp,s_gw)

            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  m2 = norb_dmft(a,m1)
                  s_gw_dc(m2,m2,:,:) = dc_tmp(m2,m2,:,:)
               enddo
            enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else if (DCtype==2) then
             write(*,'(A)') 'WARNING: DCtype = 2 has yet to be implemented !!!'
             write(*,'(A)') ' For now we just read the data from s_gw_dc.dat !!!'

      open(unit=iounit,file="s_gw_dc.dat",status="old")
      do w=1,nomega
         read(iounit,*) readindata
         do s=0,nspin-1
            do m1=0,norb-1
               do m2=0,norb-1
                 s_gw_dc(m1+1,m2+1,s+1,w)          &
                         &        =     readindata(1 + s*norb*norb*2 + m1*norb*2 + m2*2 +1)  &
                         &         + ci*readindata(1 + s*norb*norb*2 + m1*norb*2 + m2*2 +2)
               enddo
           enddo
        enddo ! s lopp
      enddo ! w loop
             close(iounit)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else
             write(*,'(A,I2,A)') 'ERROR: Doublecounting type ',DCtype,' not recognized !!!'
             stop 1

          endif ! DCtype

          deallocate( dc_tmp )

       else
          write(*,'(A)') ' Use the DFT Hartree+Fock for the LDA+DMFT Doublecounting...'
          ! This is subtracted already when creating Gloc, so we dont need to do anything
       endif ! if (useSigK==1) then
       write(*,'(A)')

   end subroutine get_sigma_gw_doublecounting

!  ============================================================
!  == Calculate the Doublecounting for the Polarization
!  ============================================================
   subroutine get_pol_doublecounting(p_gw_dc,P_gw,gloc)
      complex(kr),intent(inout) :: p_gw_dc(:,:,:,:)   ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: P_gw(:,:,:,:,:)    ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: gloc(:,:,:,:)      ! norb,norb,nspin,nomega

      integer                 :: w,k,s,m1,m2,a
      real(kr)                :: coeffs(5)
      complex(kr),allocatable :: dc_tmp(:,:,:,:)   ! norb**2,norb**2,nspin,nnu

!!!!!!!!!!!!!
             integer(ki) :: m
             integer(ki),parameter :: iounit = 10
             real(kr)              :: readindata(1+2*nspin*norb**2) ! temparray for readin
             complex(kr)           :: tmp
!!!!!!!!!!!!!

      p_gw_dc = (0.0_kr, 0.0_kr)

      if (usePolK/=0) then

         allocate( dc_tmp(norb**2,norb**2,nspin,nnu) )
         dc_tmp  = (0.0_kr,0.0_kr)

         ! which type of DC?
         if (DCtype==1) then
            ! This is the [GG]_loc doublecounting

            write(*,'(A)') ' Calculate the doublecounting as [GG]_loc ...'
            write(*,'(A)') ' We assume the DC to be orbital-diagonal, since p_imp is ...'
            write(*,'(A)') ' WARNING!!! THIS IS NOT CORRECT FOR TRUE MULTI-ORBITAL!!!!!!!!!!!!!'

            call get_locpart_mats(dc_tmp,P_gw)

            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  m2 = norb_dmft(a,m1)-1
                  p_gw_dc(m2*norb+m2+1,m2*norb+m2+1,:,:) = dc_tmp(m2*norb+m2+1,m2*norb+m2+1,:,:)
               enddo
            enddo

          else if (DCtype==2) then
             write(*,'(A)') 'WARNING: DCtype = 2 has yet to be implemented for the polarization !!!'
             write(*,'(A)') ' the doublecounting p_gw_dc will be zero! !!!'


          else
             write(*,'(A,I2,A)') 'ERROR: Doublecounting type ',DCtype,' not recognized !!!'
             stop 1

          endif ! DCtype

          deallocate( dc_tmp )

       endif ! if (usePolK==1) then
       write(*,'(A)')
   end subroutine get_pol_doublecounting


!  ============================================================
!  == Calculate the local part of a Matsubara Function
!  ============================================================
   subroutine get_locpart_mats(f_loc,f_gw)
      complex(kr),intent(inout) :: f_loc(:,:,:,:)     ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: f_gw(:,:,:,:,:)    ! norb,norb,nspin,nkpts,nomega

      integer :: k,n
      f_loc = (0.0_kr, 0.0_kr)

      n = size(f_gw(1,1,1,1,:))
      do k=1,size(f_gw(1,1,1,:,1))
         f_loc(:,:,:,1:n) = f_loc(:,:,:,1:n) + f_gw(:,:,:,k,:)/size(f_gw(1,1,1,:,1))
      enddo
   end subroutine get_locpart_mats

!  ============================================================
!  == Generate the lattice Green's function at given freq,spin,k-vec
!  ============================================================
   function gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      integer(ki),intent(in)    :: s,k,w               ! spin,kpoint,freq
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

      complex(kr) :: gfk(norb,norb)

      complex(8),allocatable :: gloc_tmp(:,:)  ! norb,norb
      real(kr),allocatable   :: unity(:,:)     ! norb,norb
      integer(ki)            :: i

  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gloc_tmp(norb,norb) )
      allocate( unity(norb,norb) )

!  ====================================================
!  == Create identity matrix for use with iw_n + mu ==
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

!  ====================================================

     ! First create the inverse DFT Green's function
     gloc_tmp = (ci*wn(w-1) + mu + hfield(s) )*unity - h_dft(:,:,s,k) 
     !gloc_tmp = ( -5.0 + w*15.0/nomega  + ci*0.01 + mu)*unity - h_dft(:,:,s,k)

     ! Then add all GW, DC and DMFT contributions
     if (.not. onlyDFT) then
        gloc_tmp = gloc_tmp               &
          &      + dft_hartree(:,:,s)     & ! Subtract the DFT Hartree term (no Hartree and Exchange left)
          &      + dft_exchange(:,:,s)    & ! also subtract exchange
          &      - simp(:,:,s,w)            ! Add the impurity DMFT Selfenergy including Hartree-Fock

       if ( useSigK/=0 )then
        gloc_tmp = gloc_tmp               &
          &      + v_xc(:,:,s,k)          & ! Subtract V_XC from DFT
          &      - s_gw(:,:,s,k,w)        & ! Add the GW ex-cor. Selfenergy (add Exchange back)
          &      + s_gw_dc(:,:,s,w)         ! remove the doublecounting between GW and DMFT (remove local exchange)
!        else
!           ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!           gloc_tmp = gloc_tmp + dft_exchange(:,:,s)   ! Subtract the DFT Exchange term, it's already in simp
        endif
     endif ! (.not. onlyDFT)

!      call check_sym_orb("gloc_tmp before inv",gloc_tmp)

     ! Use this if you want G0 on real frequencies
     !gloc_tmp = ( -5.0+w*15.0/nomega+0.01*ci +mu )*unity - h_dft(:,:,s,k)
     ! Use this if you want G on real frequencies
!     gloc_tmp = ( -10.0+w*13.0/nomega+0.01*ci +mu )*unity - h_dft(:,:,s,k) &
!                                 &    + dft_hartree(:,:,s)+ dft_exchange(:,:,s) - simp(:,:,s,w) 


     ! Then invert with LAPACK
     call ZGETRF(norb, norb, gloc_tmp, norb, ipiv, info_lpck)
     if (info_lpck /= 0) then
        write(*,'(A,I3)') 'ERROR: Greensfunction matrix is numerically singular! Return value', &
                  & info_lpck
        stop 
     end if
     call ZGETRI(norb, gloc_tmp, norb, ipiv, work, norb, info_lpck)
    if (info_lpck /= 0) then
       stop 'Matrix inversion failed!'
    end if
!     call check_sym_orb("gloc_tmp after inv",gloc_tmp)

    gfk = gloc_tmp
!     write(*,'(I6,I3,I6)') k,s,w
!     call check_sym("sigma_k",s_gw(:,:,:,k,:))
!     call check_sym("gfk",gloc)


      deallocate( unity )
      deallocate( gloc_tmp )
   end function gfk

!  ============================================================
!  == Calculate the local Green's function
!  ============================================================
   subroutine get_gloc(gloc,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      complex(kr),intent(inout) :: gloc(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

      integer(ki)            :: w,k,s

!  ====================================================
!  == Do the summation over all k ==
      gloc = (0.0_kr,0.0_kr)

!     Here one should parallelize over frequencies!
      do w=1,nomega
         call progress(w,nomega)
         
         do k=1,nkpts
            do s=1,nspin
              gloc(:,:,s,w) = gloc(:,:,s,w)  &
               &  + gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)/real(nkpts,kind=kr)

            enddo ! s
         enddo ! k
      enddo ! 
!  ====================================================
   end subroutine get_gloc



!  ============================================================
!  == Calculate the filling of a Green's function
!  ============================================================
   function filling(gf)
      real(kr)               :: filling
      complex(kr),intent(in) :: gf(:)    ! nomega
      integer(ki)            :: w
      real(kr)               :: coeffs(5)
      
      filling = 0.0_kr
      ! Use high-frequency correction up to 4-th order      

      coeffs = get_highfreq_coeff(gf,1)
      filling = 0.5_kr - coeffs(3)*beta/4.0_kr + coeffs(5)*(beta**3)/48.0_kr
      do w=1,nomega
         filling = filling + (2.0_kr/beta)*( real(gf(w))        &
                                 &   + coeffs(3)/wn(w-1)**2     &
                                 &   - coeffs(5)/wn(w-1)**4 )
      enddo        
   end function filling
   
!  ============================================================
!  == Calculate the total charge of a Green's function
!  ============================================================
   function total_filling(gf)
      complex(kr),intent(in) :: gf(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr)                  :: total_filling
      integer                   :: s,m
      
      total_filling = 0.0_kr
      do s=1,nspin
         do m=1,norb
            total_filling = total_filling + filling(gf(m,m,s,:))
         enddo
      enddo   
   end function total_filling

!  ============================================================
!  == Return the density matrix
!  ============================================================
   subroutine density_matrix(density,gf)
      complex(kr),intent(inout) :: density(:,:)   ! norb,norb
      complex(kr),intent(in)    :: gf(:,:,:)      ! norb,norb,nomega

      integer                :: m1,m2,w
      real(kr)               :: coeffs(5)
      complex(kr)            :: fill
      
      density = (0.0_kr,0.0_kr)
         do m1=1,norb
            ! Diagonal
            density(m1,m1) = filling(gf(m1,m1,:))

            ! Offdiagonal
            do m2=1,norb
               coeffs = get_highfreq_coeff(gf(m1,m2,:),0)
               fill = 0.5_kr*coeffs(2) - coeffs(3)*beta/4.0_kr + coeffs(5)*(beta**3)/48.0_kr
               do w=1,nomega
                  fill = fill + (1.0_kr/beta)*( gf(m1,m2,w) + conjg(gf(m2,m1,w))      &
                                          &   + 2*coeffs(3)/wn(w-1)**2                     &
                                          &   - 2*coeffs(5)/wn(w-1)**4 )
               enddo ! w
               density(m1,m2) = fill
            
            enddo ! m2
         enddo ! m1
   end subroutine density_matrix

!  ============================================================
!  == Calculate the Weiss field and local orbital levels
!  == For the non-DMFT orbitals Gbath will be Gloc
!  ============================================================
   subroutine get_weiss_field_loc_levels_dyson(gbath,mu_loc,gloc,simp)
      complex(kr),intent(inout) :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr),intent(inout)    :: mu_loc(:,:,:)       ! norb,norb,nspin
      complex(kr),intent(in)    :: gloc(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega

      complex(8),allocatable  :: gloc_tmp(:,:)        ! norb,norb BUT reduced to onsite
      complex(8),allocatable  :: simp_tmp(:,:)        ! norb,norb BUT reduced to onsite
      ! use this for loc.level calculation
      complex(kr),allocatable :: gbath_inv(:,:,:,:)   ! norb,norb,nspin,nomega
      real(kr)                :: coeffs(5)  !,coeffsGinv(5),coeffsS(5)
      integer(ki)             :: w,s,i,j,a,m1,m2,mm1,mm2
  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gloc_tmp(norb,norb) )
      allocate( simp_tmp(norb,norb) )
      allocate( gbath_inv(norb,norb,nspin,nomega) )
      gbath =  (0.0_kr,0.0_kr)
      mu_loc = (0.0_kr)      

      do w=1,nomega
         do s=1,nspin

     ! Remove all inter-atom terms   ! Take only the submatrix of the correlated orbitals
            gloc_tmp =  (0.0_kr,0.0_kr)
            simp_tmp =  (0.0_kr,0.0_kr)

            ! use failsafe only the diagonals as default
            do m1=1,norb
               gloc_tmp(m1,m1) = gloc(m1,m1,s,w)
               simp_tmp(m1,m1) = simp(m1,m1,s,w)
            enddo

            ! copy only remaining intra-atom terms of the dmft orbitals
            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  do m2=0,norbPerAtom(a)
                     mm1 = norb_dmft(a,m1)
                     mm2 = norb_dmft(a,m2)
                     gloc_tmp(mm1,mm2) = gloc(mm1,mm2,s,w)
                     simp_tmp(mm1,mm2) = simp(mm1,mm2,s,w)
                 enddo ! m2
               enddo ! m1
            enddo ! a

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gloc_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: local Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gloc_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            
            ! then get gbath^-1 by Dyson equation: gloc^-1 + sigma
            gloc_tmp = gloc_tmp + simp_tmp

            ! now gloc_tmp is equal to gbath^-1 !
            
            ! save for later for the local orbital levels
            gbath_inv(:,:,s,w) = gloc_tmp

            ! Then invert gbath^-1 !
            call ZGETRF(norb, norb, gloc_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: local Weiss field is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gloc_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
         
             gbath(:,:,s,w) = gloc_tmp
         enddo
      enddo

      
      ! now get the local orbital levels
      do s=1,nspin
         do a=1,noAtoms
            do m1=1,norbPerAtom(a)
               do m2=1,norbPerAtom(a)
                  i = norb_dmft(a,m1)
                  j = norb_dmft(a,m2)
                  coeffs = get_highfreq_coeff(gbath_inv(i,j,s,:),0)
                  mu_loc(i,j,s) = coeffs(1)

            !coeffs = get_highfreq_coeff(gloc(i,i,s,:),1)
            !mu_loc(i,s) = dft_hartree(i,i,s) + dft_exchange(i,i,s) - coeffs(3)
               enddo ! m2
            enddo ! m1
         enddo ! a
      enddo ! s


      deallocate( gloc_tmp )
      deallocate( simp_tmp )
      deallocate( gbath_inv )
   end subroutine get_weiss_field_loc_levels_dyson


!  ============================================================
!  == Calculate the Weiss field and local orbital levels
!  == With the general formula allowing for nonlocal Selfenergies
!  == For the non-DMFT orbitals Gbath will be Gloc
!  ============================================================
   subroutine get_weiss_field_loc_levels(gbath,mu_loc,gloc,  &
                     & h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      complex(kr),intent(inout) :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr),intent(inout)    :: mu_loc(:,:,:)       ! norb,norb,nspin
      complex(kr),intent(in)    :: gloc(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      complex(8),allocatable  :: g_tmp(:,:)       ! norb,norb
      real(kr),allocatable    :: unity(:,:)       ! norb,norb
      complex(kr),allocatable :: epsGeps(:,:)     ! norb,norb
      complex(kr),allocatable :: epsG(:,:)        ! norb,norb
      complex(kr),allocatable :: Geps(:,:)        ! norb,norb
      complex(kr),allocatable :: Glocinv(:,:,:,:) ! norb,norb,nspin,nomega
      complex(kr),allocatable :: epsGeps_proj(:,:) ! norb,norb
      complex(kr),allocatable :: epsG_proj(:,:)    ! norb,norb
      complex(kr),allocatable :: Geps_proj(:,:)    ! norb,norb
      real(kr),allocatable    :: mu_loc_proj(:,:)

      integer(ki)            :: w,s,i,j,k,m1,m2,mm1,mm2, a

      allocate( g_tmp(norb,norb) )
      allocate( unity(norb,norb) )
      allocate( mu_loc_proj(norb,norb) )
      allocate( epsGeps(norb,norb) )
      allocate( epsG(norb,norb) )
      allocate( Geps(norb,norb) )
      allocate( Glocinv(norb,norb,nspin,nomega) )
      allocate( epsGeps_proj(norb,norb) )
      allocate( epsG_proj(norb,norb) )
      allocate( Geps_proj(norb,norb) )
      g_tmp =  (0.0_kr,0.0_kr)
      gbath =  (0.0_kr,0.0_kr)
      mu_loc = (0.0_kr)     
      mu_loc_proj = (0.0_kr)     
      epsGeps =  (0.0_kr,0.0_kr)
      epsG =  (0.0_kr,0.0_kr)
      Geps =  (0.0_kr,0.0_kr)
      Glocinv =  (0.0_kr,0.0_kr)
      epsGeps_proj =  (0.0_kr,0.0_kr)
      epsG_proj =  (0.0_kr,0.0_kr)
      Geps_proj =  (0.0_kr,0.0_kr)

!  ====================================================
!  == Create identity matrix for use with iw_n + mu ==
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

      ! calculate mu_loc here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,'(A)') ' Calculate the local orbital levels...'
      mu_loc = (0.0_kr)
      do k=1,nkpts
         do s=1,nspin
            g_tmp = ( mu + hfield(s) )*unity - h_dft(:,:,s,k) 

            ! Then add all GW, DC and DMFT contributions
            if (.not. onlyDFT) then
               g_tmp = g_tmp                 &
   &                    + dft_hartree(:,:,s) &     ! Subtract the DFT Hartree term (no Hartree and Exchange left)
   &                    + dft_exchange(:,:,s)   ! Subtract the DFT Exchange term

               if ( useSigK/=0 )then
                  g_tmp = g_tmp     + v_xc(:,:,s,k)       ! Subtract V_XC from DFT to get the true noninteracting H0
!               else
!                  ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                  g_tmp = g_tmp + dft_exchange(:,:,s)   ! Subtract the DFT Exchange term, it's already in simp
               endif
            endif ! (.not. onlyDFT)
            mu_loc(:,:,s) = mu_loc(:,:,s) + g_tmp/nkpts
         enddo ! s
      enddo ! k


      ! calculate Glocinv !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do s=1,nspin
         do w=1,nomega
!            g_tmp = unity
            g_tmp = (0.0_kr,0.0_kr)

            !use failsafe only the diagonals as default
            do m1=1,norb
               g_tmp(m1,m1) = gloc(m1,m1,s,w)
            enddo

            ! copy only remaining intra-atom terms of the dmft orbitals
            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  do m2=1,norbPerAtom(a)
                     mm1 = norb_dmft(a,m1)
                     mm2 = norb_dmft(a,m2)
                     g_tmp(mm1,mm2) = gloc(mm1,mm2,s,w)
                 enddo ! m2
               enddo ! m1
            enddo ! a


            ! This is missing the other equivalent atoms
            !do m1=1,norbPerAtom
            !   do m2=1,norbPerAtom
            !      mm1 = dmftOrbsIndex(m1)
            !      mm2 = dmftOrbsIndex(m2)
            !      g_tmp(mm1,mm2) = gloc(mm1,mm2,s,w)
            !   enddo ! m2
            !enddo ! m1
 
            call ZGETRF(norb, norb, g_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,'(A,I3)') 'ERROR: Greensfunction matrix is numerically singular! Return value', &
                                                   & info_lpck
               stop 
            end if
            call ZGETRI(norb, g_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
               stop 'Matrix inversion failed!'
            end if
            Glocinv(:,:,s,w) = g_tmp
         enddo ! w
      enddo ! s

      ! calculate Bath green's function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,'(A)') ' Calculate general bath...'
      do w=1,nomega
         call progress(w,nomega)

         do s=1,nspin
            epsGeps = (0.0_kr,0.0_kr)
            epsG = (0.0_kr,0.0_kr)
            Geps = (0.0_kr,0.0_kr)
            g_tmp = (0.0_kr,0.0_kr)

            do k=1,nkpts
               g_tmp = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)

               ! First create the inverse DFT Green's function
!               g_tmp = (ci*wn(w-1) + mu + hfield(s) )*unity - h_dft(:,:,s,k) 
!
!               ! Then add all GW, DC and DMFT contributions
!               if (.not. onlyDFT) then
!                  g_tmp = g_tmp               &
!      &                 + dft_hartree(:,:,s)     & ! Subtract the DFT Hartree term (no Hartree and Exchange left)
!      &                 - simp(:,:,s,w)            ! Add the impurity DMFT Selfenergy including Hartree-Fock
!
!                  if ( useSigK==1 )then
!                     g_tmp = g_tmp               &
!      &                    + v_xc(:,:,s,k)          & ! Subtract V_XC from DFT
!      &                    - s_gw(:,:,s,k,w)        & ! Add the GW ex-cor. Selfenergy (add Exchange back)
!      &                    + s_gw_dc(:,:,s,w)         ! remove the doublecounting between GW and DMFT (remove local exchange)
!                  else
!                     ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                     g_tmp = g_tmp + dft_exchange(:,:,s)   ! Subtract the DFT Exchange term, it's already in simp
!                  endif
!               endif ! (.not. onlyDFT)
!
!               ! Then invert with LAPACK
!               call ZGETRF(norb, norb, g_tmp, norb, ipiv, info_lpck)
!               if (info_lpck /= 0) then
!                  write(*,'(A,I3)') 'ERROR: Greensfunction matrix is numerically singular! Return value', &
!                            & info_lpck
!                  stop 
!               end if
!               call ZGETRI(norb, g_tmp, norb, ipiv, work, norb, info_lpck)
!              if (info_lpck /= 0) then
!                 stop 'Matrix inversion failed!'
!              end if

              ! Then add up for every k-point the different combinations
              epsGeps = epsGeps + matmul(matmul( h_dft(:,:,s,k), g_tmp), h_dft(:,:,s,k) )/nkpts
              epsG = epsG + matmul( h_dft(:,:,s,k), g_tmp)/nkpts
              Geps = Geps + matmul( g_tmp, h_dft(:,:,s,k) )/nkpts

!              epsGeps = epsGeps + matmul(matmul( h_dft(:,:,s,k)-v_xc(:,:,s,k), g_tmp), h_dft(:,:,s,k)-v_xc(:,:,s,k) )/nkpts
!              epsG = epsG + matmul( h_dft(:,:,s,k)-v_xc(:,:,s,k), g_tmp)/nkpts
!              Geps = Geps + matmul( g_tmp, h_dft(:,:,s,k)-v_xc(:,:,s,k) )/nkpts

            enddo ! k

            epsGeps_proj = (0.0_kr,0.0_kr)
            Geps_proj    = (0.0_kr,0.0_kr)
            epsG_proj    = (0.0_kr,0.0_kr)
            mu_loc_proj  = (0.0_kr,0.0_kr)
            !use failsafe only the diagonals as default
            do m1=1,norb
               epsGeps_proj(m1,m1) = epsGeps(m1,m1)
               Geps_proj(   m1,m1) = Geps(   m1,m1) 
               epsG_proj(   m1,m1) = epsG(   m1,m1)
               mu_loc_proj( m1,m1) = mu_loc( m1,m1,s)
            enddo

            ! copy only remaining intra-atom terms of the dmft orbitals
            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  do m2=1,norbPerAtom(a)
                     mm1 = norb_dmft(a,m1)
                     mm2 = norb_dmft(a,m2)

                     epsGeps_proj(mm1,mm2) = epsGeps(mm1,mm2)
                     Geps_proj(mm1,mm2) = Geps(mm1,mm2) 
                     epsG_proj(mm1,mm2) = epsG(mm1,mm2)
                     mu_loc_proj(mm1,mm2) = mu_loc(mm1,mm2,s)
                 enddo ! m2
               enddo ! m1
            enddo ! a


            ! This is missing the other equivalent atoms
            ! Now project on the first equiv. atom and DMFT orbitals
!            epsGeps_proj = unity
!            Geps_proj = unity
!            epsG_proj = unity
!            mu_loc_proj = unity
!            do m1=1,norbPerAtom
!               do m2=1,norbPerAtom
!                  mm1 = dmftOrbsIndex(m1)
!                  mm2 = dmftOrbsIndex(m2)
!                  epsGeps_proj(mm1,mm2) = epsGeps(mm1,mm2)
!                  Geps_proj(mm1,mm2) = Geps(mm1,mm2) 
!                  epsG_proj(mm1,mm2) = epsG(mm1,mm2)
!                  mu_loc_proj(mm1,mm2) = mu_loc(mm1,mm2,s)
!               enddo ! m2
!            enddo ! m1

            ! now we can build Gbath^-1
            g_tmp = ci*wn(w-1)*unity + mu_loc_proj - ( epsGeps_proj - matmul(matmul(epsG_proj,Glocinv(:,:,s,w)),Geps_proj) )

            ! Then invert gbath^-1 !
            call ZGETRF(norb, norb, g_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: local Weiss field is numerically singular! Return value',&
                            & info_lpck
               stop 
            end if
            call ZGETRI(norb, g_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
               stop 'Matrix inversion failed!'
            end if
         
            gbath(:,:,s,w) = g_tmp

         enddo ! s
      enddo ! w

!      if (usefreqU==1) then
!         write(*,'(A,F8.3)') " We shift the impurity levels AFTER GENERATING Gbath by +K'(0) =", Kp0
!         do m1=1,norb_dmft
!            m2 = dmftOrbsIndex(m1)
!            do s=1,nspin
!               mu_loc(m2,m2,s) = mu_loc(m2,m2,s) + Kp0
!            enddo
!         enddo
!      endif


      deallocate( g_tmp )
      deallocate( unity )
      deallocate( mu_loc_proj )
      deallocate( epsGeps )
      deallocate( epsG )
      deallocate( Geps )
      deallocate( glocinv )
      deallocate( epsGeps_proj )
      deallocate( epsG_proj )
      deallocate( Geps_proj )

   end subroutine get_weiss_field_loc_levels

!  ============================================================
!  == Calculate the effective local interaction U(w)
!  ============================================================
   subroutine get_uloc(uloc,P_gw,pimp,p_gw_dc,W_gw,V_gw,bare)
      complex(kr),intent(inout) :: uloc(:,:,:,:)      ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: P_gw(:,:,:,:,:)    ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: pimp(:,:,:,:)      ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: p_gw_dc(:,:,:,:)   ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)    :: W_gw(:,:,:,:,:)    ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(in)    :: V_gw(:,:,:,:,:)    ! norb**2,norb**2,nspin,nkpts,nnu
      logical,intent(in)    :: bare               ! just get the bare U

      real(kr)    :: unity(norb**2,norb**2)
      complex(kr) :: udensw(nnu)
      real(kr)    :: coeffs(2)
      complex(kr) :: wloc(norb**2,norb**2)   , wproj(norb**2,norb**2)     
      complex(kr) :: ploc(norb**2,norb**2)   , pproj(norb**2,norb**2)    
      complex(kr) :: pwloc(norb**2,norb**2)  , pwproj(norb**2,norb**2)      
      complex(kr) :: wploc(norb**2,norb**2)  , wpproj(norb**2,norb**2)      
      complex(kr) :: pwploc(norb**2,norb**2) , pwpproj(norb**2,norb**2)       
      complex(kr) :: tmp(norb**2,norb**2) 

      integer(ki) :: s,i,j,k,l,a,b,c,d,n,nwV, aa
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(2) ! temparray for readin

!  ====================================================
!  == Create identity matrix
      unity = (0.0_kr)
      do i=1,norb**2
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================
      uloc = (0.0_kr,0.0_kr)

      if (bare .or. (usePolK==0 .and. useUw==0) ) then
         write(*,'(A)') ' Use the bare static interaction for U(w)...'
         nwV = size(V_gw(1,1,1,1,:))

         do n=1,size(uloc(1,1,1,:))
            do k=1,size(V_gw(1,1,1,:,n))
               uloc(:,:,:,n) = uloc(:,:,:,n) + V_gw(:,:,:,k,nwV)/size(V_gw(1,1,1,:,nwV))
            enddo ! k
         enddo ! n

      else if (useUw==1) then
         write(*,'(A)') ' Read U(w) from file udensw.dat...'
         !!!!!!!!!!!!!!
         udensw = (0.0_kr,0.0_kr)
         open(unit=iounit,file="udensw.dat",status="old")
         do n=1,size(udensw)
            read(iounit,*) readindata  ! frequency, U(w) value
            udensw(n) = readindata(2)  ! this function has no imaginary part
         enddo ! w loop
         close(iounit)
         !!!!!!!!!!!

         ! Now we have the frequency dependent density part, construct the matrix Uloc(w)
         coeffs = get_highfreq_coeff_bosonic(udensw)
         nwV = size(V_gw(1,1,1,1,:))

         do n=1,size(uloc(1,1,1,:))           ! just take the bare part from Vbare
            do k=1,size(V_gw(1,1,1,:,n))
               uloc(:,:,:,n) = uloc(:,:,:,n) + V_gw(:,:,:,k,nwV)/size(V_gw(1,1,1,:,nwV))
            enddo ! k

            ! now we screen the density-density terms
            do aa=1,noAtoms
               do k=1,norbPerAtom(aa)
                  do l=1,norbPerAtom(aa)
                     a = norb_dmft(aa,k)-1
                     b = norb_dmft(aa,l)-1
                     uloc(a*norb+a+1, b*norb+b+1,:,n) = uloc(a*norb+a+1, b*norb+b+1,:,n) - (coeffs(1)-udensw(n))
                  enddo
               enddo 
            enddo ! a

         enddo ! n

      else if (usePolK/=0) then
         write(*,'(A)') ' Calculate effective interaction U(w) from P and W...'
         do n=1,nnu
            call progress(n,nnu)

            do s=1,nspin
               wloc = (0.0_kr,0.0_kr)
               ploc = (0.0_kr,0.0_kr)
               pwloc = (0.0_kr,0.0_kr)
               wploc = (0.0_kr,0.0_kr)
               pwploc = (0.0_kr,0.0_kr)
   
               do k=1,nkpts
                  ! This is the corrected polarizaton
                  tmp = P_gw(:,:,s,k,n) + pimp(:,:,s,n) - p_gw_dc(:,:,s,n)

                  wloc   =   wloc + W_gw(:,:,s,k,n)/nkpts
                  ploc   =   ploc + tmp/nkpts
                  pwloc  =  pwloc + matmul(tmp,W_gw(:,:,s,k,n))/nkpts
                  wploc  =  wploc + matmul(W_gw(:,:,s,k,n),tmp)/nkpts
                  pwploc = pwploc + matmul(tmp,matmul(W_gw(:,:,s,k,n),tmp))/nkpts
               enddo ! k

!               if (norb/=norbPerAtom) then
!                  stop " ERROR: Calculation of selfconsistent U(w) for norb/=norbPerAtom not implemented!!!"
!               endif

               ! Now project on the first equiv. atom and DMFT orbitals
               wproj = (0.0_kr,0.0_kr)
               pproj = unity  ! All others can be zero because they are not inverted
               pwproj = (0.0_kr,0.0_kr)
               wpproj = (0.0_kr,0.0_kr)
               pwpproj = (0.0_kr,0.0_kr)
               do aa=1,noAtoms
                  do i=1,norbPerAtom(aa)
                     do j=1,norbPerAtom(aa)
                        do k=1,norbPerAtom(aa)
                           do l=1,norbPerAtom(aa)
                              a = norb_dmft(aa,i)-1
                              b = norb_dmft(aa,j)-1
                              c = norb_dmft(aa,k)-1
                              d = norb_dmft(aa,l)-1
   
                              wproj(   a*norb+c+1, b*norb+d+1 ) =   wloc(a*norb+c+1, b*norb+d+1)
                              pproj(   a*norb+c+1, b*norb+d+1 ) =   ploc(a*norb+c+1, b*norb+d+1)
                              pwproj(  a*norb+c+1, b*norb+d+1 ) =  pwloc(a*norb+c+1, b*norb+d+1)
                              wpproj(  a*norb+c+1, b*norb+d+1 ) =  wploc(a*norb+c+1, b*norb+d+1)
                              pwpproj( a*norb+c+1, b*norb+d+1 ) = pwploc(a*norb+c+1, b*norb+d+1)                           
                           enddo ! l
                        enddo ! k
                     enddo ! j
                  enddo ! i
               enddo ! aa 

               ! Then just copy back
               wloc = wproj
               ploc = pproj !+ 0.0001*unity  ! the inversion can be a bit tricky
               pwloc = pwproj
               wploc = wpproj
               pwploc = pwpproj
   
               ! include correction term for causal result
               if (bathMethod==0) then
                  tmp = ploc + pwploc  
                  call get_inverse(tmp, tmp, norb**2)
                  uloc(:,:,s,n) = wloc - matmul(pwloc,matmul(tmp,wploc))
               else
                  ! standard Dyson U = W * (1 + PV)^-1
                  tmp = unity + matmul(ploc,wloc)
                  call get_inverse(tmp, tmp, norb**2)
                  uloc(:,:,s,n) = matmul(wloc,tmp)
               endif
   
            enddo ! s
         enddo ! n

      else
         stop 'Combination of usePolK and useUw not possible in get_uloc'
      endif
   end subroutine get_uloc


!  ============================================================
!  == Calculate the Selfenergy from the Dyson equation
!  ============================================================
   subroutine get_sigma_dyson(sigma,gimp,gbath)
      complex(kr),intent(inout) :: sigma(:,:,:,:)     ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: gimp(:,:,:,:)      ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega

      complex(8),allocatable  :: gbath_tmp(:,:)       ! norb,norb
      complex(8),allocatable  :: gimp_tmp(:,:)        ! norb,norb
      complex(8),allocatable  :: simp_tmp(:,:)        ! norb,norb
      integer(ki)             :: w,s,i,m,m1,m2,mm1,mm2,a
      real(kr)                :: coeffs_gb(norb,2,5),coeffs_gi(norb,2,5), lincoeff
  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gbath_tmp(norb,norb) )
      allocate( gimp_tmp(norb,norb) )
      allocate( simp_tmp(norb,norb) )
      

      write(*,'(A)') 'Obtain Selfenergy from Dyson equation...'


      write(*,'(A)') 'We only read the impurity Self-energy for the DMFT orbitals, others are not modified'
!      sigma = (0.0_kr,0.0_kr)

      ! When reading in Gimp we already copied the imp.GF from the
      ! first to all other equiv. atoms, i.e. only the onsite
      ! terms of Gbath and Gimp are nonzero
      
      ! STILL, let's restrict to only the first atom to obtain the Selfenergy,
      ! since this is what went into the impurity solver 

      coeffs_gb = (0.0_kr)
      coeffs_gi = (0.0_kr)
      lincoeff = (0.0_kr)
      do a=1,noAtoms
         do m=1,norbPerAtom(a)
            do s=1,nspin
               m1 = norb_dmft(a,m)
               coeffs_gb(m1,s,:) = get_highfreq_coeff(gbath(m1,m1,s,:),0)
               coeffs_gi(m1,s,:) = get_highfreq_coeff(gimp(m1,m1,s,:),0)
            enddo
         enddo
      enddo

      do w=1,nomega
         do s=1,nspin
            gbath_tmp = (0.0_kr,0.0_kr)
            gimp_tmp = (0.0_kr,0.0_kr)
            simp_tmp = (0.0_kr,0.0_kr)

            do m1=1,norb  ! Failsafe
               gbath_tmp(m1,m1) = (1.0_kr,0.0_kr)
               gimp_tmp(m1,m1) = (1.0_kr,0.0_kr)
            enddo ! m1

            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  do m2=1,norbPerAtom(a)
                     mm1 = norb_dmft(a,m1)
                     mm2 = norb_dmft(a,m2)
                     gbath_tmp(mm1,mm2) = gbath(mm1,mm2 ,s,w)
                     gimp_tmp(mm1,mm2)  = gimp(mm1,mm2 ,s,w)
                  enddo ! m2
               enddo ! m1
            enddo ! a

            ! Fix tails
            do a=1,noAtoms
               do m=1,norbPerAtom(a)
                  mm1 = norb_dmft(a,m)
!                  gbath_tmp(mm1,mm1) = gbath_tmp(mm1,mm1)/coeffs_gb(mm1,s,2)  ! Turns out it doesn't work so well!
!                  gimp_tmp(mm1,mm1)  = gimp_tmp(mm1,mm1) /coeffs_gi(mm1,s,2)
                  if (w==1) then
                     if (abs(coeffs_gb(mm1,s,2)-1.0_kr) > 0.05_kr) then 
                        write(*,'(A,I3,A,I2,A,F9.5)') 'ERROR: Gbath tail: m=',mm1,", s=",s,' = ', coeffs_gb(mm1,s,2)
                        stop(1)
                     endif
                     if (abs(coeffs_gi(mm1,s,2)-1.0_kr) > 0.05_kr) then 
                        write(*,'(A,I3,A,I2,A,F9.5)') 'ERROR: Gimp tail: m=',mm1,", s=",s,' = ', coeffs_gi(mm1,s,2)
                        stop(1)
                     endif
                  endif
               enddo ! m
            enddo ! a
            !if (w==1) then
            !   write(*,*) 'Gbath = ',gbath_tmp
            !   write(*,*) 'Gimp = ',gimp_tmp
            !endif

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gbath_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gbath_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            ! gbath_tmp is now the inverse of g_bath

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gimp_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: impurity Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gimp_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            ! gimp_tmp is now the inverse of g_imp
            
            ! Calculate selfenergy
            simp_tmp = gbath_tmp - gimp_tmp

!            ! consistency check
!            do m1=1,norbPerAtom
!               if ( AIMAG( simp_tmp(m1,m1)) > 0.0_kr ) then
!                  simp_tmp(m1,m1) = REAL( simp_tmp(m1,m1) )
!               end if
!            end do

            ! copy to full Selfenergy but don't overwrite other stuff
            do a=1,noAtoms
               do m1=1,norbPerAtom(a)
                  do m2=1,norbPerAtom(a)
                     mm1 = norb_dmft(a,m1)
                     mm2 = norb_dmft(a,m2)
                     sigma(mm1,mm2,s,w) = simp_tmp(mm1,mm2)

!                     if (a==0) then
!                        sigma(mm1,mm2,s,w) = simp_tmp(m1+1, m2+1)
!                     else
!                        sigma(mm1,mm2,3-s,w) = simp_tmp(m1+1, m2+1)
!                     endif

                  enddo ! m2
               enddo ! m1
            enddo ! a

         enddo ! s
      enddo ! w

      ! Do final tail check and correct if there's a linar component
!      do s=1,nspin
!         do m1=1,norb_dmft
!            mm1 = dmftOrbsIndex(m1)
!
!            lincoeff = get_linear_imag_coeff( sigma(mm1,mm1,s,:) )
!
!            do w=1,nomega
!               sigma(mm1,mm1,s,w) = sigma(mm1,mm1,s,w) - lincoeff*ci*(2*w-1)*pi/beta
!            enddo ! w
!
!         enddo ! m1
!      enddo ! s

      deallocate(gbath_tmp)
      deallocate(gimp_tmp)
      deallocate(simp_tmp)
   end subroutine get_sigma_dyson
   
!  ============================================================
!  == Copy the selfenergy to other atoms for AFM order ! HIGHLY HARDCODED !!!!!!!
!  ============================================================
   subroutine copy_sigma_afm(f)
      complex(kr),intent(inout) :: f(:,:,:,:)     ! norb,norb,nspin,nomega

      integer(ki)             :: s, m1,m2
         ! now copy to the other equivalent atoms

         ! ONSITE TERMS
         do m1=1,2
            do m2=1,2
               do s=1,2

                  ! A site            
                  f( m1+6,m2+6,s,:)   = f( m1,m2 ,s,:)
                  f( m1+10,m2+10,s,:) = f( m1,m2 ,s,:)
                  f( m1+12,m2+12,s,:) = f( m1,m2 ,s,:)

                  ! B site
                  f( m1+4,m2+4,s,:)   = f( m1+2,m2+2 ,s,:)
                  f( m1+8,m2+8,s,:)   = f( m1+2,m2+2 ,s,:)
                  f( m1+14,m2+14,s,:) = f( m1+2,m2+2 ,s,:)

               enddo
            enddo 
         enddo 

         ! Intersite TERMS
         do m1=1,2
            do m2=1,2
               do s=1,2

                  f( m1+0,m2+4,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+4,m2+0,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+0,m2+8,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+8,m2+0,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+6,m2+2,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+2,m2+6,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+10,m2+2,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+2,m2+10,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+6,m2+4,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+4,m2+6,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+12,m2+4,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+4,m2+12,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+6,m2+14,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+14,m2+6,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+10,m2+8,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+8,m2+10,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+12,m2+8,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+8,m2+12,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+10,m2+14,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+14,m2+10,s,:) = f( m1+2,m2+0 ,s,:)

                  f( m1+12,m2+14,s,:) = f( m1+0,m2+2 ,s,:)
                  f( m1+14,m2+12,s,:) = f( m1+2,m2+0 ,s,:)

               enddo
            enddo 
         enddo 

   end subroutine copy_sigma_afm

!  ============================================================
!  == Reperiodize sigma - HIGHLY HARDCODED !!!
!  ============================================================
   subroutine reperiod_sigma(s_gw,simp,h_dft,v_xc,s_gw_dc,dft_hartree,dft_exchange)
      complex(kr),intent(inout) :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(inout) :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin


      integer(ki)               :: ikx,iky,i, ik, m1,m2, n, s
      real(kr)                  :: kx,ky
      complex(kr),allocatable   :: sloc(:,:,:,:)       !norb,norb,nspin,nomega
      complex(kr)               :: gtmp(norb,norb), gunf(norb,norb), g0unf(norb,norb)
      complex(kr)               :: U(norb,norb), Ud(norb,norb)

      allocate( sloc(norb,norb,nspin,nomega) )

      write(*,*) ''
      write(*,*) '!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!!! REPERIODIZE THE SELFENERGY! HAVE YOU CHECKED THE TRAFO? !!!!!'
      write(*,*) '!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) ''

      s_gw = (0.0_kr,0.0_kr)
      sloc = (0.0_kr,0.0_kr)

      do ikx = 0,nkx-1
         do iky = 0,nky-1

            kx = ikx*2*pi/nkx
            ky = iky*2*pi/nky
            ik = ikx*nky+iky+1


                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

                  ! dimer in x-direction FOR 1D MODEL , only take nonlocal parts
               if (1==1) then
                  s_gw(1,2,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*kx) )
                  s_gw(2,1,:,ik,:) = simp(2,1,:,:)*( 0.0_kr + exp(-ci*kx) )
              endif

              if (1==0) then ! Dimer for GF reperiodization
                 U(1,1) = 1/sqrt(2.0)
                 U(2,2) = 1/sqrt(2.0)
                 U(1,2) = -exp(+ci*kx/2) /sqrt(2.0)
                 U(2,1) = +exp(-ci*kx/2) /sqrt(2.0)

                 Ud(1,1) = 1/sqrt(2.0)
                 Ud(2,2) = 1/sqrt(2.0)
                 Ud(1,2) = +exp(+ci*kx/2) /sqrt(2.0)
                 Ud(2,1) = -exp(-ci*kx/2) /sqrt(2.0)

                 do n=1,nomega
                  do s=1,nspin
                    gtmp = gfk(s,ik,n,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,.false.)
                    gunf = matmul( Ud, matmul(gtmp, U) )

                    ! construct G0
                    call get_inverse(gtmp,gtmp,norb)  ! G^-1
                    gtmp = gtmp + simp(:,:,s,n)       ! G0^-1
                    call get_inverse(gtmp,gtmp,norb)  ! G0
                    g0unf = matmul( Ud, matmul(gtmp, U) )

                    gtmp = (0.0_kr,0.0_kr)
                    gtmp(1,1) = 1.0_kr/g0unf(1,1) - 1.0_kr/gunf(1,1)
                    gtmp(2,2) = 1.0_kr/g0unf(2,2) - 1.0_kr/gunf(2,2)

                    s_gw(:,:,s,ik,n) = matmul( U, matmul(gtmp, Ud) )

                   enddo
                 enddo

                 sloc = sloc + s_gw(:,:,:,ik,:)/nkpts
              endif

                  ! 3site in x-direction FOR 1D MODEL , only take nonlocal parts
                !  s_gw(1,3,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*kx) )
                !  s_gw(3,1,:,ik,:) = simp(2,1,:,:)*( 0.0_kr + exp(-ci*kx) )

                  ! 4site in x-direction FOR 1D MODEL , only take nonlocal parts
!                  s_gw(1,4,:,ik,:) = simp(2,3,:,:)*( 0.0_kr + exp(+ci*kx) )
!                  s_gw(4,1,:,ik,:) = simp(3,2,:,:)*( 0.0_kr + exp(-ci*kx) )
!                  simp(1,2,:,:) = simp(2,3,:,:)
!                  simp(2,1,:,:) = simp(3,2,:,:)
!                  simp(3,4,:,:) = simp(2,3,:,:)
!                  simp(4,3,:,:) = simp(3,2,:,:)
!                  simp(1,1,:,:) = simp(2,2,:,:)
!                  simp(4,4,:,:) = simp(2,2,:,:)

                  ! 5site in x-direction FOR 1D MODEL , only take nonlocal parts
!                  s_gw(1,5,:,ik,:) = simp(2,3,:,:)*( 0.0_kr + exp(+ci*kx) )
!                  s_gw(5,1,:,ik,:) = simp(3,2,:,:)*( 0.0_kr + exp(-ci*kx) )
!                  simp(1,2,:,:) = simp(2,3,:,:)
!                  simp(2,1,:,:) = simp(3,2,:,:)
!                  simp(4,5,:,:) = simp(2,3,:,:)
!                  simp(5,4,:,:) = simp(3,2,:,:)
!                  simp(1,1,:,:) = simp(3,3,:,:)
!                  simp(2,2,:,:) = simp(3,3,:,:)
!                  simp(4,4,:,:) = simp(3,3,:,:)
!                  simp(5,5,:,:) = simp(3,3,:,:)

                   ! Just copy the end terms for 1D cluster
!                   s_gw(1,norb,:,ik,:) = simp(1,2,:,:)     *exp(+ci*kx)
!                   s_gw(norb,1,:,ik,:) = simp(norb,norb-1,:,:)*exp(-ci*kx)


                  ! dimer in xy-direction
                  !s_gw(1,1,:,ik,:) = 0.0
                  !s_gw(2,2,:,ik,:) = 0.0
                  !s_gw(1,2,:,ik,:) = (exp(+ci*kx) + exp(+ci*ky) + exp(+ci*(kx+ky)) ) * simp(1,2,:,:)
                  !s_gw(2,1,:,ik,:) = (exp(-ci*kx) + exp(-ci*ky) + exp(-ci*(kx+ky)) ) * simp(2,1,:,:)

               ! dimer in x-direction, only take nonlocal parts
               if (1==0) then
                  s_gw(1,1,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*ky) + exp(-ci*ky) )
                  s_gw(2,2,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*ky) + exp(-ci*ky) )

                  s_gw(1,2,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*kx) )
                  s_gw(2,1,:,ik,:) = simp(2,1,:,:)*( 0.0_kr + exp(-ci*kx) )
               endif

               ! plaquette: We leave out the local part so s_gw is purely nonlocal
               if (1==0) then
                  s_gw(1,2,:,ik,:) = simp(1,2,:,:)*( 0.0_kr + exp(+ci*kx)  )
                  s_gw(2,1,:,ik,:) = simp(2,1,:,:)*( 0.0_kr + exp(-ci*kx)  )

                  s_gw(1,3,:,ik,:) = simp(1,3,:,:)*( 0.0_kr + exp(+ci*ky)  )
                  s_gw(3,1,:,ik,:) = simp(3,1,:,:)*( 0.0_kr + exp(-ci*ky)  )

                  s_gw(2,4,:,ik,:) = simp(2,4,:,:)*( 0.0_kr + exp(+ci*ky)  )
                  s_gw(4,4,:,ik,:) = simp(4,2,:,:)*( 0.0_kr + exp(-ci*ky)  )

                  s_gw(3,4,:,ik,:) = simp(3,4,:,:)*( 0.0_kr + exp(+ci*kx)  )
                  s_gw(4,3,:,ik,:) = simp(4,3,:,:)*( 0.0_kr + exp(-ci*kx)  )

                  s_gw(1,4,:,ik,:) = simp(1,4,:,:)*( 0.0_kr + exp(+ci*kx) + exp(+ci*ky) + exp(+ci*(kx+ky)) )
                  s_gw(4,1,:,ik,:) = simp(4,1,:,:)*( 0.0_kr + exp(-ci*kx) + exp(-ci*ky) + exp(-ci*(kx+ky)) )

                  s_gw(2,3,:,ik,:) = simp(2,3,:,:)*( 0.0_kr + exp(-ci*kx) + exp(+ci*ky) + exp(+ci*(-kx+ky)) )
                  s_gw(3,2,:,ik,:) = simp(3,2,:,:)*( 0.0_kr + exp(+ci*kx) + exp(-ci*ky) + exp(-ci*(-kx+ky)) )
               endif

         enddo ! iky   
      enddo ! ikx

!      write(*,*) "WARNING: WE ALSO REPLACE S_IMP BY THE LOCAL REPERIODIZED SIGMA !!!"
!      simp = sloc
!      do ik = 1,nkpts
!         s_gw(:,:,:,ik,:) = s_gw(:,:,:,ik,:) - sloc
!      enddo

      deallocate( sloc )

   end subroutine reperiod_sigma


!  ============================================================
!  == Calculate the hybridization function
!  ============================================================
   subroutine get_hybrid(hybrid,gbath,mu_loc)
      complex(kr),intent(inout) :: hybrid(:,:,:,:)     ! norb,norb,nspin,nomega
      complex(kr),intent(inout) :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: mu_loc(:,:,:)       ! norb,norb,nspin

      complex(8),allocatable  :: gbath_tmp(:,:)        ! norb,norb
      real(kr),allocatable    :: unity(:,:)            ! norb,norb
      integer(ki)             :: w,s,i
      real(kr)                :: coeffs(5)
  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gbath_tmp(norb,norb) )
      allocate( unity(norb,norb) )

      hybrid = (0.0_kr,0.0_kr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!write(*,*) "WARNING: I am modifying Gbath in get_hybrid !!! Set tail to 1 !!!"
!!do i=1,norb
!!   do s=1,nspin
!!      coeffs = get_highfreq_coeff(gbath(i,i,s,:),0)
!!      write(*,*) 'Coeffs for m=',i,", s=",s," : ", coeffs
!!
!!      gbath(i,i,s,:) = gbath(i,i,s,:)/coeffs(2)
!!   enddo
!!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*) '!!! Reminder: We replace the matrix inversion in get_hybrid by 1/Gbath !!!'

!  ====================================================
!  == Create identity matrix for use with iw_n ==
!  == Convert the local orbital levels into a matrix =
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

      do w=1,nomega
         do s=1,nspin

            ! use a copy of the local Greensfunction
            gbath_tmp = gbath(:,:,s,w)

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gbath_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gbath_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            ! gbath_tmp is now the inverse of g_bath
            
            ! Calculate hybridization function
            hybrid(:,:,s,w) = ci*wn(w-1)*unity + mu_loc(:,:,s) - gbath_tmp

            ! on real frequencies!
!            hybrid(:,:,s,w) = ( -5.0+w*13.0/nomega+0.01*ci )*unity + mu_loc(:,:,s) - gbath_tmp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! not necessary actually
!            do i=1,norb
!               hybrid(i,i,s,w) = ci*wn(w-1) + mu_loc(i,i,s) - 1.0_kr/gbath(i,i,s,w)
!            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         enddo
      enddo

!      write(*,'(A)') 'WARNING: We set the c0 coeff of the Hybrid function to zero, but calculate mu_loc'
!      write(*,'(A)') 'WARNING: from the formula: HartreeFock - g2 '
!      write(*,'(A)') 'WARNING: If these values differ this usually means that the Tail in Gbath has a problem'
!      write(*,'(A)') 'WARNING: This keeps all other coeffs untouched but improves the Fourier transform'
!      write(*,'(A)') 'WARNING: and is probably the best compromise'
!      write(*,'(A)',advance='no') 'Hybrid coeff: '

!      do s=1,nspin
!         do i=1,norb
!            coeffs = get_highfreq_coeff(hybrid(i,i,s,:),0)
!            write(*,'((F7.3),2X)',advance='no') coeffs(1)
!            do w=1,nomega
!               hybrid(i,i,s,w) = hybrid(i,i,s,w) - coeffs(1)
!            end do
!         end do
!      end do
      write(*,'(A)') ' '

      deallocate(gbath_tmp)
      deallocate(unity)
   end subroutine get_hybrid

!  ============================================================
!  == Calculate the hybridization function with a different method
!  ============================================================
   subroutine get_hybrid2(hybrid,mu_loc,gloc,gimp)
      complex(kr),intent(inout) :: hybrid(:,:,:,:)     ! norb,norb,nspin,nomega
      real(kr),intent(inout)    :: mu_loc(:,:,:)       ! norb,norb,nspin
      complex(kr),intent(in)    :: gloc(:,:,:,:)      ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: gimp(:,:,:,:)      ! norb,norb,nspin,nomega

      complex(8),allocatable  :: gloc_tmp(:,:)        ! norb,norb
      complex(8),allocatable  :: gimp_tmp(:,:)        ! norb,norb
      integer(ki)             :: w,s,i
      real(kr)                :: coeffs(5)
  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gloc_tmp(norb,norb) )
      allocate( gimp_tmp(norb,norb) )

      do w=1,nomega
         do s=1,nspin

            ! use a copy of the local Greensfunction
            gloc_tmp = gloc(:,:,s,w)
            gimp_tmp = gimp(:,:,s,w)

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gloc_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gloc_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gimp_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gimp_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            
            ! Calculate hybridization function
            hybrid(:,:,s,w) = hybrid(:,:,s,w) - gloc_tmp + gimp_tmp
         enddo
      enddo


      do s=1,nspin
         do i=1,norb
            coeffs = get_highfreq_coeff(hybrid(i,i,s,:),0)
            do w=1,nomega
               hybrid(i,i,s,w) = hybrid(i,i,s,w) - coeffs(1)
            end do
            !mu_loc(i,i,s) = -coeffs(1)
         end do
      end do
      write(*,'(A)') ' '

      deallocate(gloc_tmp)
      deallocate(gimp_tmp)
   end subroutine get_hybrid2

!  ============================================================
!  == Calculate the bath from hybrid
!  ============================================================
   subroutine get_weiss_field_from_hybrid(gbath,mu_loc,hybrid)
      complex(kr),intent(inout) :: gbath(:,:,:,:)      ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: mu_loc(:,:,:)       ! norb,norb,nspin
      complex(kr),intent(in)    :: hybrid(:,:,:,:)     ! norb,norb,nspin,nomega

      complex(8),allocatable  :: gbath_tmp(:,:)        ! norb,norb
      real(kr),allocatable    :: unity(:,:)            ! norb,norb
      integer(ki)             :: w,s,i
      real(kr)                :: coeffs(5)
  !   Lapack stuff
      complex(8), dimension(norb) :: work  ! work array for LAPACK
      integer, dimension(norb)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( gbath_tmp(norb,norb) )
      allocate( unity(norb,norb) )

!  ====================================================
!  == Create identity matrix for use with iw_n ==
!  == Convert the local orbital levels into a matrix =
      unity = (0.0_kr)
      do i=1,norb
         unity(i,i) = 1.0_kr
      enddo
!  ====================================================

      do w=1,nomega
         do s=1,nspin

            ! this is gbath inverse
            gbath_tmp = ci*wn(w-1)*unity + mu_loc(:,:,s) - hybrid(:,:,s,w)

            ! Then invert with LAPACK
            call ZGETRF(norb, norb, gbath_tmp, norb, ipiv, info_lpck)
            if (info_lpck /= 0) then
               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
                         & info_lpck
               stop 
            end if
            call ZGETRI(norb, gbath_tmp, norb, ipiv, work, norb, info_lpck)
            if (info_lpck /= 0) then
              stop 'Matrix inversion failed!'
            end if
            
            ! Calculate hybridization function
            gbath(:,:,s,w) = gbath_tmp

         enddo
      enddo

      deallocate(gbath_tmp)
      deallocate(unity)
   end subroutine get_weiss_field_from_hybrid
   
!  ============================================================
!  == Calculate and write hybrid function on imaginary time
!  == BUT: Write only data for first equivalent DMFT atom !
!  ============================================================
   subroutine write_imag_time_hybrid(filename,hybrid)
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: hybrid(:,:,:,:)     ! norbPerAtom,norbPerAtom,nspin,nomega

      complex(kr)             :: delta
      complex(kr),allocatable :: deltam(:,:,:)                       ! deltam(norb,norb,nspin)
      real(kr)                :: tau,c1,c2,c3,c4,wnn,coswnt,sinwnt,coeffs(5)
      integer(ki),parameter   :: iounit = 10, iounit2=11, iounit3=12
      integer(ki)             :: w,s,s2,m1,m2,t, a, i,j, ndim
      character(len=1024)     :: path1,path2,path3

      ! Be careful of the ordering ! 
      ! For the CTHYB spin runs first, then orbitals !

      ndim = size(hybrid(:,1,1,1))

      allocate( deltam(ndim,ndim,nspin) )

      write(*,'(A,A)') 'Writing hybridization(tau) to ',filename

      open(unit=iounit, file=filename//".dat",              status="unknown")
      open(unit=iounit2,file=filename//"_matrix.dat",       status="unknown")
      open(unit=iounit3,file=filename//"_matrix_solver.dat",status="unknown")
      do t=1,ntau+1
         call progress(t,ntau+1)

         tau = (t-1)*beta/ntau
         write(iounit,'((ES23.15E3),(4X))',advance='no') tau
         write(iounit2,'((ES23.15E3),(4X))',advance='no') tau
         do m1=1,ndim
            do m2=1,ndim
               do s=1,nspin
                  delta = 0.0_kr
                  if (m1==m2) then  ! only use HF correction for diagonals
                     coeffs = get_highfreq_coeff(hybrid(m1,m2,s,:),0)
                     c1 = coeffs(2)
                     c2 = coeffs(3)
                     c3 = coeffs(4)
                     c4 = coeffs(5)

                     ! Use high-frequency correction up to 4-th order
                     ! assume positive tau for high freq. correction
                     delta = -0.5_kr*c1  + (2*tau-beta)*c2/4  &
                            & +tau*(beta-tau)*c3/4            &
                            & +(2*tau-beta)*(2*tau**2-2*tau*beta-beta**2)*c4/48
                           
                        ! now sum up the function - correction itself
                        do w=1,nomega
                           wnn = wn(w-1)
                           sinwnt = sin(wnn*tau)
                           coswnt = cos(wnn*tau)
   
                           delta = delta + (real(hybrid(m1,m2,s,w))*coswnt  &
                         &               + aimag(hybrid(m1,m2,s,w))*sinwnt  &
                         &               + c1*sinwnt/wnn                    &
                         &               + c2*coswnt/wnn**2                 &
                         &               - c3*sinwnt/wnn**3                 &
                         &               - c4*coswnt/wnn**4 )*2/beta                               
                        enddo  

                     ! keep it causal when there is Monte Carlo noise
                      if (real(delta)>=0.0_kr) delta = -1.0E-10
                      write(iounit,'(ES23.15E3,4X)',advance='no') real(delta)  ! used for SEGMENT solver

                  else ! i!=j
                     delta = 0.0_kr
                     do w=1,nomega
                        wnn = wn(w-1)
                        delta = delta + (       hybrid(m1,m2,s,w) *exp(-ci*wnn*tau)       &
                                    &    +conjg(hybrid(m2,m1,s,w))*exp(+ci*wnn*tau) )/beta
                     enddo
                  endif ! i==j
   
                  ! general matrix format for plotting
                  write(iounit2,'(2(ES23.15E3,4X))',advance='no') real(delta),aimag(delta)

                  ! save for later
                  deltam(m1,m2,s) = delta
               enddo ! s loop
            enddo ! m2 loop
         enddo ! m1 loop
         write(iounit,'(A1)') ' '
         write(iounit2,'(A1)') ' '

         ! matrix solver format
         do m1=0,ndim-1
            do s=0,nspin-1
               do m2=0,ndim-1
                  do s2=0,nspin-1
                     if (s==s2) then
                        write(iounit3,'((I6,4X,I3,4X,I3,4X,ES17.10,4X,ES17.10,4X))') t-1,(m1*2+s),(m2*2+s2), & 
                                                            &  real(deltam(m1+1,m2+1,s+1)),aimag(deltam(m1+1,m2+1,s+1))
                     else
                        write(iounit3,'((I6,4X,I3,4X,I3,4X,ES17.10,4X,ES17.10,4X))') t-1,(m1*2+s),(m2*2+s2),0.0,0.0
                     endif
                  enddo
               enddo
            enddo
         enddo

      enddo ! t loop
      write(*,'(A)') ' '
      close(iounit)
      close(iounit2)
      close(iounit3)

      deallocate(deltam)

   end subroutine write_imag_time_hybrid

!  ============================================================
!  == Calculate and write bath Green's function on imaginary time
!  ============================================================
   subroutine transform_to_imag_time(filename,func)
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:)     ! ndim,ndim,nspin,nomega

      complex(kr)             :: delta
      real(kr)                :: tau,c1,c2,c3,c4,wnn,coswnt,sinwnt,coeffs(5)
      integer(ki),parameter   :: iounit = 10, iounit2=11, iounit3=12
      integer(ki)             :: w,s,s2,i,m1,m2,t
      complex(kr),allocatable :: deltam(:,:,:,:) ! (ndim,ndim,nspin,ntau+1)
      integer(ki) :: ndim


      ndim = size(func(:,1,1,1))
      allocate( deltam(ndim,ndim,nspin,ntau+1) )

      open(unit=iounit,file=filename//".dat",status="unknown")
      open(unit=iounit2,file=filename//"_matrix.dat",status="unknown")

      do t=1,ntau+1
         call progress(t,ntau+1)
      
         tau = (t-1)*beta/ntau
         write(iounit,'((ES23.15E3),(4X))',advance='no') tau
         write(iounit2,'((ES23.15E3),(4X))',advance='no') tau
         do m1=1,ndim
            do m2=1,ndim
               do s=1,nspin
                  if (m1==m2) then
                     coeffs = get_highfreq_coeff(func(m1,m2,s,:),1)
                     c1 = coeffs(2)
                     c2 = coeffs(3)
                     c3 = coeffs(4)
                     c4 = coeffs(5)
   
                     ! Use high-frequency correction up to 4-th order
                     ! assume positive tau for high freq. correction
                     delta = -0.5_kr*c1  + (2*tau-beta)*c2/4  &
                            & +tau*(beta-tau)*c3/4            &
                            & +(2*tau-beta)*(2*tau**2-2*tau*beta-beta**2)*c4/48
                           
                        ! now sum up the function - correction itself
                        do w=1,nomega
                           wnn = wn(w-1)
                           sinwnt = sin(wnn*tau)
                           coswnt = cos(wnn*tau)
   
                           delta = delta + (real(func(m1,m2,s,w))*coswnt  &
                         &               + aimag(func(m1,m2,s,w))*sinwnt  &
                         &               + c1*sinwnt/wnn                    &
                         &               + c2*coswnt/wnn**2                 &
                         &               - c3*sinwnt/wnn**3                 &
                         &               - c4*coswnt/wnn**4 )*2/beta                               
                        enddo  
   
                     ! keep it causal when there is Monte Carlo noise
                     if (m1==m2) then
                         if (real(delta)>=0.0_kr) delta = -1.0E-10
                         write(iounit,'(ES23.15E3,4X)',advance='no') real(delta)  ! used for SEGMENT solver
                     endif
                  
                  else ! m1!=m2
                     delta = 0.0_kr
                     do w=1,nomega
                        wnn = wn(w-1)
                        delta = delta + (       func(m1,m2,s,w) *exp(-ci*wnn*tau)       &
                                    &    +conjg(func(m2,m1,s,w))*exp(+ci*wnn*tau) )/beta
                     enddo
                  endif

                  ! general matrix format for plotting
                  write(iounit2,'(2(ES23.15E3,4X))',advance='no') real(delta),aimag(delta)

                  ! save for later
                  deltam(m1,m2,s,t) = delta
               enddo ! s loop
            enddo ! m2 loop
         enddo ! m1 loop
         write(iounit,'(A1)') ' '
         write(iounit2,'(A1)') ' '

      enddo ! t loop
      write(*,'(A)') ' '
      close(iounit)
      close(iounit2)

      open(unit=iounit3,file=filename//"_matrix_solver.dat",status="unknown")
      write(iounit3,'((I1),(3X),(I3),(3X),(I7),(3X),(F10.4))') 2, ndim, ntau+1, beta
      ! matrix solver format
      do s=0,nspin-1
         do m1=0,ndim-1
            do m2=0,ndim-1
               do t=1,ntau+1
                  write(iounit3,'((I1,4X,I4,4X,I4,4X,I7,4X,ES17.10,4X,ES17.10,4X))') s, m1, m2, t-1, & 
                                                             &  real(deltam(m1+1,m2+1,s+1,t)),aimag(deltam(m1+1,m2+1,s+1,t))
               enddo ! t
            enddo ! m2
         enddo ! m1
      enddo ! s


      close(iounit3)

      deallocate( deltam )

   end subroutine transform_to_imag_time

!  ============================================================
!  == Calculate and write bath Green's function on imaginary time
!  == This is for CT-INT 
!  ============================================================
   subroutine write_imag_time_bath(filename,hybrid,mu_loc)
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: hybrid(:,:,:,:)     ! ndim,ndim,nspin,nomega
      complex(kr),intent(in)       :: mu_loc(:,:,:)      ! ndim,ndim,nspin

      complex(kr)             :: delta
      real(kr)                :: tau,c1,c2,c3,c4,wnn,coswnt,sinwnt,coeffs(5)
      integer(ki),parameter   :: iounit = 10, iounit2=11, iounit3=12
      integer(ki)             :: w,s,s2,i,j,m1,m2,t,a, ndim
      complex(kr),allocatable :: ginv(:,:)         ! ndim,ndim
      complex(kr),allocatable :: deltam(:,:,:,:) ! (ndim,ndim,nspin,ntau+1)
      complex(kr),allocatable :: gbath(:,:,:,:)  ! (ndim,ndim,nspin,nomega)
      real(kr),allocatable    :: unity(:,:)      ! ndim,ndim
      character(len=1024)     :: path1,path2,path3

      ! first generate gbath with shifted mu_loc

      ndim = size(hybrid(:,1,1,1))

      allocate (gbath(ndim,ndim,nspin,nomega) )
      allocate( deltam(ndim,ndim,nspin,ntau+1) )
      allocate (ginv(ndim,ndim) )
      allocate (unity(ndim,ndim) )

      unity = (0.0_kr)
      do i=1,ndim
         unity(i,i) = 1.0_kr
      enddo

      do w=1,nomega
         do s=1,nspin
            ! We can use Uinput here, CT-INT only works with constant U and single orb
            ginv = unity*( ci*wn(w-1) - 0.5*Uinput ) + mu_loc(:,:,s) - hybrid(:,:,s,w)
            
            call get_inverse( gbath(:,:,s,w), ginv, ndim )  
         enddo ! s
      enddo ! w

      ! Be careful of the ordering ! 
      ! For the CTHYB spin runs first, then orbitals !

      deltam = (0.0_kr,0.0_kr)

      open(unit=iounit, file=filename//".dat",              status="unknown")
      open(unit=iounit2,file=filename//"_matrix.dat",       status="unknown")
      open(unit=iounit3,file=filename//"_matrix_solver.dat",status="unknown")

      do t=1,ntau+1
         call progress(t,ntau+1)
   
         tau = (t-1)*beta/ntau
         write(iounit,'((ES23.15E3),(4X))',advance='no') tau
         write(iounit2,'((ES23.15E3),(4X))',advance='no') tau
         do m1=1,ndim
            do m2=1,ndim
               do s=1,nspin
                  if (m1==m2) then
                     coeffs = get_highfreq_coeff(gbath(m1,m2,s,:),1)
                     c1 = coeffs(2)
                     c2 = coeffs(3)
                     c3 = coeffs(4)
                     c4 = coeffs(5)

                     ! Use high-frequency correction up to 4-th order
                     ! assume positive tau for high freq. correction
                     delta = -0.5_kr*c1  + (2*tau-beta)*c2/4  &
                            & +tau*(beta-tau)*c3/4            &
                            & +(2*tau-beta)*(2*tau**2-2*tau*beta-beta**2)*c4/48
                           
                        ! now sum up the function - correction itself
                        do w=1,nomega
                           wnn = wn(w-1)
                           sinwnt = sin(wnn*tau)
                           coswnt = cos(wnn*tau)
   
                           delta = delta + (real(gbath(m1,m2,s,w))*coswnt  &
                         &               + aimag(gbath(m1,m2,s,w))*sinwnt  &
                         &               + c1*sinwnt/wnn                    &
                         &               + c2*coswnt/wnn**2                 &
                         &               - c3*sinwnt/wnn**3                 &
                         &               - c4*coswnt/wnn**4 )*2/beta                               
                        enddo ! w
   
                     ! keep it causal when there is Monte Carlo noise
                     if (real(delta)>=0.0_kr) delta = -1.0E-10
                     write(iounit,'(ES23.15E3,4X)',advance='no') real(delta)  ! used for SEGMENT solver
                  
                  else ! i/=j
                     delta = 0.0_kr
                     do w=1,nomega
                        wnn = wn(w-1)
                        delta = delta + (       gbath(m1,m2,s,w) *exp(-ci*wnn*tau)       &
                                    &    +conjg(gbath(m2,m1,s,w))*exp(+ci*wnn*tau) )/beta
                     enddo !w
                  endif !i==j?

                  ! general matrix format for plotting
                  write(iounit2,'(2(ES23.15E3,4X))',advance='no') real(delta),aimag(delta)

                  ! save for later
                  deltam(m1,m2,s,t) = delta
               enddo ! s loop
            enddo ! m2 loop
         enddo ! m1 loop
         write(iounit,'(A1)') ' '
         write(iounit2,'(A1)') ' '
      enddo ! t loop
      write(*,'(A)') ' '
      close(iounit)
      close(iounit2)
   
      write(iounit3,'((I1),(3X),(I3),(3X),(I7),(3X),(F10.4))') 2, ndim, ntau+1, beta
      ! matrix solver format
      do s=0,nspin-1
         do m1=0,ndim-1
            do m2=0,ndim-1
               do t=1,ntau+1
                  write(iounit3,'((I1,4X,I4,4X,I4,4X,I7,4X,ES17.10,4X,ES17.10,4X))') s, m1, m2, t-1, & 
                                                             &  real(deltam(m1+1,m2+1,s+1,t)),aimag(deltam(m1+1,m2+1,s+1,t))
               enddo ! t
            enddo ! m2
         enddo ! m1
      enddo ! s
      close(iounit3)

      deallocate( deltam )
      deallocate( gbath )
      deallocate (ginv )
      deallocate (unity )

   end subroutine write_imag_time_bath

!  ============================================================
!  == Calculate the high-freq. expansion coeffs. with least squares
!  ============================================================
   function get_highfreq_coeff(f,gf) result(coeffs)
        complex(kr),intent(in) :: f(:)    ! nomega
        integer,intent(in)     :: gf      ! 0=do general fit, 1=assume that c0=0 and c1=1
        real(kr)               :: coeffs(5)
        integer  :: i,n,imagOrder,realOrder

        ! matrices for least square fit
        real(8),allocatable :: XmatR(:,:), XmatinvR(:,:), betaMatR(:), YmatR(:)
        real(8),allocatable :: XmatI(:,:), XmatinvI(:,:), betaMatI(:), YmatI(:)
        real(8),allocatable :: ls_workR(:),ls_workI(:)         ! work array for LAPACK
        integer,allocatable :: ls_ipivR(:),ls_ipivI(:)         ! pivot indices
        integer :: info_lpck
        external DGETRF
        external DGETRI

!if (gf == 666) then
        realOrder = 2
        imagOrder = 2
!else
!        realOrder = 3
!        imagOrder = 2
!endif
!        if (gf==1) then
!            realOrder = 2
!            imagOrder = 1
!        end if
        allocate( XmatR   (nfit     ,realOrder) )
        allocate( XmatinvR(realOrder,realOrder) )
        allocate( betaMatR(realOrder) )
        allocate( YmatR   (nfit) )
        allocate( ls_workR(realOrder) )
        allocate( ls_ipivR(realOrder) )

        allocate( XmatI   (nfit     ,imagOrder) )
        allocate( XmatinvI(imagOrder,imagOrder) )
        allocate( betaMatI(imagOrder) )
        allocate( YmatI   (nfit) )
        allocate( ls_workI(imagOrder) )
        allocate( ls_ipivI(imagOrder) )

        ! Use least square fit
        do i=1,nfit ! frequencies
            YmatR(i) =  real(f(nomega-nfit+i))
            YmatI(i) = aimag(f(nomega-nfit+i))
            do n=0,realOrder-1 ! order of the monome
                XmatR(i,n+1) = real( (ci*wn(nomega-nfit+i-1))**(-2*n) )
            enddo
            do n=0,imagOrder-1 ! order of the monome
                XmatI(i,n+1) = aimag( (ci*wn(nomega-nfit+i-1))**(-2*n-1) )
            enddo
        enddo

        XmatinvR = matmul( transpose(XmatR),XmatR )
        call DGETRF(realOrder, realOrder, XmatinvR, realOrder, ls_ipivR, info_lpck)
        call DGETRI(realOrder, XmatinvR, realOrder, ls_ipivR, ls_workR, realOrder, info_lpck)
        betaMatR = matmul( XmatinvR, matmul(transpose(XmatR), YmatR) )

        XmatinvI = matmul( transpose(XmatI),XmatI )
        call DGETRF(imagOrder, imagOrder, XmatinvI, imagOrder, ls_ipivI, info_lpck)
        call DGETRI(imagOrder, XmatinvI, imagOrder, ls_ipivI, ls_workI, imagOrder, info_lpck)
        betaMatI = matmul( XmatinvI, matmul(transpose(XmatI), YmatI) )

        coeffs = (0.0_kr)
        do n=0,realOrder-1
            coeffs(2*n+1) = betaMatR(n+1)
        enddo
        do n=0,imagOrder-1
            coeffs(2*n+2) = betaMatI(n+1)
        enddo

        deallocate( XmatR )
        deallocate( XmatinvR )
        deallocate( betaMatR )
        deallocate( YmatR   )
        deallocate( ls_workR )
        deallocate( ls_ipivR )

        deallocate( XmatI )
        deallocate( XmatinvI )
        deallocate( betaMatI )
        deallocate( YmatI )
        deallocate( ls_workI)
        deallocate( ls_ipivI )

   end function get_highfreq_coeff

!  ============================================================
!  == Get the linear component in the imaginary part of a GF with least squares
!  ============================================================
   function get_linear_imag_coeff(f) result(coeff)
        complex(kr),intent(in) :: f(:)    ! nomega
        real(kr)               :: coeff
        integer  :: i,n,imagOrder

        ! matrices for least square fit
        real(8),allocatable :: XmatI(:,:), XmatinvI(:,:), betaMatI(:), YmatI(:)
        real(8),allocatable :: ls_workI(:)         ! work array for LAPACK
        integer,allocatable :: ls_ipivI(:)         ! pivot indices
        integer :: info_lpck
        external DGETRF
        external DGETRI

        imagOrder = 2

        allocate( XmatI   (nfit     ,imagOrder+1) )
        allocate( XmatinvI(imagOrder+1,imagOrder+1) )
        allocate( betaMatI(imagOrder+1) )
        allocate( YmatI   (nfit) )
        allocate( ls_workI(imagOrder+1) )
        allocate( ls_ipivI(imagOrder+1) )

        ! Use least square fit
        do i=1,nfit ! frequencies
            YmatI(i) = aimag(f(nomega-nfit+i))

            ! The linear part
            XmatI(i,1) = aimag( (ci*wn(nomega-nfit+i-1))**(1) )

            do n=0,imagOrder-1 ! order of the monome
                XmatI(i,n+2) = aimag( (ci*wn(nomega-nfit+i-1))**(-2*n-1) )
            enddo            

        enddo

        XmatinvI = matmul( transpose(XmatI),XmatI )
        call DGETRF(imagOrder+1, imagOrder+1, XmatinvI, imagOrder+1, ls_ipivI, info_lpck)
        call DGETRI(imagOrder+1, XmatinvI, imagOrder+1, ls_ipivI, ls_workI, imagOrder+1, info_lpck)
        betaMatI = matmul( XmatinvI, matmul(transpose(XmatI), YmatI) )

        coeff = betaMatI(1)

        deallocate( XmatI )
        deallocate( XmatinvI )
        deallocate( betaMatI )
        deallocate( YmatI )
        deallocate( ls_workI)
        deallocate( ls_ipivI )

   end function get_linear_imag_coeff

!  ============================================================
!  == Calculate the low-freq. expansion coeffs. with least squares
!  == This returns a complex array of coefficients:
!  == The real part are the coeffs for the real part, and imag. vice versa
!  ============================================================
   function get_lowfreq_coeff(f,ncoeff,nfit) result(coeffs)
        complex(kr),intent(in) :: f(:)    ! nomega
        integer,intent(in)     :: ncoeff,nfit      ! how many coefficients, how many freqs for fitting
        complex(kr)               :: coeffs(ncoeff)
        integer  :: i,n

        ! matrices for least square fit
        real(8),allocatable :: Xmat(:,:), Xmatinv(:,:)
        real(8),allocatable :: betaMatI(:), YmatI(:), betaMatR(:), YmatR(:)
        real(8),allocatable :: ls_workR(:),ls_workI(:)         ! work array for LAPACK
        integer,allocatable :: ls_ipivR(:),ls_ipivI(:)         ! pivot indices
        integer :: info_lpck
        external DGETRF
        external DGETRI

        allocate( Xmat   (nfit     ,ncoeff) )
        allocate( Xmatinv(ncoeff,ncoeff) )
        allocate( betaMatR(ncoeff) )
        allocate( YmatR   (nfit) )
        allocate( ls_workR(ncoeff) )
        allocate( ls_ipivR(ncoeff) )

        allocate( betaMatI(ncoeff) )
        allocate( YmatI   (nfit) )
        allocate( ls_workI(ncoeff) )
        allocate( ls_ipivI(ncoeff) )

        ! Use least square fit
        do i=1,nfit ! frequencies
            YmatR(i) =  real(f(i))
            YmatI(i) = aimag(f(i))
            do n=0,ncoeff-1 ! order of the monome
                Xmat(i,n+1) = wn(i-1)**n
            enddo
        enddo

        Xmatinv = matmul( transpose(Xmat),Xmat )
        call DGETRF(ncoeff, ncoeff, Xmatinv, ncoeff, ls_ipivR, info_lpck)
        call DGETRI(ncoeff, Xmatinv, ncoeff, ls_ipivR, ls_workR, ncoeff, info_lpck)
        betaMatR = matmul( Xmatinv, matmul(transpose(Xmat), YmatR) )
        betaMatI = matmul( Xmatinv, matmul(transpose(Xmat), YmatI) )

        coeffs = (0.0_kr,0.0_kr)
        do n=1,ncoeff
            coeffs(n) = betaMatR(n) + ci*betaMatI(n)
        enddo

        deallocate( Xmat )
        deallocate( Xmatinv )
        deallocate( betaMatR )
        deallocate( YmatR   )
        deallocate( ls_workR )
        deallocate( ls_ipivR )

        deallocate( betaMatI )
        deallocate( YmatI )
        deallocate( ls_workI)
        deallocate( ls_ipivI )

   end function get_lowfreq_coeff

!  ============================================================
!  == Calculate the low-freq. expansion coeffs. with least squares FOR BOSONIC FREQUENCIES
! WE ASSUME n=0 is not given !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FIRST ELEMENT IS N=1
!  == This returns a complex array of coefficients:
!  == The real part are the coeffs for the real part, and imag. vice versa
!  ============================================================
   function get_lowfreq_coeff_bos(f,ncoeff,nfit) result(coeffs)
        complex(kr),intent(in) :: f(:)    ! nomega
        integer,intent(in)     :: ncoeff,nfit      ! how many coefficients, how many freqs for fitting
        complex(kr)               :: coeffs(ncoeff)
        integer  :: i,n

        ! matrices for least square fit
        real(8),allocatable :: Xmat(:,:), Xmatinv(:,:)
        real(8),allocatable :: betaMatI(:), YmatI(:), betaMatR(:), YmatR(:)
        real(8),allocatable :: ls_workR(:),ls_workI(:)         ! work array for LAPACK
        integer,allocatable :: ls_ipivR(:),ls_ipivI(:)         ! pivot indices
        integer :: info_lpck
        external DGETRF
        external DGETRI

        allocate( Xmat   (nfit     ,ncoeff) )
        allocate( Xmatinv(ncoeff,ncoeff) )
        allocate( betaMatR(ncoeff) )
        allocate( YmatR   (nfit) )
        allocate( ls_workR(ncoeff) )
        allocate( ls_ipivR(ncoeff) )

        allocate( betaMatI(ncoeff) )
        allocate( YmatI   (nfit) )
        allocate( ls_workI(ncoeff) )
        allocate( ls_ipivI(ncoeff) )

        ! Use least square fit
        do i=1,nfit ! frequencies
            YmatR(i) =  real(f(i))
            YmatI(i) = aimag(f(i))
            do n=0,ncoeff-1 ! order of the monome
                Xmat(i,n+1) = vn(i-0)**n
            enddo
        enddo

        Xmatinv = matmul( transpose(Xmat),Xmat )
        call DGETRF(ncoeff, ncoeff, Xmatinv, ncoeff, ls_ipivR, info_lpck)
        call DGETRI(ncoeff, Xmatinv, ncoeff, ls_ipivR, ls_workR, ncoeff, info_lpck)
        betaMatR = matmul( Xmatinv, matmul(transpose(Xmat), YmatR) )
        betaMatI = matmul( Xmatinv, matmul(transpose(Xmat), YmatI) )

        coeffs = (0.0_kr,0.0_kr)
        do n=1,ncoeff
            coeffs(n) = betaMatR(n) + ci*betaMatI(n)
        enddo

        deallocate( Xmat )
        deallocate( Xmatinv )
        deallocate( betaMatR )
        deallocate( YmatR   )
        deallocate( ls_workR )
        deallocate( ls_ipivR )

        deallocate( betaMatI )
        deallocate( YmatI )
        deallocate( ls_workI)
        deallocate( ls_ipivI )

   end function get_lowfreq_coeff_bos
   
!  ============================================================
!  == Calculate the high-freq. expansion coeffs. for general Matsubara function
!  ============================================================
   function get_highfreq_coeff_fit(f,gf) result(coeffs)
      complex(kr),intent(in) :: f(:)    ! nomega
      integer,intent(in)     :: gf      ! 0=do general fit, 1=assume that c0=0 and c1=1 
      real(kr)               :: coeffs(5)
      integer  :: i,navg=20,dist=30      ! average the coeffs, take last points seperated by dist
      real(kr) :: wbeg,wmid,wend,frbeg,frmid,frend,fimid,fiend,denom,denom1,denom23


      if ( SIZE(coeffs) /= 5) stop 'ERROR: High freq. coeff. of order >4 not supported'
      
      coeffs = (0.0_kr)
      do i=1,navg
         wbeg = wn(nomega-1-2*dist-i)
         wmid = wn(nomega-1-1*dist-i)
         wend = wn(nomega-1-0*dist-i)
         
         frbeg = real(f(nomega-2*dist-i))
         frmid = real(f(nomega-1*dist-i))
         frend = real(f(nomega-0*dist-i))
         
         fimid = aimag(f(nomega-1*dist-i))
         fiend = aimag(f(nomega-0*dist-i))
         
         denom = wend**2 - wmid**2;
         
         if (gf==0) then ! general fit, adjust all coeffs
            denom1 = ( wend*wend-wmid*wmid )                       &
               &    *( wend*wend*wmid*wmid - wend*wend*wbeg*wbeg   &
               &     + wbeg*wbeg*wbeg*wbeg - wmid*wmid*wbeg*wbeg )
               
            denom23 = (wbeg**4)*( wend**2 - wmid**2 )              &
               &    + (wmid**4)*( wbeg**2 - wend**2 )              &
               &    + (wend**4)*( wmid**2 - wbeg**2 )              
               
            coeffs(1) = coeffs(1) + ( (wend**4)*frend*(wmid**2-wbeg**2)      &
                      &             + (wmid**4)*frmid*(wbeg**2-wend**2)      &
                      &             + (wbeg**4)*frbeg*(wend**2-wmid**2) )/denom1    
                      
            coeffs(2) = coeffs(2) + ( fimid*wmid**3 - fiend*wend**3 )/denom    
            
            coeffs(3) = coeffs(3) + ( ((wbeg*wend)**4)*(frbeg-frend)         &
                      &             + ((wbeg*wmid)**4)*(frmid-frbeg)         &
                      &             + ((wmid*wend)**4)*(frend-frmid) )/denom23
                      
            coeffs(4) = coeffs(4) + ((wmid*wend)**2)*( fimid*wmid - fiend*wend )/denom
            
            coeffs(5) = coeffs(5) + ((wbeg*wmid*wend)**2)                   &
                      &      *( ((wbeg*wend)**2)*( frbeg-frend )            &
                      &       + ((wbeg*wmid)**2)*( frmid-frbeg )            &
                      &       + ((wmid*wend)**2)*( frend-frmid ) )/denom23
                      
         else ! assume we have a Greensfunction with c0=0 and c1=1, then things are easy
            coeffs(1) = coeffs(1) + 0.0_kr;
            coeffs(2) = coeffs(2) + 1.0_kr;
            coeffs(3) = coeffs(3) + ( (wmid**4)*frmid - (wend**4)*frend )/denom;
            coeffs(4) = coeffs(4) + (fiend*wend + 1.0_kr)*wend*wend;
            coeffs(5) = coeffs(5) + ((wend*wmid)**2)*( frmid*wmid**2 - frend*wend**2 ) &
                                  &                   /denom;
         end if
      enddo    
      coeffs = coeffs/navg  
   end function get_highfreq_coeff_fit

!  ============================================================
!  == Calculate the high-freq. expansion coeffs. for general bosonic Matsubara function
! == Will return c0 and c2 only !
!  ============================================================
   function get_highfreq_coeff_bosonic(f) result(coeffs)
        complex(kr),intent(in) :: f(:)    ! nnu
!        integer,intent(in)     :: tozero      ! 0=do general fit, 1=assume that c0=0
        real(kr)               :: coeffs(2)
        integer  :: i,n,imagOrder,realOrder

        ! matrices for least square fit
        real(8),allocatable :: XmatR(:,:), XmatinvR(:,:), betaMatR(:), YmatR(:)
        real(8),allocatable :: XmatI(:,:), XmatinvI(:,:), betaMatI(:), YmatI(:)
        real(8),allocatable :: ls_workR(:),ls_workI(:)         ! work array for LAPACK
        integer,allocatable :: ls_ipivR(:),ls_ipivI(:)         ! pivot indices
        integer :: info_lpck
        external DGETRF
        external DGETRI

        realOrder = 2
        imagOrder = 1 ! just use one, but don't return it

!        if (gf==1) then
!            realOrder = 2
!            imagOrder = 1
!        end if
        allocate( XmatR   (nfit     ,realOrder) )
        allocate( XmatinvR(realOrder,realOrder) )
        allocate( betaMatR(realOrder) )
        allocate( YmatR   (nfit) )
        allocate( ls_workR(realOrder) )
        allocate( ls_ipivR(realOrder) )

        allocate( XmatI   (nfit     ,imagOrder) )
        allocate( XmatinvI(imagOrder,imagOrder) )
        allocate( betaMatI(imagOrder) )
        allocate( YmatI   (nfit) )
        allocate( ls_workI(imagOrder) )
        allocate( ls_ipivI(imagOrder) )

        ! Use least square fit
        do i=1,nfit ! frequencies
            YmatR(i) =  real(f(nomega-nfit+i))
            YmatI(i) = aimag(f(nomega-nfit+i))
            do n=0,realOrder-1 ! order of the monome
                XmatR(i,n+1) = real( (ci*vn(nomega-nfit+i-1))**(-2*n) )
            enddo
            do n=0,imagOrder-1 ! order of the monome
                XmatI(i,n+1) = aimag( (ci*vn(nomega-nfit+i-1))**(-2*n-1) )
            enddo
        enddo

        XmatinvR = matmul( transpose(XmatR),XmatR )
        call DGETRF(realOrder, realOrder, XmatinvR, realOrder, ls_ipivR, info_lpck)
        call DGETRI(realOrder, XmatinvR, realOrder, ls_ipivR, ls_workR, realOrder, info_lpck)
        betaMatR = matmul( XmatinvR, matmul(transpose(XmatR), YmatR) )

        XmatinvI = matmul( transpose(XmatI),XmatI )
        call DGETRF(imagOrder, imagOrder, XmatinvI, imagOrder, ls_ipivI, info_lpck)
        call DGETRI(imagOrder, XmatinvI, imagOrder, ls_ipivI, ls_workI, imagOrder, info_lpck)
        betaMatI = matmul( XmatinvI, matmul(transpose(XmatI), YmatI) )

        coeffs = (0.0_kr)
!        do n=0,realOrder-1
!            coeffs(2*n+1) = betaMatR(n+1)
!        enddo
!        do n=0,imagOrder-1
!            coeffs(2*n+2) = betaMatI(n+1)
!        enddo
         coeffs(1) = betaMatR(1)
         coeffs(2) = betaMatR(2)

        deallocate( XmatR )
        deallocate( XmatinvR )
        deallocate( betaMatR )
        deallocate( YmatR   )
        deallocate( ls_workR )
        deallocate( ls_ipivR )

        deallocate( XmatI )
        deallocate( XmatinvI )
        deallocate( betaMatI )
        deallocate( YmatI )
        deallocate( ls_workI)
        deallocate( ls_ipivI )
   end function get_highfreq_coeff_bosonic

   
!  ============================================================
!  == Obtain the charge of a local Green's function
!  ============================================================
   subroutine get_local_charge(charge,gloc)
      real(kr),intent(inout)    :: charge(:,:)        ! norb,nspin
      complex(kr),intent(in)    :: gloc(:,:,:,:)      ! norb,norb,nspin,nomega
      integer(ki)               :: s,m
      
      do s=1,nspin
         do m=1,norb
            charge(m,s) = filling( gloc(m,m,s,:) )
         enddo
      enddo
   end subroutine get_local_charge
  
!  ============================================================
!  == Evaluate the GW selfenergy at w=0 and make it static
!  ============================================================
   subroutine make_nonloc_qp(s_gw)
      complex(kr),intent(inout) :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega

      integer     :: k,n,m1,m2
      complex(kr),allocatable :: sig_loc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),allocatable :: sig_loc0(:,:,:)     ! norb,norb,nspin
      complex(kr),allocatable :: sig_nonloc0(:,:,:)  ! norb,norb,nspin
      allocate( sig_loc(norb,norb,nspin,nomega) )
      allocate( sig_nonloc0(norb,norb,nspin) )
      allocate( sig_loc0(norb,norb,nspin) )

      if ( useSigK/=0 ) then

         write(*,*) '!!! WARNING !!!!!!!!!!!!!!!!!!!!!'
         write(*,*) '!!! QP-lize the Selfenergy !!!!!!'
         write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'     

         ! first get local Selfenergy 
         sig_loc = (0.0_kr, 0.0_kr)
         do k=1,nkpts
            sig_loc = sig_loc + s_gw(:,:,:,k,:)/nkpts
         enddo

         ! Linearly interpolate the local Selfenergy to zero
         sig_loc0 = sig_loc(:,:,:,1)                              &
                &   - wn(1)*( sig_loc(:,:,:,2)-sig_loc(:,:,:,1) ) &
                &                                /( wn(2)-wn(1) )

         ! then modify k-dependent Selfenergy
         do k=1,nkpts
            ! interpolate the k-dep. Sigma to zero
            sig_nonloc0 = s_gw(:,:,:,k,1)            &
                    &   - wn(1)*( s_gw(:,:,:,k,2)-s_gw(:,:,:,k,1) ) &
                    &                              /( wn(2)-wn(1) )

            do n=1,nomega
               ! Set Sigma to the usual local part and 
               ! add the real shift from w=0
               s_gw(:,:,:,k,n) = sig_loc(:,:,:,n)             &
                           & + real( sig_nonloc0 - sig_loc0 )

               ! Also eavluate the offdiagonal terms at w=0
               ! Only the local diagonal part is thus freq dep.
               do m1=1,norb
                  do m2=1,norb
                     if (m1 /= m2) s_gw(m1,m2,:,k,n) = real( sig_nonloc0(m1,m2,:) )
                  enddo
               enddo
            enddo
         enddo
      endif

      deallocate( sig_loc )
      deallocate( sig_loc0 )
      deallocate( sig_nonloc0 )

   end subroutine make_nonloc_qp

!  ============================================================
!  == Initialize the impurity Selfenergy with insulating solution
!  ============================================================
   subroutine init_insulating(simp)
      complex(kr),intent(inout) :: simp(:,:,:,:)        ! norb,norb,nspin,nomega
      integer(ki)               :: w,m,m1,s,a
      do w=1,nomega
         do a=1,noAtoms
            do s=1,nspin
               do m=1,norbPerAtom(a)
                  m1 = norb_dmft(a,m)
                  simp(m1,m1,s,w) = real(simp(m1,m1,s,w)) + 1.0_kr/( ci*wn(w-1) )
               enddo
            enddo
         enddo
      enddo
   end subroutine init_insulating

 
!  ============================================================
!  == Calculate the Hartree and Fock term for a given Greensfunction
!  == But only for the first equivalent DMFT atom, them copy to other
!  ============================================================
   subroutine get_hartree_fock(hartree,exchange,gloc,uloc)
      real(kr),intent(inout)  :: hartree(:,:,:)       ! norb,norb,nspin
      real(kr),intent(inout)  :: exchange(:,:,:)      ! norb,norb,nspin
      complex(kr),intent(in)  :: gloc(:,:,:,:)        ! norb,norb,nspin,nomega
      complex(kr),intent(in)  :: uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu

      integer(ki)             :: s1,s2,m1,m2, m1d,m2d, a, wend
      real(kr),allocatable    :: charge(:,:)           ! norb,nspin
      real(kr)                :: totCharge

      allocate( charge(norb,nspin) )
      charge   = (0.0_kr)
      hartree  = (0.0_kr)
      exchange = (0.0_kr)
      totCharge = 0.0_kr
      wend = size(uloc(1,1,1,:))

      call get_local_charge(charge,gloc)

      write(*,'(A)') ' Calculate the Hartree and Fock terms for each DMFT atom...'

      do a=1,noAtoms
         ! This is Hartree
         do m1=1,norbPerAtom(a)
            do s1=1,nspin
              m1d = norb_dmft(a,m1)
   
              ! Sum up total Charge for later, this is the charge on the DMFT atom
              totCharge = totCharge + charge( m1d, s1)

              ! This is U_0*( n(m,up) + n(m,down) )
               hartree(m1d,m1d,s1) = hartree(m1d,m1d,s1)            &
                            &       + uloc((m1d-1)*norb+m1d,(m1d-1)*norb+m1d,1,wend)*( charge(m1d,1) + charge(m1d,2) )
   
               do m2=1,norbPerAtom(a)
                  m2d = norb_dmft(a,m2)
                  if (m2d /= m1d) then
                     ! this is (U0-2J_lm)*( n(m2,up) + n(m2,down) )
   
                     hartree(m1d,m1d,s1) = hartree(m1d,m1d,s1)             &
                              &      + uloc((m1d-1)*norb+m1d,(m2d-1)*norb+m2d,1,wend)*( charge(m2d,1) + charge(m2d,2) )
                  endif
               enddo ! m2
           enddo ! s1
        enddo ! m1
     enddo ! a

     ! If we use a frequency dependent interaction the 
     ! Hartree part is modified
     ! We need to calculate Kp(0+) here = ( Uinput-U(0) )/2
!     if (usefreqU==1) then
!
!        do m1=1,norbPerAtom
!           m1d = dmftOrbsIndex(m1)
!           do s1=1,nspin
!              hartree(m1d,m1d,s1) = hartree(m1d,m1d,s1) - 2*(totCharge-0)*(Uinput-Uatzero)/2
!           enddo
!        enddo
!     endif

     ! This is Fock exchange
      do a=1,noAtoms
         do m1=1,norbPerAtom(a)
            do s1=1,nspin
              m1d = norb_dmft(a,m1)

             ! This is -U_0 * n(m1,samespin)
              exchange(m1d,m1d,s1) = exchange(m1d,m1d,s1)            &
                           &       - uloc((m1d-1)*norb+m1d,(m1d-1)*norb+m1d,1,wend)*charge(m1d,s1)

               do m2=1,norbPerAtom(a)
                  m2d = norb_dmft(a,m2)
                 if (m2 /= m1) then
                    ! this is -J * n(m2,samespin)

                    exchange(m1d,m1d,s1) = exchange(m1d,m1d,s1)             &
                   &      - uloc((m1d-1)*norb+m2d,(m2d-1)*norb+m1d,1,wend) * charge(m2d,s1)
                 endif
              enddo ! m2
          enddo ! s1
        enddo ! m1
     enddo ! a Atom loop

     deallocate( charge )


   end subroutine get_hartree_fock

!  ============================================================
!  == Calculate the local orbital levels from
!  == Hyb = iw + mu_loc - Gimp^-1 - Sigma^-1
!  ============================================================
!   subroutine get_loc_levels(mu_loc,gimp,uloc)
!      real(kr),intent(inout)  :: mu_loc(:,:,:)        ! norb,norb,nspin
!      complex(kr),intent(in)  :: gimp(:,:,:,:)        ! norb,norb,nspin,nomega
!      complex(kr),intent(in)     :: uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu
!
!      integer(ki)    :: s,m1,m2,m,w
!      real(kr)       :: coeffs(5)
!      real(kr)       :: hartree(norb,norb,nspin)
!      real(kr)       :: exchange(norb,norb,nspin)
!      complex(8)     :: g_tmp(norb,norb)
!      complex(kr),allocatable :: gimpinv(:,:,:,:)        ! norb,norb,nspin,nomega
!
!  !   Lapack stuff
!      complex(8), dimension(norb) :: work  ! work array for LAPACK
!      integer, dimension(norb)    :: ipiv  ! pivot indices
!      integer :: info_lpck
!      external ZGETRF
!      external ZGETRI
!
!      allocate( gimpinv( norb,norb,nspin,nomega ) )
!
!      ! First invert Gimp for all spins and frequencies
!      do w=1,nomega
!         do s=1,nspin
!            g_tmp = gimp(:,:,s,w)
!
!            ! Then invert with LAPACK
!            call ZGETRF(norb, norb, g_tmp, norb, ipiv, info_lpck)
!            if (info_lpck /= 0) then
!               write(*,*)' ERROR: bath Greensfunction is numerically singular! Return value',&
!                         & info_lpck
!               stop 
!            end if
!            call ZGETRI(norb, g_tmp, norb, ipiv, work, norb, info_lpck)
!            if (info_lpck /= 0) then
!              stop 'Matrix inversion failed!'
!            end if
!
!            gimpinv(:,:,s,w) = g_tmp
!
!         enddo ! s
!      enddo ! w
!
!      ! get the hartree fock terms
!      call get_hartree_fock(hartree,exchange,gimp,uloc)
!
!      write(*,'(A)') ' Hartree-Fock calculated from Gloc for local mu-levels:'
!      do m1=1,norb
!         do s=1,nspin
!            write(*,'(A,I2,A,I2,A,F8.5)') 'orb=',m1,', s=',s,':',hartree(m1,m1,s) + exchange(m1,m1,s)
!         enddo
!      enddo
!
!      write(*,'(A)') ' c0 coeff. of Gimp^-1:'
!      do m1=1,norb
!         do m2=1,norb
!            do s=1,nspin
!               coeffs = get_highfreq_coeff(gimpinv(m1,m2,s,:),0)
!               if (m1==m2) then
!                  write(*,'(A,I2,A,I2,A,F8.5)') 'orb=',m1,', s=',s,': ',coeffs(1)
!               endif
!      
!               mu_loc(m1,m2,s) = coeffs(1) + hartree(m1,m2,s) + exchange(m1,m2,s)
!
!            enddo
!         enddo
!      enddo
!
!      deallocate( gimpinv )
!
!   end subroutine get_loc_levels

!  ============================================================
!  == Calculate the Hartree and Fock Doublecounting terms
!  == Here we only decide whether we keep the DFT HF or
!  == average it over all orbitals (closer to DFT)
!  ============================================================
   subroutine get_hartree_fock_dc(hartree_dc,exchange_dc,hartree,exchange)
      real(kr),intent(inout)  :: hartree_dc(:,:,:)       ! norb,norb,nspin
      real(kr),intent(inout)  :: exchange_dc(:,:,:)      ! norb,norb,nspin
      real(kr),intent(in)     :: hartree(:,:,:)          ! norb,norb,nspin
      real(kr),intent(in)     :: exchange(:,:,:)         ! norb,norb,nspin

      integer(ki)  :: s,m,a,i
      real(kr)     :: tmp

      hartree_dc  = (0.0_kr)
      exchange_dc = (0.0_kr)

      if ( avgHartreeFock == 1 ) then
         write(*,'(A)') ' Average Hartree and Exchange Terms for DFT doublecounting (each spin and atom seperately)...'
         do a=1,noAtoms
   
            do s=1,nspin
               ! Hartree
               tmp = 0.0_kr
               do m=1,norbPerAtom(a)
                  i = norb_dmft(a,m)
                  tmp = tmp + hartree(i,i,s)
               enddo
               do m=1,norbPerAtom(a)
                  i = norb_dmft(a,m)
                  hartree_dc(i,i,s) = tmp/norbPerAtom(a)
               enddo
   
               ! Exchange
               tmp = 0.0_kr
               do m=1,norbPerAtom(a)
                  i = norb_dmft(a,m)
                  tmp = tmp + exchange(i,i,s)
               enddo
               do m=1,norbPerAtom(a)
                  i = norb_dmft(a,m)
                  exchange_dc(i,i,s) = tmp/norbPerAtom(a)
               enddo
            enddo
         enddo ! a atom loop
      else
         hartree_dc  = hartree
         exchange_dc = exchange
      endif
   
   end subroutine get_hartree_fock_dc

   
!  ============================================================
!  == Adjust the real part of the Selfenergy to the full impurity Hartree+Fock term
!  == But only for the DMFT orbitals
!  ============================================================
   subroutine adjust_to_hartree_fock(simp,hartree_dc,exchange,exchange_dc)
      complex(kr),intent(inout) :: simp(:,:,:,:)        ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: hartree_dc(:,:,:)    ! norb,norb,nspin
      real(kr),intent(in)       :: exchange(:,:,:)      ! norb,norb,nspin
      real(kr),intent(in)       :: exchange_dc(:,:,:)   ! norb,norb,nspin

      integer(ki)             :: s1,m1,s2,m2,w,a
      real(kr)                :: coeffs(5)
      real(kr),allocatable    :: exchange_to_use(:,:,:)      ! norb,norb,nspin
          
      allocate( exchange_to_use(norb,norb,nspin) )

      ! Here we set the first guess of the impurity selfenergy:
      ! In the beginning the Hartree term should exactly cancel the DFT Hartree
      ! contribution (usually averaged if avgHartreeFock is set), which is stored in hartree_dc
      ! Therefore, we always set the Hartree part to this (usually averaged) value
      ! For the exchange we have to decide:

       write(*,'(A)') ' Use the s_gw_dc as first s_imp but adjust high-freq. value to'
       if ( avgHartreeFock == 1 ) then
          write(*,'(A)',advance='no') ' the orbitally averaged Hartree '
       else
          write(*,'(A)',advance='no') ' the orbital-dependent Hartree'
       endif

!      if ( useSigK==1) then
         ! If we use the GW Selfenergy then we have the correct orbital dependent
         ! exchange in the GW Selfenergy and in the GW Doublecounting
         ! Therefore, we have to use the correct orbital dependent exchange,
         ! so that Sigma_GW - Sigma_DC + Sigma_imp hast again the GW exchange

!         write(*,'(A)') ' and orbital-dependent Exchange:'
!         exchange_to_use = exchange

!      else
         ! If not, we do standard LDA+DMFT, and for the DC we just subtract the averaged
         ! HF term, so we also use this as the first guess for s_imp

         write(*,'(A)') ' and orbitally averaged Exchange:'
         exchange_to_use = exchange_dc
!      endif

      do a=1,noAtoms
        do s1=1,nspin
           write(*,'(A,I2,A,I1,A,4X)',advance='no') ' Atom=',a,', Spin=',s1,':'
           do m1=1,norbPerAtom(a)
              m2 = norb_dmft(a,m1)
   
              coeffs = get_highfreq_coeff( simp(m2,m2,s1,:) ,0)
   
              write(*,'((F8.3),2X)',advance='no') hartree_dc(m2,m2,s1)+exchange_to_use(m2,m2,s1)
              do w=1,nomega
                 simp(m2,m2,s1,w) = simp(m2,m2,s1,w) - coeffs(1)   &
                                &   + hartree_dc(m2,m2,s1)+exchange_to_use(m2,m2,s1)
              enddo
           enddo ! m1
           write(*,'(A1)') ' '
         enddo ! s1
      enddo ! a Atoms

      deallocate( exchange_to_use )

   end subroutine adjust_to_hartree_fock


!  ============================================================
!  == Check and adjust the real part of the GW doublecounting to the DMFT exchange
!  == But only for the DMFT orbitals
!  ============================================================
   subroutine adjust_dc_to_exchange(s_gw_dc,exchange)
      complex(kr),intent(inout) :: s_gw_dc(:,:,:,:)     ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: exchange(:,:,:)      ! norb,norb,nspin

      integer(ki)             :: s1,m1,s2,m2,w,a
      real(kr)                :: coeffs(5)

      if (useSigK/=0) then  ! otherwise we do LDA+DMFT and s_gw_dc is zero
         do a=1,noAtoms
            do s1=1,nspin
               do m1=1,norbPerAtom(a)
                  m2 = norb_dmft(a,m1)

                  coeffs = get_highfreq_coeff( s_gw_dc(m2,m2,s1,:) ,0)

                  if ( abs(coeffs(1)-exchange(m2,m2,s1))>0.001 ) then
                     write(*,'(A,I2,A,I2,A,(F8.4),A,(F8.4),A)') ' Exchange in the s_gw_dc for DMFT orbital ',m2,   &
                                 & ', spin ',s1, ' differs to DMFT exchange by ',                           &
                                 &  abs(coeffs(1)-exchange(m2,m2,s1)),' eV, readjust to ',exchange(m2,m2,s1),' eV'
                      do w=1,nomega
                          s_gw_dc(m2,m2,s1,w) = s_gw_dc(m2,m2,s1,w) - coeffs(1)   &
                                      &   + exchange(m2,m2,s1)
                      enddo ! w
                   endif

                enddo ! m1
             enddo ! s1
         enddo ! a Atoms
      endif ! (useSigK==1)

   end subroutine adjust_dc_to_exchange
  
!  ============================================================
!  == Check total filling and adjust chemical potenital if necessary
!  ============================================================
   subroutine adjust_chem_potential(gloc,h_dft,v_xc,s_gw,s_gw_dc, &
                               &    simp,dft_hartree,dft_exchange, onlyDFT)
      complex(kr),intent(inout) :: gloc(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only the DFT Green's function
      
      real(kr)  :: fill,fillright,fillleft,muleft,muright,dmu, coeffs(5)
      integer   :: sgn,i,j, maxmuiter
      maxmuiter = 100

      ! Do a regula-falsi approach to find the mu which gives the right filling
      
      fill = total_filling(gloc)
      write(*,'(A)') ' Current total charge <=> desired charge'
      write(*,'(5X,F7.3,A,F7.3)') fill,' <=> ',nel
      
      if (abs(fill-nel)>charge_error) then
         write(*,'(A)') ' Readjust chemical potential...'

         ! determine direction where to go
         if (fill > nel) then
            dmu = -max( Uinput/2.0_kr, 1.0_kr )
            fillright = fill
            muright = mu
         else
            dmu = max( Uinput/2.0_kr, 1.0_kr )
            fillleft = fill
            muleft = mu
         endif
         
         ! get other bound 
         do i=1,maxmuiter ! set an upper limit
            mu = mu + dmu

            call get_gloc(gloc,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
            fill = total_filling(gloc)
            write(*,'(A,F10.4,A,F8.4)') ' mu= ',mu,' => total charge= ',fill
            call flush(6)

            if ( (fill-nel)*dmu > 0.0_kr ) exit  ! just look for a sign change
            if (i==maxmuiter) stop 'ERROR: Could not determine correct chemical potential!'
         enddo
         
         if (fill < nel) then ! now it's the other way round
            muleft = mu
            fillleft = fill
         else
            muright = mu
            fillright = fill
         endif
         
         ! Start regula-falsi method
         do i=1,maxmuiter ! set an upper limit
            mu = (muleft*(fillright-nel) - muright*(fillleft-nel))  &
               &                       /(fillright-fillleft)
            
            call get_gloc(gloc,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
            fill = total_filling(gloc)
            write(*,'(A,F10.4,A,F8.4)') 'mu= ',mu,' => total charge= ',fill
            call flush(6)
         
            if (abs(fill-nel)<=charge_error) exit
            
            if (fill > nel) then
               muright = mu
               fillright = fill
            else
               muleft = mu
               fillleft = fill
            endif
            if (i==maxmuiter) stop 'ERROR: Could not determine correct chemical potential!'
         enddo
         write(*,'(A)')             ' Successfully corrected the total charge:'
         write(*,'(A)')             ' Total charge <=> desired charge'
         write(*,'(6X,F7.3,A,F7.3)') fill,' <=> ',nel
         write(*,'(A,F10.4)')        ' New chemical potential: ',mu
         write(*,'(A)') ' '
      endif
   end subroutine adjust_chem_potential

!  ============================================================
!  == Symmetrize certain orbital/spin components
!  ============================================================
   subroutine symmetrize_matsubara(myname,f)
      character(len=*),intent(in) :: myname
      complex(kr),intent(inout)   :: f(:,:,:,:)        ! norb,norb,nspin,nomega
      integer(ki)                 :: s,grp,m, n, i,j,a,b, nunique, iounit=13
      complex(kr),allocatable     :: tmp(:,:,:),     & ! norb,norb,nomega
                                  &  tmp2(:,:)         ! nspin,nomega
      real(kr),allocatable       :: maxdiff(:,:,:,:)   ! norb,norb : norb,norb
      integer(ki),allocatable    :: unique(:,:)       ! norb,norb 
      integer(ki),allocatable    :: sym(:,:), symmat_in(:,:)           ! norb,norb 
      complex(kr)                :: symtmp(nomega)
      
      ! average or flip spin depending on parameters
      ! PM average
      if (spinorder==1) then
         write(*,'(A,A,A)') ' Symmetrize spin components of ',myname,' (PM)...'
         allocate( tmp(norb,norb,nomega) )
         tmp = (0.0_kr,0.0_kr)
         do s=1,nspin
            tmp = tmp + f(:,:,s,:)
         enddo
         do s=1,nspin
            f(:,:,s,:) = tmp/nspin
         enddo
         deallocate(tmp)
      ! AFM spin flip
      else if (spinorder==-100) then
         write(*,'(A,A,A)') ' Flip spin components of ',myname,' (AFM)...'
         allocate( tmp(norb,norb,nomega) )

         tmp = f(:,:,1,:)
         f(:,:,1,:) = f(:,:,2,:)
         f(:,:,2,:) = tmp

         deallocate(tmp)
      endif

! old method for only diagonals !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! orbital symmetries
!      if ( degOrbs%noGroups > 0 ) then
!         write(*,'(A,A,A)') ' Symmetrize orbital components of ',myname,' (only the diagonals) ...'
!         allocate( tmp2(nspin,nomega) )
!
!         do grp=1,degOrbs%noGroups
!
!            do m=1,degOrbs%noDegenOrbs(grp)-1
!               write(*,'(I2,A)', advance='no') degOrbs%degenOrbsIndex(grp,m),' <->'
!            enddo
!            write(*,'(I2)') degOrbs%degenOrbsIndex(grp, degOrbs%noDegenOrbs(grp) )
!
!            tmp2 = (0.0_kr,0.0_kr)
!            do m=1,degOrbs%noDegenOrbs(grp)
!               tmp2 = tmp2 + f( degOrbs%degenOrbsIndex(grp,m),degOrbs%degenOrbsIndex(grp,m),:,:)
!            enddo
!            do m=1,degOrbs%noDegenOrbs(grp)
!               f( degOrbs%degenOrbsIndex(grp,m),degOrbs%degenOrbsIndex(grp,m),:,:) = tmp2/degOrbs%noDegenOrbs(grp)
!            enddo
!         enddo
!
!         deallocate( tmp2 )
!      endif


      if (orbSym==0) then
         ! no symmetrization
      else
         allocate( maxdiff(norb,norb,norb,norb) )
         allocate( unique(norb,norb) )
         allocate( sym(norb,norb) )
         allocate( symmat_in(norb,norb) )
         maxdiff = (0.0_kr)
      
         do s=1,nspin
            if (orbSym==1) then  ! determine the symmetry yourself from looking at the differences      
      
               ! get maxdiff
               do i=1,norb
                  do j=1,norb
                     do a=1,norb
                        do b=1,norb
                           do n=1,nomega
                              if ( abs( f(i,j,s,n) - f(a,b,s,n) ) .gt. maxdiff(i,j,a,b) ) then
                                 maxdiff(i,j,a,b) = abs( f(i,j,s,n) - f(a,b,s,n) )
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

            else if (orbSym==2) then ! read the orbSym matrix and proceed as before
               open(unit=iounit,file="orbsym.dat",status="old")
               do i=1,norb
                  read(iounit,*) symmat_in(i,:)
               enddo
               rewind(iounit)
               close(iounit)

               do i=1,norb
                  do j=1,norb
                     do a=1,norb
                        do b=1,norb
                              if ( abs( symmat_in(i,j) - symmat_in(a,b) ) .gt. maxdiff(i,j,a,b) ) then
                                 maxdiff(i,j,a,b) = abs( symmat_in(i,j) - symmat_in(a,b) )
                              endif
                        enddo
                     enddo
                  enddo
               enddo
            endif ! orbSym==1 or 2

            unique = (1_ki)
            sym = (0)
            nunique = 1 ! count number of nonunique entries
   
            do i=1,norb
               do j=1,norb
   
                  if ( unique(i,j)==1_ki ) then
                     sym(i,j) = nunique
                     nunique = nunique+1
                  else
                     cycle
                  endif
   
                  do a=1,norb
                     do b=1,norb
                  
                        if ( a.eq.i .and. b.eq.j) cycle
   
                        if ( maxdiff(i,j,a,b) .lt. 0.001 .and. unique(i,j)==1_ki ) then
                           sym(a,b) = sym(i,j)
                           unique(a,b) = 0_ki
                        endif
   
                     enddo ! b
                  enddo ! a 
               enddo ! j
            enddo ! i
            
            ! there are nunique-1 unique elements because we started counting at 1
            nunique = nunique-1
   
            write(*,'(A,I4,A,A,A,I1)') ' Found ',nunique,' unique elements in ',myname,' for s = ',s
            write(*,'(A)') ' Symmetries:'
            do i=1,norb
               do j=1,norb
                  write(*,'(I4,1X)',advance='no') sym(i,j) 
               enddo
               write(*,'(A)') ' '
            enddo
   
            ! now symmetrize
            do m=1,nunique
   
               symtmp = (0.0_kr,0.0_kr)
               n = 0           ! count how many identical elements for normalization
   
               do i=1,norb
                  do j=1,norb
                     if (sym(i,j) == m) then
                        symtmp = symtmp + f(i,j,s,:)
                        n = n+1
                     endif
                  enddo
               enddo
   
               do i=1,norb
                  do j=1,norb
                     if (sym(i,j) == m) then
                        f(i,j,s,:) = symtmp / n
                     endif
                  enddo
               enddo
   
            enddo ! m=1,nunique
   
         enddo ! s
   
         deallocate( maxdiff )
         deallocate( unique )
         deallocate( sym )
         deallocate( symmat_in )

      endif

   end subroutine symmetrize_matsubara


!  ============================================================
!  == Neglect offidagonal terms (in Gloc)
!  ============================================================
   subroutine neglect_offdiag(gloc)
      complex(kr),intent(inout) :: gloc(:,:,:,:)       ! norb,norb,nspin,nomega
      integer(ki)            :: m1,m2


      if (neglectOffdiagBath==1) then

         write(*,'(A)') 'We neglect any possible offdiagonal components in Gloc before obtaining the bath!'
 
         do m1=1,norb
            do m2=1,norb

               if (m1/=m2) then
                  gloc(m1,m2,:,:) = (0.0_kr,0.0_kr)
               endif
 
            enddo
         enddo
      endif

   end subroutine neglect_offdiag


!  ============================================================
!  == Calculate the Monopole U(w) term from uloc
!  ============================================================
   subroutine get_Uavg(Uavg,uloc)
      complex(kr),intent(inout) :: Uavg(:,:)       ! noAtoms,nnu
      complex(kr),intent(in)    :: uloc(:,:,:,:) ! norb**2,norb**2,nspin,nnu

      integer(ki) :: i,j,a,m1,m2
      complex(kr) :: tmp

      Uavg = (0.0_kr,0.0_kr)
      do a=1,noAtoms
         do m1=1,norbPerAtom(a)
            do m2=1,norbPerAtom(a)
               i = norb_dmft(a,m1)
               j = norb_dmft(a,m2)
               Uavg(a,:) = Uavg(a,:) + uloc((i-1)*norb+i,(j-1)*norb+j,1,:)/norbPerAtom(a)**2
            enddo ! m2
         enddo ! m1
      enddo ! a

   end subroutine get_Uavg


!  ============================================================
!  == Correct ALPS Selfenergy
!  ============================================================
!   subroutine correct_alps_simp(simp,Uavg)
!      complex(kr),intent(inout) :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
!      real(kr),intent(in)       :: Uavg(:)        ! nnu
!      integer(ki)               :: w,m,m2,s
!      write(*,*) 'DO NOT USE correct_alps_simp !!!!!'
!      STOP 1
!
!      if (useUw/=0) then
!         write(*,'(A,F8.3)') '!!! We correct the ALPS selfenergy by -(U_bare-U(0)) =',( Uinput - Uavg(1) )
!         do m=1,norb_dmft
!            m2 = dmftOrbsIndex(m)
!            do s=1,nspin
!               do w=1,nomega
!                  simp(m2,m2,s,w) = simp(m2,m2,s,w) - ( Uinput - Uavg(1) )
!               enddo
!            enddo
!         enddo
!      endif
!   end subroutine correct_alps_simp

!  ============================================================
!  == Correct Selfenergy tail to HF from Gimp
!  ============================================================
!   subroutine correct_simp_tail(simp,gimp,uloc)
!      complex(kr),intent(inout) :: simp(:,:,:,:)        ! norb,norb,nspin,nomega
!      complex(kr),intent(in)    :: gimp(:,:,:,:)        ! norb,norb,nspin,nomega
!      complex(kr),intent(in)       :: uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu
!
!      integer(ki)    :: s,m1,m2,m,w
!      real(kr)       :: coeffs(5)
!      real(kr)       :: hartree(norb,norb,nspin)
!      real(kr)       :: exchange(norb,norb,nspin)
!
!      write(*,'(A)') '!!! We correct the ALPS selfenergy tail to HF !'
!      call get_hartree_fock(hartree,exchange,gimp,uloc)
!
!      do m=1,norb
!         do s=1,nspin
!            do w=nomega/3,nomega
!               simp(m,m,s,w) = simp(m,m,s,w) - real(simp(m,m,s,w)) + hartree(m,m,s) + exchange(m,m,s)
!            enddo
!         enddo
!      enddo
!   end subroutine correct_simp_tail

!  ==============================================================================================
!  ==============================================================================================
!  ==============================================================================================
!  == ONLY INTERACTION ROUTINES: W,U,P,Umatrix etc
!  ==============================================================================================
!  ==============================================================================================






!  ============================================================
!  == Calculate the DCA patches
!  ============================================================
   subroutine get_g0K_DCA(g0K,g0KUshift,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      complex(kr),intent(inout) :: g0K(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),intent(inout) :: g0KUshift(:,:,:)    ! npx*npy,nspin,nomega
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

      complex(kr) :: g_tmp(norb,norb)
      complex(kr) :: epsGeps(npx*npy)
      complex(kr) :: epsG(npx*npy)
      complex(kr) :: Geps(npx*npy)
      complex(kr) :: Gloc(npx*npy)
      complex(kr) :: Sloc(npx*npy)
      real(kr)    :: mu_loc(npx*npy)

      integer(ki)            :: w,s,i,j,k,m1,m2,mm1,mm2,ikx,iky, ip


      g_tmp =  (0.0_kr,0.0_kr)

      ! calculate Bath green's function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,'(A)') ' Calculate DCA coarse-grained Greens function...'
      do w=1,nomega
         call progress(w,nomega)

         do s=1,nspin
            Gloc = (0.0_kr,0.0_kr)
            Sloc = (0.0_kr,0.0_kr)
            epsGeps = (0.0_kr,0.0_kr)
            epsG = (0.0_kr,0.0_kr)
            Geps = (0.0_kr,0.0_kr)
            mu_loc = (0.0_kr)      

            do ikx=0,nkx-1
               do iky=0,nky-1
                  k = ikx*nky + iky +1 
                  !ip = patchindex(ikx,iky)
                  ip = kPatchIndex(k)
               !   write(*,'(I2,A)',advance='no') ip,' '

                  g_tmp = gfk(s,k,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
                  ! First create the inverse DFT Green's function
!                  g_tmp = ci*wn(w-1) + mu + hfield(s)  - h_dft(1,1,s,k) 
   
                  ! Then add all GW, DC and DMFT contributions
                  if (.not. onlyDFT) then
!                     g_tmp = g_tmp                  &
!         &                 + dft_hartree(1,1,s)     & ! Subtract the DFT Hartree term (no Hartree and Exchange left)
!         &                 + dft_exchange(1,1,s)    & ! Subtract the DFT Exchange term
!         &                 - simp(1,1,s,w)            ! Add the impurity DMFT Selfenergy including Hartree-Fock
                     Sloc(ip) = Sloc(ip) + simp(1,1,s,w)/nKpatch(ip)
   
                     if ( useSigK/=0 )then
!                        g_tmp = g_tmp               &
!         &                    + v_xc(1,1,s,k)          & ! Subtract V_XC from DFT
!         &                    - s_gw(1,1,s,k,w)        & ! Add the GW ex-cor. Selfenergy (add Exchange back)
!         &                    + s_gw_dc(1,1,s,w)         ! remove the doublecounting between GW and DMFT (remove local exchange)
                        Sloc(ip) = Sloc(ip) + (s_gw(1,1,s,k,w) - s_gw_dc(1,1,s,w))/nKpatch(ip)
!                     else
!                        ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                        g_tmp = g_tmp + dft_exchange(1,1,s)   ! Subtract the DFT Exchange term, it's already in simp
                     endif
                  endif ! (.not. onlyDFT)
!                 g_tmp = 1.0_kr/g_tmp
   
                 ! Then add up for every k-point the different combinations
                 Gloc(ip) = Gloc(ip) + g_tmp(1,1)/nKpatch(ip)
                 epsGeps(ip) = epsGeps(ip) + h_dft(1,1,s,k) * g_tmp(1,1) * h_dft(1,1,s,k) /nKpatch(ip)
                 epsG(ip) = epsG(ip) + h_dft(1,1,s,k) * g_tmp(1,1) /nKpatch(ip)
                 Geps(ip) = Geps(ip) + g_tmp(1,1) * h_dft(1,1,s,k) /nKpatch(ip)
   
   
                  ! calculate mu_loc here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  g_tmp(1,1) =  mu + hfield(s)  - h_dft(1,1,s,k) 
   
                  ! Then add all GW, DC and DMFT contributions
                  if (.not. onlyDFT) then
                     g_tmp(1,1) = g_tmp(1,1)                    &
         &                    + dft_hartree(1,1,s)    &  ! Subtract the DFT Hartree term (no Hartree and Exchange left)
         &                    + dft_exchange(1,1,s)       ! Subtract the DFT Exchange term
      
                     if ( useSigK/=0 )then
                        g_tmp(1,1) = g_tmp(1,1)     + v_xc(1,1,s,k)       ! Subtract V_XC from DFT to get the true noninteracting H0
!                     else
!                        ! if we do DMFT without GW, we also want to replace the local exchange by DMFT
!                        g_tmp = g_tmp + dft_exchange(1,1,s)   ! Subtract the DFT Exchange term, it's already in simp
                     endif
                  endif ! (.not. onlyDFT)
                  mu_loc(ip) = mu_loc(ip) + g_tmp(1,1)/nKpatch(ip)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
               enddo ! iky
          !     write(*,'(A)') ' '
            enddo ! ikx
          !  stop 1

            ! now we can build the noninteracting DCA gf
            do k=1,npx*npy
               if (bathMethod==0) then
                  g0KUshift(k,s,w) = 1.0_kr/( ci*wn(w-1) + mu_loc(k)-0.5*Uinput - ( epsGeps(k) - epsG(k)*Geps(k) / Gloc(k) )  )
                  g0K(k,s,w)       = 1.0_kr/( ci*wn(w-1) + mu_loc(k)           - ( epsGeps(k) - epsG(k)*Geps(k) / Gloc(k) )  )
               else
                  g0K(k,s,w)       = 1.0_kr/( 1.0_kr/Gloc(k) + Sloc(k)              )
                  g0KUshift(k,s,w) = 1.0_kr/( 1.0_kr/Gloc(k) + Sloc(k)  -0.5*Uinput  )
               endif
            enddo

         enddo ! s
      enddo ! w

   end subroutine get_g0K_DCA


!  ============================================================
!  == End of the module
!  ============================================================
end module matsfunc_ops
