!  ============================================================
!  == Handles all I/O operations
!  ============================================================

module io_ops
   use constants
   use mathfunc
   implicit none

! Return Eigenvectors of a hermitian complex matrix
 contains
!  ============================================================
!  == Read the basic parameters
!  ============================================================
   subroutine read_params()
      use params
      integer  :: iounit=10, tmp, i,j
      real(kr) :: hfield_tmp

      !  =================================================
      ! First read the parameters from the GW calculation
      !  =================================================
      open(unit=iounit,file="gw_input.dat",status="old")
      read(iounit,*) norb
!      read(iounit,*) nprodstates
      read(iounit,*) gw_nspin
      nspin = 2                ! In DMFT we always use 2
      read(iounit,*) tmp, nomega
      read(iounit,*) nnu
      read(iounit,*) gw_nkpts,nqx,nqy,nqz
!      read(iounit,*) nqpts
      read(iounit,*) beta
      read(iounit,*) mu
      close(iounit)  

      beta = beta/HartreeToEV
      mu = mu*HartreeToEV
      
      write(*,*) 'Read GW input: [eV]'
      write(*,'(A8,I6)') 'norb =',norb 
!      write(*,'(A8,I6)') 'nprodstates =',nprodstates
      write(*,'(A8,I6)') 'nspin =',gw_nspin
      write(*,'(A8,I6)') 'nomega =',nomega
      write(*,'(A8,I6)') 'nnu =',nnu
      write(*,'(A8,F8.3)') 'beta =',beta
      write(*,'(A8,F8.3)') 'mu =',mu
      write(*,'(A8,I6)') 'nqx =',nqx
      write(*,'(A8,I6)') 'nqy =',nqy
      write(*,'(A8,I6)') 'nqz =',nqz
      write(*,'(A8,I6)') 'GW no. k-points =',gw_nkpts
!      write(*,'(A8,I6)') 'GW no. q-points =',nqpts
      write(*,*) ''
      

      !  =================================================
      ! Then read the parameters for the DMFT part
      !  =================================================
      ! In the input file we read U,J and they can be interpreted in two ways:
      ! If we only use a static U, U,J are used as the bare values and the screened ones,
      ! since they are constant
      ! If we use a freq.dep. interaction, they will be used as the bare values, while
      ! the screened ones at w=0 are deduced from the U(inu_n) function that is read later

      open(unit=iounit,file="dmft_input.dat",status="old")
      read(iounit,*) Uinput,Jinput,Unninput
      read(iounit,*) readUmatrix
      read(iounit,*) U0scaleFac
      read(iounit,*) nel
      read(iounit,*) ntau
      read(iounit,*) niter
      read(iounit,*) charge_error
      read(iounit,*) mixing
      read(iounit,*) nkx
      read(iounit,*) nky
      read(iounit,*) nkz
      read(iounit,*) useSigK
      read(iounit,*) usePolK
      read(iounit,*) useUw
      read(iounit,*) DCtype
      read(iounit,*) noEquivAtoms
      call read_dmft_orbitals(iounit)
!      call read_degenerate_orbitals(iounit)
      read(iounit,*) orbSym
      read(iounit,*) bathMethod
      read(iounit,*) diagBath
      read(iounit,*) neglectOffdiagBath
      read(iounit,*) avgHartreeFock
      read(iounit,*) updateHFdc
      read(iounit,*) spinorder
      read(iounit,*) hfield_tmp
      read(iounit,*) npx
      read(iounit,*) npy
      read(iounit,*) interpDCA
      read(iounit,*) contcalc
      nkpts = nkx*nky*nkz
      close(iounit)
      
      write(*,*) 'Read DMFT input: [eV]'
      write(*,'(A,F8.3)') 'bare U_avg (F0)       = ',Uinput
      write(*,'(A,F8.3)') 'bare J_avg (F0+F2/14) = ',Jinput
      write(*,'(A,F8.3)') 'bare U_nn (dens.) = ',Unninput
      
      write(*,'(A18,I2)') 'readUmatrix =',readUmatrix
      write(*,'(A18,F8.3)') 'U0scaleFac =',U0scaleFac
      write(*,'(A18,F8.3)') 'nel =',nel
      write(*,'(A18,I6)')  'ntau =',ntau
      write(*,'(A18,I6)')  'niter =',niter
      write(*,'(A18,F8.3)') 'charge error =',charge_error
      write(*,'(A18,F7.4)') 'mixing factor =',mixing
      write(*,'(A18,I3)') 'N_kx =',nkx
      write(*,'(A18,I3)') 'N_ky =',nky
      write(*,'(A18,I3)') 'N_kz =',nkz
      write(*,'(A18,I2)') 'useSigK =',useSigK
      write(*,'(A18,I2)') 'usePolK =',usePolK
      write(*,'(A18,I2)') 'useUw =',useUw
      write(*,'(A18,I2)') 'doublecounting =',DCtype
      write(*,'(A18,I2)') 'no. equiv. atoms =',noEquivAtoms
      write(*,'(A18,I2)') 'orbital Symmetrie =',orbSym
      write(*,'(A18,I2)') 'bath method =',bathMethod
      write(*,'(A18,I2)') 'diagonalize bath =',diagBath
      write(*,'(A18,I2)') 'neglect offdiag. Elements =',neglectOffdiagBath
      write(*,'(A18,I2)') 'avgHartreeFock =',avgHartreeFock
      write(*,'(A18,I2)') 'updateHFdc =',updateHFdc
      write(*,'(A18,I2)') 'spinorder =',spinorder
      write(*,'(A18,F7.4)') 'hfield =',hfield_tmp
      write(*,'(A18,I2)') 'nKxpatches =',npx
      write(*,'(A18,I2)') 'nKypatches =',npy
      write(*,'(A18,I2)') 'interpDCA =',interpDCA
      write(*,'(A18,I2)') 'continue =',contcalc
      write(*,'(A)') ''

      ! set external magnetic field
      hfield(1) = +hfield_tmp
      hfield(2) = -hfield_tmp

      ! some basic checks for DCA !!!!!!!!!!!!!

      if ( useSigK /= 1 .and. npx*npy .gt. 1) then
         stop 'ERROR: DCA only works with useSigK=1 !'
      endif

      if (npx*npy .gt. 1) then 
         if (nkz/=1) then
            stop 'ERROR: DCA option only works with nkz=1 !'
         endif
         if (norb/=1) then
            stop 'ERROR: DCA option only works with norb=1 !'
         endif

         if ( npx>npy) then
            stop 'ERROR: DCA only works with npy>=npx for tilted patches !'
         endif

         allocate(kPatchIndex(nkpts))
         allocate(nKpatch(npx*npy))

         call fillKpatchIndex(kPatchIndex,nKpatch)
         call write_kPatchIndex("kPatchIndex.dat")

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Print info about orbitals used in DMFT
      write(*,'(A)') ' == Orbitals for DMFT: =========================================='
      write(*,'(A,I2,A,I2,A)',advance='no') ' == ',norb_dmft,' orbitals out of ',norb,' are treated in DMFT: '
      do i=1,norb_dmft-1
         write(*,'(I2,A)',advance='no') dmftOrbsIndex(i),', '
      enddo
      write(*,'(I2)') dmftOrbsIndex(norb_dmft)

      ! check if norb_dmft can be divided by noEquivAtoms
      if ( mod(norb_dmft,noEquivAtoms)==0 ) then
         norbPerAtom = norb_dmft/noEquivAtoms
         write(*,'(A,I2,A)') ' == We have ',noEquivAtoms,' equivalent atoms, so we construct a'
         write(*,'(A,I2,A)') ' == ',norbPerAtom,'-orbital impurity model for the orbitals:'
         write(*,'(A)',advance='no') ' == '
         do i=1,norbPerAtom-1
            write(*,'(I2,A)',advance='no') dmftOrbsIndex(i),', '
         enddo
         write(*,'(I2)') dmftOrbsIndex(norbPerAtom)

         do j=2,noEquivAtoms
            write(*,'(A,I2,A)',advance='no') ' == and copy the impurity Selfenergy to atom ',j,' with orbitals: '
            do i=1,norbPerAtom-1
               write(*,'(I2,A)',advance='no') dmftOrbsIndex((j-1)*norbPerAtom+i),', '
            enddo
            write(*,'(I2)') dmftOrbsIndex((j-1)*norbPerAtom+norbPerAtom)
         enddo
      else
         write(*,'(A)') 'ERROR: Number of DMFT orbitals cannot be devided by number of equiv. Atoms !!!'
         stop 1
      endif

!      write(*,'(A)') ' == -------------------------------------------------------------'
!      ! Print info about degenerate orbitals
!      write(*,'(A,I2,A)') ' == Found ', degOrbs%noGroups ,' groups of degenerate orbitals:'
!      do i=1,degOrbs%noGroups
!         write(*,'(A)',advance='no') ' == '
!         do j=1,degOrbs%noDegenOrbs(i)-1
!            write(*,'(I2,A)', advance='no') degOrbs%degenOrbsIndex(i,j),' <->'
!         enddo
!         write(*,'(I2)') degOrbs%degenOrbsIndex(i, degOrbs%noDegenOrbs(i) )
!      enddo
!      write(*,'(A)') ' ================================================================'
      write(*,'(A)') ''

      ! check if offdiagonals and diagonlization make sense
!      if ( neglectOffdiagBath==1 .and. diagBath==1  ) then
!         write(*,'(A)') 'ERROR: neglectOffdiagBath==1 and diagBath==1 does not make sense!'
!         stop 1
!      endif


      !  =================================================
      !  Adjust parameters in case of freq.dep.interactions
      !  =================================================
!      if ( useSigK /= 1 ) then
!         write(*,'(A)') ' Ignore value for usePolK, since we do not use any GW input'
!         usePolK = 0
!      endif
!      ! If we use a static U we can set U(inu_n) to be a single element, i.e. nnu=1
!      if ( useUw /= 1 ) then
!         write(*,'(A)') ' Ignore value for usePolK, since we use a static U'
!         usePolK = 0
!
!         write(*,'(A)') ' Ignore value for nnu, since we use a static U'
!         nnu=1
!      endif
      if ( usePolK==0 .and. useUw==0 .and. useSigK/=2 ) then
         write(*,'(A)') ' We set nnu=1 since no polarization or frequency dependent interaction is used'
         nnu=1

         if (readUmatrix==1) then
            write(*,'(A)') ' We will read the interaction parameters from umatrix_in.dat (will not be overwritten)'
         else
            write(*,'(A)') ' We will use a Slater/Kanamori form for the interactions, generated from Uinput,Jinput'
         endif
      else
         write(*,'(A)') ' The interaction parameters given are treated as bare values, screening will be calculated.'
      endif
      write(*,'(A)') ' '

      ! =================================================================
      ! Write a final summary of what we are actually doing =============
      ! =================================================================
      write(*,'(A)')                   ' == Final Summary: =============================================='
      if (useSigK==1) then
        write(*,'(A)',advance='yes')   ' == This will be a Sigma(k)+DMFT calculation with              =='
        write(*,'(A)',advance='yes')   ' == Sigma(k) read from gw_sigma.dat (not updated), and         =='
      else if (useSigK==2) then
        write(*,'(A)',advance='yes')   ' == This will be a Sigma(k)+DMFT calculation with              =='
        write(*,'(A)',advance='yes')   ' == Sigma(k) calculated selfconsistently from GW, and          =='
      else
        write(*,'(A)',advance='yes')   ' == This will be a normal DMFT calculation with                =='
      endif

      if (usePolK==1) then
        write(*,'(A)',advance='yes')   ' == momentum dependent Polarization Pol(k)                     =='
        write(*,'(A)',advance='yes')   ' == read from gw_pol.dat (not updated), and                    =='
      else if (usePolK==2) then
        write(*,'(A)',advance='yes')   ' == momentum dependent Polarization Pol(k)                     =='
        write(*,'(A)',advance='yes')   ' == Pol(k) calculated selfconsistently from GW, and            =='
      else
        write(*,'(A)',advance='yes')   ' == no Polarization used, and                                  =='
      endif
      if ( useUw == 1 ) then
        write(*,'(A)',advance='yes')   ' == frequency dependent interactions U(w)   with               =='
        write(*,'(A)',advance='yes')   ' == with U(w) read from Unninu2.dat                            =='
      else if ( useUw == 2 ) then
        write(*,'(A)',advance='yes')   ' == frequency dependent interactions U(w)   with               =='
        write(*,'(A)',advance='yes')   ' == with U(w) generated from GW+EDMFT                          =='
      else
        write(*,'(A)',advance='yes')   ' == static interactions U,J (at w=0 or input file)             =='
      endif
      write(*,'(A)') ' ================================================================'

      write(*,'(A)') ' '
      ! =================================================================
      ! =================================================================
      
   end subroutine read_params


!  ============================================================
!  == Read the orbitals that are used in DMFT, they might contain equivalent atoms
!  ============================================================
   subroutine read_dmft_orbitals(iounit)
        use params
        integer,intent(in)  :: iounit
        character           :: singlec
        integer(ki)         :: n,m,i, tmpi, justreadinteger
        integer,allocatable :: intarray(:)

        ! Right now we are in the line in the dmft_input.dat file with the DMFT orbitals

        read(iounit,'(A)', advance='no') singlec
        justreadinteger = -1
        m=1
        allocate( intarray(norb) )
        intarray = (0)
        norb_dmft = 0

        do while (singlec /= '#')
         ! we have read one character, and it's not #

         ! check if its an integer
         tmpi = ichar(singlec) - ichar('0')
         if ( tmpi .ge. 0 .and. tmpi .le. 9 ) then

            ! If we have just read an integer before without any separation,
            ! it must be something larger than 9, so add this
            if (justreadinteger == 1) then
               intarray(m) = intarray(m)*10 + tmpi
            else
               intarray(m) = tmpi
            endif
            norb_dmft = m

            justreadinteger = 1

         ! if not, then it's something else
         else
            ! did we read an integer before? then we can increase the index
            if(justreadinteger==1) m = m + 1

            ! now we can reset this
            justreadinteger = -1
         endif

         ! Read the next character in the line anbd repeat
         read(iounit,'(A)', advance='no') singlec
        enddo
        ! now we read all orbitals that should be treated in DMFT
        ! finish line for next read
        read(iounit,'(A)')

        ! Check if more orbitals are specified than provided by the GW/DFT input
        if (norb_dmft>norb) then
           write(*,'(A)') "ERROR: More DMFT orbitals specified than available !!!"
           stop 1
        endif
        ! Check if there are at least one orbital
        if (norb_dmft==0) then
           write(*,'(A)') "ERROR: No orbitals specified for DMFT !!!"
           stop 1
        endif
        ! check if all orbitals actually exist
        do n=1,norb_dmft
           if ( intarray(n)>norb .or. intarray(n)<=0 ) then
              write(*,'(A,I3,A)') "ERROR: Specified correlated orbital ",intarray(n)," does not exist !!!"
              stop 1
           endif
        enddo

        ! allocate and copy into proper array
        allocate( dmftOrbsIndex( norb_dmft ) )

      !  ! first sort the array by doing a stupid loop from 1 to norbs (maximum number possible)
      !  ! this scales as norb*norb_dmft, which is bad but doesnt matter
      !  i = 1
      !  do m=1,norb
      !     do n=1,norb_dmft
      !        if ( intarray(n)==m ) then
      !           dmftOrbsIndex(i) = m
      !           i = i+1
      !        endif
      !     enddo
      !  enddo

        dmftOrbsIndex = intarray(1:norb_dmft)
        deallocate( intarray )

      !  ! check if the user specified an orbital twice
      !  do n=2,norb_dmft
      !     if ( dmftOrbsIndex(n-1) == dmftOrbsIndex(n) ) then
      !        write(*,'(A,I3,A)') "ERROR: Specified orbital ",dmftOrbsIndex(n)," more than once as correlated!!!"
      !        stop 1
      !     endif
      !  enddo

        ! Finally done
    end subroutine read_dmft_orbitals


!  ============================================================
!  == Read the degenerate DMFT orbitals
!  ============================================================
!   subroutine read_degenerate_orbitals(iounit)
!        use params
!        integer,intent(in)  :: iounit
!        character           :: singlec
!        integer(ki)         :: n,m,i, tmpi, justreadinteger
!        integer(ki)         :: degenindex, degengroup
!        integer,allocatable :: degenarraytmp(:,:,:)
!
!        allocate( degenarraytmp(norb,norb,2) )
!        degenarraytmp = (0)
!        degenindex = 1
!        degengroup = 1
!        justreadinteger = -1
!
!        read(iounit,'(A)', advance='no') singlec
!        do while (singlec /= '#')
!           ! we have read one character, and it's not #
!
!           ! check if its an integer
!           tmpi = ichar(singlec) - ichar('0')
!           if ( tmpi .ge. 0 .and. tmpi .le. 9 ) then
!
!              ! If we have just read an integer before without any separation,
!              ! it must be something larger than 9, so add this
!              if (justreadinteger == 1) then
!                 degenarraytmp(degengroup,degenindex,1) = degenarraytmp(degengroup,degenindex,1)*10 + tmpi
!              else
!                 degenarraytmp(degengroup,degenindex,1) = tmpi
!              endif
!              degenarraytmp(degengroup,1,2) = degenindex
!
!              justreadinteger = 1
!
!           ! if not, then it's something else
!           else
!              ! did we read an integer before? then we can increase the index
!              if(justreadinteger==1) degenindex = degenindex + 1
!
!              ! now we can reset this
!              justreadinteger = -1
!
!              ! check if doublecolon, then we have new group
!              ! if degenindex 1, no orbitals have been specified yet
!              if (singlec == ':' .and. degenindex>1) then
!                 degengroup = degengroup +1
!                 degenindex = 1
!              endif
!           endif
!
!           ! Read the next character in the line anbd repeat
!           read(iounit,'(A)', advance='no') singlec
!        enddo
!
!        ! Check if the last group is empty
!        if (degenindex==1) degengroup = degengroup-1
!
!        ! Check if specified orbitals acutally exist
!        do n=1,degengroup
!           do m=1,degenarraytmp(n,1,2)
!              tmpi = 0
!              ! is this orbital also appearing in the DMFT orbitals?
!            !  do i=1,norb_dmft
!            !     if ( degenarraytmp(n,m,1) == dmftOrbsIndex(i) ) tmpi = 1
!              do i=1,norb
!                 if ( degenarraytmp(n,m,1) == i ) tmpi = 1
!              enddo
!              if ( tmpi == 0 ) then
!                 write(*,'(A,I2,A)') 'ERROR: Specified degenerate orbital ',degenarraytmp(n,m,1),  &
!                              &      ', which is not in the orbital list !!!'
!                 stop 1
!              endif
!           enddo
!        enddo
!
!        ! are there degenerate orbitals at all?
!        if ( degengroup==0 ) then
!           write(*,'(A)') 'No degenerate orbitals specified'
!           degOrbs%noGroups = 0
!        else
!            ! Yes there are
!            ! Allocate with norb_dmft, this might be too much but we dont care
!            degOrbs%noGroups = degengroup
!            allocate( degOrbs%noDegenOrbs(degengroup) )         ! groups
!            allocate( degOrbs%degenOrbsIndex(degengroup,norb_dmft) ) ! groups, orbs of group
!
!            do n=1,degOrbs%noGroups
!               ! copy all orbitals in group n into the array, their number is given in degenarraytmp(n,1,2)
!               degOrbs%noDegenOrbs(n) = degenarraytmp(n,1,2)
!               degOrbs%degenOrbsIndex(n,1:degOrbs%noDegenOrbs(n)) = degenarraytmp(n, 1:degOrbs%noDegenOrbs(n), 1)
!            enddo
!        endif
!        read(iounit,'(A)')
!        deallocate( degenarraytmp )
!
!   end subroutine read_degenerate_orbitals



!  ============================================================
!  == Read the input from GW
!  ============================================================
   subroutine read_gw_input(h_dft,v_xc,S_gw,V_gw,P_gw)
      use params
      use fourier_interpolation
      complex(kr),intent(inout) :: h_dft(:,:,:,:)        ! norb,norb,nspin,nkpts
      complex(kr),intent(inout) :: v_xc(:,:,:,:)         ! norb,norb,nspin,nkpts
      complex(kr),intent(inout) :: S_gw(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(inout) :: V_gw(:,:,:,:,:)       ! norb**2,norb**2,nspin,nkpts,nnu
      complex(kr),intent(inout) :: P_gw(:,:,:,:,:)       ! norb**2,norb**2,nspin,nkpts,nnu

      integer              :: m1,m2,s,w,i,iounit=10
      real(kr),allocatable :: tmp(:)             ! 3+gw_nspin*norb*norb
      real(kr),allocatable :: tmp2(:)            ! 1+2*gw_nspin*norb*norb
      real(kr),allocatable :: tmpprod1(:)        ! 1+gw_nspin*norb**2*norb**2
      real(kr),allocatable :: tmpprod2(:)        ! 2+gw_nspin*norb**2*norb**2
      real(kr),allocatable :: tmpprod3(:)        ! norb**2
      real(kr)             :: tmptmp             ! Just temp
      complex(kr),allocatable :: eigvectors(:,:)        ! for diagonalizing the tight-binding hamiltonian

      ! The following arrays are for reading in the GW data on a coarse grid with gw_nkpts k-points 
      ! After reading, the data will be interpolated onto a dense grid
      ! with nkpts k-points and stored in the corresponding arrays (see above)
      ! "_q" will be attached to the objects corresponding to the coarse grid and qx,qy,qz is used
      integer                 :: iq
      complex(kr),allocatable :: h_dft_q(:,:,:,:)   ! norb,norb,nspin,gw_nkpts
      complex(kr),allocatable :: v_xc_q(:,:,:,:)    ! norb,norb,nspin,gw_nkpts
      complex(kr),allocatable :: S_gw_q(:,:,:,:,:)  ! norb,norb,nspin,gw_nkpts,nnu
      allocate( h_dft_q  (norb,norb,nspin,gw_nkpts) )
      
      allocate(tmp(3+2*gw_nspin*norb*norb))
      allocate(tmp2(1+2*gw_nspin*norb*norb))
      allocate(tmpprod1(1+gw_nspin*norb**2*norb**2))
      allocate(tmpprod2(2+gw_nspin*norb**2*norb**2))
      allocate(tmpprod3(norb**2))
      allocate(eigvectors(norb,norb))

      ! Init the arrays that contain the final data
      h_dft            = (0.0_kr,0.0_kr)

      ! Init arrays for reading in the data
      h_dft_q     = (0.0_kr,0.0_kr)

      
    !  ============================================================
    !  == Read the DFT Hamiltonian
    !  ============================================================
      write(*,'(A)') ' Read the DFT Hamiltonian'
      open(unit=iounit,file="gw_h_dft.dat",status="old")
      read(iounit,*) ! comment line

      do iq=1,gw_nkpts
         read(iounit,*) tmptmp,tmp  ! first number is just k index
         gw_kvecs(iq,1) = tmp(1)
         gw_kvecs(iq,2) = tmp(2)
         gw_kvecs(iq,3) = tmp(3)
      
         do s=0,gw_nspin-1
            do m2=0,norb-1
               do m1=0,norb-1
                  h_dft_q(m1+1,m2+1,s+1,iq) =      tmp(3 + 2*s*norb**2 + 2*m1*norb + 2*m2 +1)      &
                                          &  +ci*( tmp(3 + 2*s*norb**2 + 2*m1*norb + 2*m2 +2) )
               enddo !m2 loop        
            enddo !m1 loop
         enddo !s loop
      enddo !iq loop
      close(iounit)
      if (gw_nspin==1) h_dft_q(:,:,2,:) = h_dft_q(:,:,1,:)  ! Just copy to other spin
      h_dft_q = h_dft_q*HartreeToEV

!  ==================================================================================
!  == Only if useSigK==1 we read the V_XC, GW-Selfenergy
!  ==================================================================================
      if (useSigK==1) then
        !  ============================================================
        !  == Read the exchange-correlation potential V_XC
        !  ============================================================
          allocate( v_xc_q(norb,norb,nspin,gw_nkpts) )        ; v_xc_q  = (0.0_kr,0.0_kr)
          allocate( S_gw_q(norb,norb,nspin,gw_nkpts,nomega) ) ; S_gw_q  = (0.0_kr,0.0_kr)

          write(*,'(A)') ' Read the exchange correlation potential'
          open(unit=iounit,file="gw_v_xc.dat",status="old")
          read(iounit,*) ! comment line
          do iq=1,gw_nkpts
             read(iounit,*) tmptmp,tmp   ! first number is just k index

             do s=0,gw_nspin-1
                do m2=0,norb-1
                   do m1=0,norb-1
                      v_xc_q(m1+1,m2+1,s+1,iq) =       tmp(3 + 2*s*norb**2 + 2*m1*norb + 2*m2 +1) &
                                               &+ ci*( tmp(3 + 2*s*norb**2 + 2*m1*norb + 2*m2 +2) )
                   enddo !m2 loop
                enddo !m1 loop
             enddo !s loop
          enddo !iq loop
          close(iounit)
          if (gw_nspin==1) v_xc_q(:,:,2,:) = v_xc_q(:,:,1,:)  ! Just copy to other spin
          v_xc_q = v_xc_q*HartreeToEV

        !  ============================================================
        !  == Read the GW Selfenergy
        !  ============================================================
          write(*,'(A)') ' Read the GW Selfenergy Sigma_XC'
          open(unit=iounit,file="gw_sigma.dat",status="old")
          read(iounit,*) ! comment line
          do iq=1,gw_nkpts
             read(iounit,*) ! comment line
             do w=1,nomega
                read(iounit,*) tmp2 ! read full line

                do s=0,gw_nspin-1
                   do m2=0,norb-1
                      do m1=0,norb-1
                         S_gw_q(m1+1,m2+1,s+1,iq,w) =      &
                &       tmp2(1 + 2*s*norb**2 + 2*m1*norb + 2*m2 +1)  &
                & +ci*( tmp2(1 + 2*s*norb**2 + 2*m1*norb + 2*m2 +2) )

                      enddo !m2 loop
                   enddo !m1 loop
                enddo !s loop
             enddo ! w loop
             read (iounit, '()', advance='yes')  ! advance empty line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             do w=1,nomega
!                do s=0,gw_nspin-1
!                   do m2=0,norb-1
!                      do m1=0,norb-1
!
!                        if (m1 /= m2) then
!                           S_gw_q(m1+1,m2+1,s+1,iq,w) = S_gw_q(m1+1,m2+1,s+1,iq,w) &
!                                                      -ci*aimag(S_gw_q(m1+1,m2+1,s+1,iq,nomega))
!                        endif 
!
!                      enddo !m2 loop
!                   enddo !m1 loop
!                enddo !s loop
!             enddo ! w loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          enddo !iq loop
          close(iounit)
          if (gw_nspin==1) S_gw_q(:,:,2,:,:) = S_gw_q(:,:,1,:,:)  ! Just copy to other spin
          S_gw_q = S_gw_q*HartreeToEV

       ! ======================================================
       endif ! if (useSigK==1) then1
       ! if useSigK==2 we dont need to read since we calculate Sigma ourselves
       ! ======================================================

 ! ==========================================================
 ! The bare Coulomb interaction, polarization 
 ! ==========================================================

    if ( usePolK==1 ) then
        !  ============================================================
        !  == Read the bare Coulomb interaction V
        !  ============================================================
          write(*,'(A)') ' Read the bare Coulomb interaction V(q) from gw_vcoul.dat'
          stop " READING P AND V FOR usePolK==1 is not implemented yet!"

!          open(unit=iounit,file="gw_vcoul.dat",status="old")
!          read(iounit,*) ! comment line
!          do iq=1,nqpts
!            read(iounit,*) tmpprod1 ! read full line
!
!            do s=0,gw_nspin-1
!               do m2=0,norb**2-1
!                  do m1=0,norb**2-1
!                     v_coul_q(m1+1,m2+1,s+1,iq) =      &
!                 &          tmpprod1(1 + s*norb**2**2 + m1*norb**2 + m2 +1)
!
!                  enddo !m2 loop
!               enddo !m1 loop
!
!               ! Now we have a matrix in the product basis for given q,s
!               ! invert it to get V^-1
!               vcoul_tmp = v_coul_q(:,:,s+1,iq)
!
!               call DGETRF(norb**2, norb**2, vcoul_tmp, norb**2, ipiv, info_lpck)
!               if (info_lpck /= 0) then
!                  write(*,*)'ERROR: Coulomb matrix in product basis is numerically singular! Return value',&
!                            & info_lpck
!                  stop
!               end if
!               call DGETRI(norb**2, vcoul_tmp, norb**2, ipiv, work, norb**2, info_lpck)
!              if (info_lpck /= 0) then
!                 stop 'Coulomb matrix in product basis inversion failed!'
!              end if
!
!              v_coulinv_q(:,:,s+1,iq) = vcoul_tmp
!
!            enddo !s loop
!          enddo !iq loop
!          close(iounit)
!          if (gw_nspin==1) then
!             v_coul_q(:,:,2,:) = v_coul_q(:,:,1,:)  ! Just copy to other spin
!             v_coulinv_q(:,:,2,:) = v_coulinv_q(:,:,1,:)  ! Just copy to other spin
!          endif
!          v_coul_q = v_coul_q*HartreeToEV
!          v_coulinv_q = v_coulinv_q/HartreeToEV ! divide by, because this is V^-1
!
!        !  ============================================================
!        !  == Read the GW polarization
!        !  ============================================================
!          write(*,'(A)') ' Read the GW polarization P(q,ivn)'
!          open(unit=iounit,file="gw_pol.dat",status="old")
!          read(iounit,*) ! comment line
!          do iq=1,nqpts
!             read(iounit,*) ! comment line
!             do w=1,nnu
!                read(iounit,*) tmpprod2 ! read full line
!
!                do s=0,gw_nspin-1
!                   do m2=0,norb**2-1
!                      do m1=0,norb**2-1
!                         p_gw_q(m1+1,m2+1,s+1,iq,w) =      &
!                  &        tmpprod2(2 + s*norb**2**2 + m1*norb**2 + m2 +1)
!
!                      enddo !m2 loop
!                   enddo !m1 loop
!                enddo !s loop
!             enddo ! w loop
!             read (iounit, '()', advance='yes')  ! advance empty line
!          enddo !iq loop
!          close(iounit)
!          if (gw_nspin==1) p_gw_q(:,:,2,:,:) = p_gw_q(:,:,1,:,:)  ! Just copy to other spin
!          p_gw_q = p_gw_q/HartreeToEV
!
!        !  ============================================================
!        !  == Read the Overlap matrix from the Product to the 2particle basis
!        !  ============================================================
!          write(*,'(A)') ' Read the overlap matrix Prod2PartOverlap'
!          open(unit=iounit,file="gw_overlap.dat",status="old")
!          read(iounit,*) ! comment line
!
!          do m1=0,norb*norb-1
!             read(iounit,*) tmpprod3 ! read full line
!             do m2=0,norb**2-1
!                Prod2PartOverlap(m1+1,m2+1) = tmpprod3(m2+1)
!             enddo !npordstates loop
!          enddo !norb loop
!          close(iounit)
!
     ! ==========================================
     endif ! usePolK
     ! ==========================================

     !  ============================================================
     !  == Now to the Fourier/Cubic interpolation if necessary
     !  ============================================================
      if (gw_nkpts .ne. nkpts ) then
         ! do the interpolation to a finer grid
         write(*,'(A,I7,A,I7,A)') ' Interpolate the GW input from ',gw_nkpts,' k-points to ',nkpts,' ...'

         write(*,'(A)') '!!! =========================================================================== !!!'
         write(*,'(A)') '!!! WARNING: No interpolation implemented on V_coul and P_gw !!!'
         write(*,'(A)') '!!! =========================================================================== !!!'
!         p_gw = p_gw_q

         call flush(6)
         
         call interpolate_kmatrix_cubic(h_dft,h_dft_q,nkx,nky,nkz,nqx,nqy,nqz)
         write(*,'(A)') ' Interpolation for h_dft finished'
         
         if (useSigK==1) then ! we only need to interpolate when we read Sigma from file
             call interpolate_kmatrix_cubic(v_xc,v_xc_q,nkx,nky,nkz,nqx,nqy,nqz)
             write(*,'(A)') ' Interpolation for v_xc finished'

             do w=1,nomega
                call interpolate_kmatrix_cubic(S_gw(:,:,:,:,w), S_gw_q(:,:,:,:,w),nkx,nky,nkz,nqx,nqy,nqz )
             enddo
             write(*,'(A)') ' Interpolation for S_gw finished'

             deallocate( v_xc_q )
             deallocate( S_gw_q )

         endif ! useSigK==1

      else ! no interpolation needed
         ! otherwise, just copy the input to the proper arrys
         h_dft = h_dft_q
!         p_gw = p_gw_q
         if (useSigK==1) then 
            v_xc = v_xc_q
            S_gw = S_gw_q            
            deallocate( v_xc_q )
            deallocate( S_gw_q )
         endif
      endif
      write(*,'(A)') ' '

      !  ============================================================
      !  == Done reading input
      !  ============================================================
      
      deallocate( h_dft_q )
      deallocate(tmp)
      deallocate(tmp2)      
      deallocate(tmpprod1)
      deallocate(tmpprod2)
      deallocate(tmpprod3)
      deallocate(eigvectors)
   end subroutine read_gw_input
   
!  ============================================================
!  == Read a local Matsubara function from file
!  == Read all norb orbitals
!  ============================================================
   subroutine read_loc_matsfunc(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:,:)  ! norb,norb,nspin,nomega

      integer(ki) :: w,s,m,a
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(1+2*nspin*norb) ! temparray for readin
      complex(kr)           :: tmp
      
      ! Just read the diagonal components!
      
      f = (0.0_kr, 0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,1,:))
         read(iounit,*) readindata
         do s=0,nspin-1
            do m=0,size(f(:,1,s+1,w))-1
              f(m+1,m+1,s+1,w) = readindata(1+s*norb*2+m*2+1) + ci*readindata(1+s*norb*2+m*2+2)

              ! Do sanity check
              if ( isnan_c( f(m+1,m+1,s+1,w) ) ) then
                 write(*,'(A,A)') "ERROR: Encountered NaN when reading ",filename
                 stop 1
              endif 

            enddo ! m loop
         enddo ! s loop
      enddo ! w loop
      close(iounit)
   end subroutine read_loc_matsfunc

!  ============================================================
!  == Read a local Matsubara function from file in vector forn
!  ============================================================
   subroutine read_loc_matsfuncvec(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:)  ! ndim,nspin,nomega

      integer(ki) :: w,s,m, ndim
      integer(ki),parameter :: iounit = 10
      real(kr),allocatable  :: readindata(:) ! temparray for readin
      complex(kr)           :: tmp

      ndim = size( f(:,1,1) )
      allocate(readindata(1+2*nspin*ndim))

      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,:))
         read(iounit,*) readindata
         do s=0,nspin-1
            do m=0,ndim-1

                 f(m+1,s+1,w)          &
                         &        =     readindata(1 + s*ndim*2 +  m*2 +1)  &
                         &         + ci*readindata(1 + s*ndim*2 +  m*2 +2)

                 ! Do sanity check
                 if ( isnan_c( f(m+1,s+1,w) ) ) then
                    write(*,'(A,A)') "ERROR: Encountered NaN when reading ",filename
                    stop 1
                 endif 

            enddo !m
         enddo ! s lopp
      enddo ! w loop
      close(iounit)
      deallocate(readindata)
   end subroutine read_loc_matsfuncvec

!  ============================================================
!  == Read a local Matsubara function from file in matrix forn
!  == Copy to the other equiv. atoms when necessary
!  ============================================================
   subroutine read_loc_matsfunc_matrix(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:,:)  ! ndim,ndim,nspin,nomega

      integer(ki) :: w,s,m1,m2,a, ndim
      integer(ki),parameter :: iounit = 10
      real(kr),allocatable  :: readindata(:) ! temparray for readin
      complex(kr)           :: tmp

      ndim = size( f(:,1,1,1) )
      allocate(readindata(1+2*nspin*ndim*ndim))

      ! f = (0.0_kr, 0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,1,:))
         read(iounit,*) readindata
!         write(*,*) 'readindata at w=',w,":"
!         write(*,*) readindata
         do s=0,nspin-1
            do m1=0,ndim-1
               do m2=0,ndim-1

                 f(m1+1,m2+1,s+1,w)          &
                         &        =     readindata(1 + s*ndim*ndim*2 +  m1*ndim*2 + m2*2 +1)  &
                         &         + ci*readindata(1 + s*ndim*ndim*2 +  m1*ndim*2 + m2*2 +2)

                 ! Do sanity check
                 if ( isnan_c( f(m1+1,m2+1,s+1,w) ) ) then
                    write(*,'(A,A)') "ERROR: Encountered NaN when reading ",filename
                    stop 1
                 endif 

               enddo
            enddo


         enddo ! s lopp
      enddo ! w loop
      close(iounit)
      deallocate(readindata)
   end subroutine read_loc_matsfunc_matrix

!  ============================================================
!  == Read a local Matsubara function from file, from the DMFT solver
!  == Read only norbPerAtom orbitals !!!
!  == Copy to the other equiv. atoms when necessary
!  ============================================================
   subroutine read_loc_matsfunc_solver(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:,:)  ! norb,norb,nspin,nomega

      integer(ki) :: w,s,m,a
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(1+2*nspin*norbPerAtom) ! temparray for readin
      complex(kr)           :: tmp

      f = (0.0_kr, 0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,1,:))
         read(iounit,*) readindata
         do s=0,nspin-1
            do m=0,norbPerAtom-1
              ! For the impsolver, spins run first, then orbitals!
              f(dmftOrbsIndex(m+1),dmftOrbsIndex(m+1),s+1,w)          &
                      &        = readindata(1+m*nspin*2+s*2+1) + ci*readindata(1+m*nspin*2+s*2+2)

                 ! Do sanity check
                 if ( isnan_c( f(dmftOrbsIndex(m+1),dmftOrbsIndex(m+1),s+1,w) ) ) then
                    write(*,'(A,A)') "ERROR: Encountered NaN when reading solver output !!!"
                    stop 1
                 endif 

            enddo

            ! now copy to the other equivalent atoms
            do a=1,noEquivAtoms-1
               do m=1,norbPerAtom
                   f( dmftOrbsIndex(a*norbPerAtom+m), dmftOrbsIndex(a*norbPerAtom+m) ,s+1,w)   &
                                   &  = f( dmftOrbsIndex(m),dmftOrbsIndex(m) ,s+1,w)
               enddo
            enddo

         enddo ! s lopp
      enddo ! w loop
      close(iounit)
   end subroutine read_loc_matsfunc_solver

!  ============================================================
!  == Read a local Matsubara function from file, from the DMFT solver
!  == Read only norbPerAtom orbitals !!!
!  == Copy to the other equiv. atoms when necessary
!  ============================================================
   subroutine read_loc_matsfunc_matrix_solver(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:,:)  ! norb,norb,nspin,nomega

      integer(ki) :: w,s,m1,m2,a
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(1+2*nspin*norbPerAtom*norbPerAtom) ! temparray for readin
      complex(kr)           :: tmp, Ftmp(norbPerAtom,norbPerAtom,nspin)

      write(*,'(A,A)') 'Read solver output ',filename
      ! Do not overwrite the full matrix !!!
      ! f = (0.0_kr, 0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,1,:))
         read(iounit,*) readindata
         do s=0,nspin-1
            do m1=0,norbPerAtom-1
               do m2=0,norbPerAtom-1
                 Ftmp(m1+1,m2+1,s+1) = readindata(1 + s*norbPerAtom**2*2 + m1*norbPerAtom*2 + m2*2 + 1)  &
                                 & +ci*readindata(1 + s*norbPerAtom**2*2 + m1*norbPerAtom*2 + m2*2 + 2)

                 ! Do sanity check
                 if ( isnan_c(Ftmp(m1+1,m2+1,s+1)) ) then
                    write(*,'(A,A)') "ERROR: Encountered NaN when reading solver output !!!"
                    stop 1
                 endif

               enddo ! m2 
            enddo ! m1

            ! Now transform back, this also works when we don't want to transform, since Utrafo = identity
            m1 = dmftOrbsIndex(1)
            m2 = dmftOrbsIndex(norbPerAtom)

            f(m1:m2,m1:m2,s+1,w) = matmul( Utrafo(:,:,s+1), matmul( Ftmp(:,:,s+1), UtrafoT(:,:,s+1) ) )

            ! now copy to the other equivalent atoms
            do a=1,noEquivAtoms-1
               do m1=1,norbPerAtom
                  do m2=1,norbPerAtom
                      f( dmftOrbsIndex(a*norbPerAtom+m1), dmftOrbsIndex(a*norbPerAtom+m2) ,s+1,w)   &
                                      &  = f( dmftOrbsIndex(m1),dmftOrbsIndex(m2) ,s+1,w)
                  enddo
               enddo 
            enddo ! a loop

         enddo ! s lopp
      enddo ! w loop
      close(iounit)
   end subroutine read_loc_matsfunc_matrix_solver

!  ============================================================
!  == Read a local Matsubara function from file, from the DMFT solver
!  == THIS IS FOR A 3D HUBBARD CUBE WITHIN THE DIMER APPROXIMATION
!  ============================================================
   subroutine read_loc_matsfunc_matrix_solver_dimercube(filename,f)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: f(:,:,:,:)  ! norb,norb,nspin,nomega

      integer(ki) :: w,s,m1,m2,a
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(1+2*nspin*norbPerAtom*norbPerAtom) ! temparray for readin
      complex(kr)           :: tmp, Ftmp(norbPerAtom,norbPerAtom,nspin)

      ! Just read and OVERWRITE the diagonal components!
      ! In this case we leave all other elements as they are
      ! Because in Gimp the non-DMFT orbitals will be equal to
      ! Gloc, so we dont want to set them to zero
      ! But be careful to initialize everything correctly

      write(*,'(A,A)') 'Read solver output ',filename
      ! f = (0.0_kr, 0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(f(1,1,1,:))
         read(iounit,*) readindata
         do s=0,nspin-1
            do m1=0,norbPerAtom-1
               do m2=0,norbPerAtom-1
                 ! For the impsolver, spins run first, then orbitals!
                 Ftmp(m1+1,m2+1,s+1) = readindata(1+m1*norbPerAtom*nspin*2+m2*nspin*2+s*2+1)  &
                                 & +ci*readindata(1+m1*norbPerAtom*nspin*2+m2*nspin*2+s*2+2)

               enddo ! m2 
            enddo ! m1

            ! Now transform back
            m1 = dmftOrbsIndex(1)
            m2 = dmftOrbsIndex(norbPerAtom)

            f(m1:m2,m1:m2,s+1,w) = matmul( Utrafo(:,:,s+1), matmul( Ftmp(:,:,s+1), UtrafoT(:,:,s+1) ) )


         enddo ! s loop

         ! now copy to the other equivalent atoms

         ! ONSITE TERMS
         do m1=1,2
            do m2=1,2
               do s=1,2

                  ! A site            
                  f( m1+6,m2+6,s,w)   = f( m1,m2 ,s,w)
                  f( m1+10,m2+10,s,w) = f( m1,m2 ,s,w)
                  f( m1+12,m2+12,s,w) = f( m1,m2 ,s,w)

                  ! B site
                  f( m1+4,m2+4,s,w)   = f( m1+2,m2+2 ,s,w)
                  f( m1+8,m2+8,s,w)   = f( m1+2,m2+2 ,s,w)
                  f( m1+14,m2+14,s,w) = f( m1+2,m2+2 ,s,w)

               enddo
            enddo 
         enddo 

         ! Intersite TERMS
         do m1=1,2
            do m2=1,2
               do s=1,2

                  f( m1+0,m2+4,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+4,m2+0,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+0,m2+8,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+8,m2+0,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+6,m2+2,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+2,m2+6,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+10,m2+2,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+2,m2+10,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+6,m2+4,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+4,m2+6,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+12,m2+4,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+4,m2+12,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+6,m2+14,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+14,m2+6,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+10,m2+8,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+8,m2+10,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+12,m2+8,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+8,m2+12,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+10,m2+14,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+14,m2+10,s,w) = f( m1+2,m2+0 ,s,w)

                  f( m1+12,m2+14,s,w) = f( m1+0,m2+2 ,s,w)
                  f( m1+14,m2+12,s,w) = f( m1+2,m2+0 ,s,w)

               enddo
            enddo 
         enddo 

      enddo ! w loop
      close(iounit)
   end subroutine read_loc_matsfunc_matrix_solver_dimercube

!  ============================================================
!  == Read the frequency dependent density-density part of U(w)
!  ============================================================
   subroutine read_udens_w(filename,udensw)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(inout)    :: udensw(:)  ! nnu

      integer(ki) :: w
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(2) ! temparray for readin
      
      udensw = (0.0_kr,0.0_kr)
      open(unit=iounit,file=filename,status="old")
      do w=1,size(udensw)
         read(iounit,*) readindata  ! frequency, U(w) value
         udensw(w) = readindata(2)  ! this function has no imaginary part
      enddo ! w loop
      close(iounit)
   end subroutine read_udens_w

!  ============================================================
!  == Read the impurity polarization from file
!  ============================================================
   subroutine read_wimp_pimp(fileW,fileT,pimp,wimp,uloc)
      use params
      use matsfunc_ops      
      character(len=*), intent(in) :: fileW
      character(len=*), intent(in) :: fileT
      complex(kr),intent(inout)    :: pimp(:,:,:,:)  ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(inout)    :: wimp(:,:,:,:)  ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(in)       :: uloc(:,:,:,:)  ! norb**2,norb**2,nspin,nnu
      integer(ki) :: w,s,m,m2
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(4), chi, coeffs(5), corrfac ! temparray for readin
      character(len=10)     :: tmp
      
      pimp = (0.0_kr)
      wimp = (0.0_kr)

      ! check the layout of the solver output, spins or orbitals first?

      if ( useUw/=0 .and. usePolK/=0 ) then
         if (norbPerAtom>1) then
            stop "ERROR: read_wimp_pimp not implemented for norbPerAtom>1 !"
         endif
      
          ! First get Chi(0)
          open(unit=iounit,file=fileT,status="old")
             read(iounit,*) tmp
             read(iounit,*) readindata
             corrfac = beta*( readindata(2)+readindata(4)  )**2
          close(iounit)

          open(unit=iounit,file=fileW,status="old")
          ! skip first comment line           
             read(iounit,*) tmp

          do w=1,size(pimp(1,1,1,:))
             read(iounit,*) readindata
             chi = readindata(2) + 2*readindata(3) + readindata(4)

             if (w==1) then
                !chi = chi - beta*( filling(gimp(1,1,1,:)) + filling(gimp(1,1,2,:))  )**2  ! subtract the constant density-density part
                chi = chi - corrfac  ! subtract the constant density-density part
             endif

             pimp(1,1,1,w) = chi/( chi*uloc(1,1,1,w) - 1 )
          enddo
           ! interpolate to zero to replace the w=0 value
           !coeffs = get_lowfreq_coeff_bos(pimp(1,1,1,2:),5,6)
           !pimp(1,1,1,1) = coeffs(1)

          pimp(1,1,2,:) = pimp(1,1,1,:) ! copy to other spin

          close(iounit)

         ! Now we have p_imp
          do w=1,size(pimp(1,1,1,:))
            do s=1,nspin
               wimp(1,1,s,w) = uloc(1,1,s,w) / ( 1 - pimp(1,1,s,w) * uloc(1,1,s,w) )
            enddo ! s
          enddo ! w

         ! now copy to other atom ; We still assume norbPerAtom==1
         ! AND WE only replace the diagonals
         do m=2,noEquivAtoms
             m2 = dmftOrbsIndex(m)-1
             pimp(m2*norb+m2+1,m2*norb+m2+1,:,:) = pimp(1,1,:,:)
             wimp(m2*norb+m2+1,m2*norb+m2+1,:,:) = wimp(1,1,:,:)
         enddo

       endif

   end subroutine read_wimp_pimp

   
!  ============================================================
!  == Write a local Green's function vector into file
!  ============================================================
   subroutine write_loc_matsfuncvec(filename,func)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:)  ! norb,nspin,nomega
      integer(ki) :: w,s,m
      integer(ki),parameter :: iounit = 10

      ! Just write the diagonal components!

      open(unit=iounit,file=filename,status="unknown")
    !  write(iounit,*) '# wn: orbital runs first, then spin'
      do w=1,size(func(1,1,:))
         write(iounit,'((ES23.16),(4X))',advance='no') wn(w-1)
         do s=1,nspin
            do m=1,size(func(:,s,w))
               write(iounit,'(2(ES23.16,2X),(4X))',advance='no')  & 
                                 & real(func(m,s,w)), aimag(func(m,s,w))
            enddo
         enddo
         write(iounit,'(A1)') ' '
      enddo
      close(iounit)
   end subroutine write_loc_matsfuncvec


!  ============================================================
!  == Write a local Green's function into file
!  ============================================================
   subroutine write_loc_matsfunc(filename,func)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:)  ! norb,norb,nspin,nomega
      integer(ki) :: w,s,m
      integer(ki),parameter :: iounit = 10

      ! Just write the diagonal components!

      open(unit=iounit,file=filename,status="unknown")
    !  write(iounit,*) '# wn: orbital runs first, then spin'
      do w=1,size(func(1,1,1,:))
         write(iounit,'((ES23.16),(4X))',advance='no') wn(w-1)
         do s=1,nspin
            do m=1,size(func(:,1,s,w))
               write(iounit,'(2(ES23.16,2X),(4X))',advance='no')  & 
                                 & real(func(m,m,s,w)), aimag(func(m,m,s,w))
            enddo
         enddo
         write(iounit,'(A1)') ' '
      enddo
      close(iounit)
   end subroutine write_loc_matsfunc

!  ============================================================
!  == Write a local matrix-valued Green's function into file
!  ============================================================
   subroutine write_loc_matsfunc_matrix(filename,func)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:)  ! norb,norb,nspin,nomega
      integer(ki) :: w,s,m1,m2
      integer(ki),parameter :: iounit = 10

      open(unit=iounit,file=filename,status="unknown")
    !  write(iounit,*) '# wn: orbital runs first, then spin'
      do w=1,size(func(1,1,1,:))
         write(iounit,'((ES23.16),(4X))',advance='no') wn(w-1)
         do s=1,nspin
            do m1=1,size(func(:,1,s,w))
               do m2=1,size(func(1,:,s,w))
                  write(iounit,'(2(ES23.16,2X),(4X))',advance='no')  & 
                                 & real(func(m1,m2,s,w)), aimag(func(m1,m2,s,w))
                  !write(iounit,'(1(ES23.16,2X),(4X))',advance='no')  & 
                  !               & abs(func(m1,m2,s,w))
               enddo
            enddo
         enddo
         write(iounit,'(A1)') ' '
      enddo
      close(iounit)
   end subroutine write_loc_matsfunc_matrix

!  ============================================================
!  == Write a local bosonic matrix-valued Green's function into file
!  ============================================================
   subroutine write_loc_matsfunc_bos_matrix(filename,func)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:)  ! norb**2,norb**2,nspin,nnu
      integer(ki) :: w,s,m1,m2
      integer(ki),parameter :: iounit = 10

      open(unit=iounit,file=filename,status="unknown")
    !  write(iounit,*) '# wn: orbital runs first, then spin'
      do w=1,size(func(1,1,1,:))
         write(iounit,'((ES23.16),(4X))',advance='no') vn(w-1)
         do s=1,nspin
            do m1=1,size(func(:,1,s,w))
               do m2=1,size(func(1,:,s,w))
                  write(iounit,'(2(ES23.16,2X),(4X))',advance='no')  & 
                                 & real(func(m1,m2,s,w)), aimag(func(m1,m2,s,w))
                  !write(iounit,'(1(ES23.16,2X),(4X))',advance='no')  & 
                  !               & abs(func(m1,m2,s,w))
               enddo
            enddo
         enddo
         write(iounit,'(A1)') ' '
      enddo
      close(iounit)
   end subroutine write_loc_matsfunc_bos_matrix

!  ============================================================
!  == Write the local part of a nonlocal fermionic matrix-valued Green's function into file
!  ============================================================
   subroutine write_locpart_matsfunc(filename,func)
      use params
      use matsfunc_ops      
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:,:)  ! norb,norb,nspin,nkpts,nomega
      integer(ki),parameter :: iounit = 10

      complex(kr),allocatable  :: floc(:,:,:,:)  ! norb,norb,nspin,noemga
      allocate( floc(norb,norb,nspin,size(func(1,1,1,1,:))) )

      call get_locpart_mats(floc,func)
      call write_loc_matsfunc_matrix(filename,floc)

      deallocate(floc)
   end subroutine write_locpart_matsfunc


!  ============================================================
!  == Write the local part of a nonlocal bosonic matrix-valued Green's function into file
!  ============================================================
   subroutine write_locpart_matsfunc_bos(filename,func)
      use params
      use matsfunc_ops      
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: func(:,:,:,:,:)  ! norb**2,norb**2,nspin,nkpts,nnu
      integer(ki),parameter :: iounit = 10

      complex(kr),allocatable  :: floc(:,:,:,:)  ! norb**2,norb**2,nspin,nnu
      allocate( floc(norb**2,norb**2,nspin,size(func(1,1,1,1,:))) )

      call get_locpart_mats(floc,func)
      call write_loc_matsfunc_bos_matrix(filename,floc)

      deallocate(floc)
   end subroutine write_locpart_matsfunc_bos

!  ============================================================
!  == Write the nearest-neighbour Greensfunction and Selfenergy in x-direction
!  ============================================================
   subroutine write_gfsig_nonloc(x,y,z,filename,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT)
      use matsfunc_ops      
      integer(ki),intent(in)    :: x,y,z               ! which nonloc ?
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)    :: h_dft(:,:,:,:)      ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: v_xc(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: s_gw_dc(:,:,:,:)    ! norb,norb,nspin,nomega
      complex(kr),intent(in)    :: simp(:,:,:,:)       ! norb,norb,nspin,nomega
      real(kr),intent(in)       :: dft_hartree(:,:,:)  ! norb,norb,nspin
      real(kr),intent(in)       :: dft_exchange(:,:,:) ! norb,norb,nspin
      logical,intent(in)        :: onlyDFT             ! calculate only [iwn + mu - H_DFT]^-1

      integer(ki)             :: w,s,ikx,iky,ikz,ik
      complex(kr),allocatable :: gfnn(:,:,:,:)       ! norb,norb,nspin,nomega
      complex(kr),allocatable :: signn(:,:,:,:)      ! norb,norb,nspin,nomega

      allocate( gfnn(norb,norb,nspin,nomega) )
      allocate( signn(norb,norb,nspin,nomega) )

!  ====================================================
!  == Do the summation over all k ==
      gfnn = (0.0_kr,0.0_kr)
      signn = (0.0_kr,0.0_kr)

      do w=1,nomega
         do ikx=0,nkx-1
            do iky=0,nky-1
               do ikz=0,nkz-1
                  ik = ikx*nky*nkz + iky*nkz + ikz + 1
                  do s=1,nspin
                    gfnn(:,:,s,w) = gfnn(:,:,s,w)  &
               &  + gfk(s,ik,w,h_dft,v_xc,s_gw,s_gw_dc,simp,dft_hartree,dft_exchange,onlyDFT) &
               &       *exp(ci*ikx*x*2*pi/nkx)                              &
               &       *exp(ci*iky*y*2*pi/nky)                            &
               &       *exp(ci*ikz*z*2*pi/nkz)  /nkpts


               if ( useSigK /= 0 )then
                  signn(:,:,s,w) = signn(:,:,s,w)  &
                         &         + (simp(:,:,s,w)+s_gw(:,:,s,ik,w)-s_gw_dc(:,:,s,w))  &
                         &           *exp(ci*ikx*x*2*pi/nkx)                            &
                         &           *exp(ci*iky*y*2*pi/nky)                            &
                         &           *exp(ci*ikz*z*2*pi/nkz)  /nkpts
         
               endif

                  enddo ! s
               enddo ! kz
            enddo ! ky
         enddo ! kx
      enddo ! w
   
      call write_loc_matsfunc_matrix("g_"//filename//".dat",gfnn)
      call write_loc_matsfunc_matrix("s_"//filename//".dat",signn)

      deallocate( gfnn )
      deallocate( signn )
   end subroutine write_gfsig_nonloc

!  ============================================================
!  == Write combinations of local Selfenergies into files
!  ============================================================
   subroutine write_loc_sigmas(s_gw_loc,s_gw_dc,simp)
      use params
      use matsfunc_ops
      complex(kr),intent(in)       :: s_gw_loc(:,:,:,:)  ! norb,norb,nspin,nomega
      complex(kr),intent(in)       :: s_gw_dc(:,:,:,:)  ! norb,norb,nspin,nomega
      complex(kr),intent(in)       :: simp(:,:,:,:)  ! norb,norb,nspin,nomega
      integer(ki) :: w,s,m
      integer(ki),parameter   :: iounit1 = 10
      integer(ki),parameter   :: iounit2 = 11
      complex(kr)             :: tmp

      ! Just write the diagonal components!

      open(unit=iounit1,file="s_gw_loc_minus_dc.dat",status="unknown")
      open(unit=iounit2,file="s_loc.dat",status="unknown")
      do w=1,size(simp(1,1,1,:))
         write(iounit1,'((ES23.16),(4X))',advance='no') wn(w-1)
         write(iounit2,'((ES23.16),(4X))',advance='no') wn(w-1)
         do s=1,nspin
            do m=1,size(simp(:,1,1,1))
               tmp = s_gw_loc(m,m,s,w) - s_gw_dc(m,m,s,w)
               write(iounit1,'(2(ES23.16,2X),(4X))',advance='no') real(tmp),aimag(tmp)

               tmp = tmp + simp(m,m,s,w)
               write(iounit2,'(2(ES23.16,2X),(4X))',advance='no') real(tmp),aimag(tmp)
            enddo
         enddo
         write(iounit1,'(A1)') ' '
         write(iounit2,'(A1)') ' '
      enddo
      close(iounit1)
      close(iounit2)
   end subroutine write_loc_sigmas


!  ============================================================
!  == Print the Hartree and Exchange terms
!  ============================================================
   subroutine print_dft_hartree_exchange(hartree,exchange)
      use params
      real(kr),intent(in) :: hartree(:,:,:)   ! norb,norb,nspin
      real(kr),intent(in) :: exchange(:,:,:)  ! norb,norb,nspin
      integer :: s,m1,m2

      write(*,'(A)') 'Hartree term:'
      do s=1,nspin
         write(*,'(A,I2)') ' spin=',s
         do m1=1,norb
            do m2=1,norb
               write(*,'((ES10.3),(2X))',advance='no') hartree(m1,m2,s)
            enddo
            write(*,'(A1)') ' '
         enddo
         write(*,'(A1)') ' '
      enddo

      write(*,'(A)') 'Exchange term:'
      do s=1,nspin
         write(*,'(A,I2)') ' spin=',s
         do m1=1,norb
            do m2=1,norb
               write(*,'((ES10.3),(2X))',advance='no') exchange(m1,m2,s)
            enddo
            write(*,'(A1)') ' '
         enddo
        write(*,'(A1)') ' '
      enddo
   end subroutine print_dft_hartree_exchange
   
!  ============================================================
!  == Print the local V_xc terms
!  ============================================================
   subroutine print_vxc_local(v_xc)
      use params
      complex(kr),intent(in) :: v_xc(:,:,:,:)      ! norb,norb,nspin,nkpts
      integer :: s,m1,m2,k
      complex(kr) :: tmp

      if (useSigK==1) then
          write(*,'(A)') 'Local V_XC:'
          do s=1,nspin
             write(*,'(A,I2)') ' spin=',s
             do m1=1,norb
                do m2=1,norb
                   ! sum over all k
                   tmp = 0.0_kr
                   do k=1,nkpts
                      tmp = tmp + v_xc(m1,m2,s,k)/nkpts
                   enddo

                   write(*,'( A,2(F7.3),A,(2X))',advance='no') "(",tmp,")"
                enddo
                write(*,'(A1)') ' '
             enddo
             write(*,'(A1)') ' '
          enddo
       endif
   end subroutine print_vxc_local
   
!  ============================================================
!  == Write a k-dependent but static object into file
!  ============================================================
   subroutine write_kdep_static_f(filename,f,diag)   
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: f(:,:,:,:)       ! norb,norb,nspin,nkpts
      integer                      :: diag             ! should we diagonalize f?
      integer :: iounit=10, s,m,m2,k,ik, ndim
      real(kr),allocatable    :: eigvals(:)
      complex(kr),allocatable :: eigvectors(:,:)

      ndim = size(f(:,1,1,1))
      allocate(eigvals(ndim))
      allocate(eigvectors(ndim,ndim))
   
      open(unit=iounit,file=filename,status="unknown")

!      do k=1,nkpts
!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) 'REMOVE MODIFICATION IN write_kdep_static_f !!!!!!!!!!!!!!!!!'
do ik=0,nkx/2 + nky/2 + nky/2 + nkz/2
    ! This is really bad stuff, change this
    k = ik*nkz*nky +1

    ! X to M
    if (ik .ge. nkx/2) k = (nkx/2)*nkz*nky  + (ik-nkx/2)*nkz +1

    ! M to Gamma
    if (ik .ge. nkx/2 + nky/2) k = (nkx/2)*nkz*nky + (nky/2)*nkz - (ik-nkx/2-nky/2)*( nkz*nky + nkz ) +1

    ! Gamma to R
   ! if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = (ik - (nkx/2 + nky/2 + nky/2))* ( nkz*nky + nkz + 1 )  +1

    ! Gamma to Z
    if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = (ik - (nkx/2 + nky/2 + nky/2))  +1


    !write(*,*) ik,k
!!!!!!!!!!!!!!!!!!!!!!!!!
!         write(iounit,'(3(F10.5,2X))',advance='no') kvecs(k,1),kvecs(k,2),kvecs(k,3)
      
         do s=1,1 !nspin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! diagonalize tight-binding 
         if (diag==1) then
            call get_eigenvectors(eigvals, eigvectors, f(:,:,s,k) , ndim )
            do m=1,ndim
               write(iounit,'(I3,2XF10.5,2X)',advance='no') ik,eigvals(m)
               do m2=1,ndim
                  write(iounit,'(1(F10.5,2X))',advance='no') abs( eigvectors(m2,m) )
               enddo
               write(iounit,'(A)') ' '
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NORMAL
         else
            do m=1,ndim
               write(iounit,'(2(F10.5,2X))',advance='no') real(f(m,m,s,k)),aimag(f(m,m,s,k))
            enddo !m loop
            write(iounit,'(A1)') ' '
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo !s loop
!         write(iounit,'(A1)') ' '
      enddo !k loop
      close(iounit)

      deallocate(eigvals)
      deallocate(eigvectors)
   end subroutine write_kdep_static_f

!  ============================================================! 
!  == Write a k-dependent object on 2D fermisurface at w=0
!  ============================================================
   subroutine write_kdep_fs(filename,f,fit)   
     use params
      use matsfunc_ops
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: f(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega
      integer(ki),intent(in)       :: fit          

      integer     :: iounit=10, s,m,m2,ikx,iky,ik
      complex(kr) :: coeff(5)

      if (size(f(1,1,1,:,1))==nkpts) then

      open(unit=iounit,file=filename,status="unknown")

      do ikx=1,nkx
         do iky=1,nky
            ik = (ikx-1)*nky*nkz + (iky-1)*nkz + 1

            write(iounit,'(2(F10.5,2X))',advance='no') (ikx-1)*2*pi/nkx, (iky-1)*2*pi/nky
            do s=1,nspin
               do m=1,size(f(:,1,1,1,1))
                  if (fit==1) then
                     coeff = get_lowfreq_coeff(f(m,m,s,ik,:),5,6)
                     write(iounit,'(2(F10.5,2X))',advance='no') real(coeff(1)),aimag(coeff(1))
                  else
                     write(iounit,'(2(F10.5,2X))',advance='no') real(f(m,m,s,ik,1)),aimag(f(m,m,s,ik,1))
                  endif
               enddo !m loop
            enddo ! s loop
            write(iounit,'(A1)') ' '
         enddo ! iky
         write(iounit,'(A1)') ' '
      enddo ! ikx
      close(iounit)

      endif

   end subroutine write_kdep_fs
   
!  ============================================================! 
!  == Write a k-dependent object into file interpolated to w=0
!  ============================================================
   subroutine write_kdep_w0_f(filename,f)   
     use params
      use matsfunc_ops
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       :: f(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega

      integer     :: iounit=10, s,m,m2,k,ik
      complex(kr) :: coeff(5)
   
      open(unit=iounit,file=filename,status="unknown")

do ik=0,nkx/2 + nky/2 + nky/2 + nkz/2
    ! This is really bad stuff, change this
    k = ik*nkz*nky +1

    ! X to M
    if (ik .ge. nkx/2) k = (nkx/2)*nkz*nky  + (ik-nkx/2)*nkz +1

    ! M to Gamma
    if (ik .ge. nkx/2 + nky/2) k = (nkx/2)*nkz*nky + (nky/2)*nkz - (ik-nkx/2-nky/2)*( nkz*nky + nkz ) +1

    ! Gamma to R
   ! if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = (ik - (nkx/2 + nky/2 + nky/2))* ( nkz*nky + nkz + 1 )  +1

    ! Gamma to Z
    if (ik .ge. nkx/2 + nky/2 + nky/2 ) k = (ik - (nkx/2 + nky/2 + nky/2))  +1

         do s=1,1 !nspin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do m=1,size(f(:,1,s,k,1))
               coeff = get_lowfreq_coeff(f(m,m,s,k,:),5,6)
               write(iounit,'(2(F10.5,2X))',advance='no') real(coeff(1)),aimag(coeff(1))
               !write(iounit,'(2(F10.5,2X))',advance='no') real(f(m,m,s,k,1)),aimag(f(m,m,s,k,1))
            enddo !m loop
            write(iounit,'(A1)') ' '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo !s loop
!         write(iounit,'(A1)') ' '
      enddo !k loop
      close(iounit)

   end subroutine write_kdep_w0_f

!  ============================================================
!  == Write Hartree and Exchange terms into a file
!  ============================================================
   subroutine write_hartree_exchange(filename,hartree,exchange)
      use params
      character(len=*), intent(in) :: filename
      real(kr),intent(in)          :: hartree(:,:,:)   ! norb,norb,nspin
      real(kr),intent(in)          :: exchange(:,:,:)  ! norb,norb,nspin
      integer :: iounit=10, s,m1,m2
      
      ! write only the diagonal terms, since hartree and exchange are diagonal

      ! first Hartree
      open(unit=iounit,file=filename,status="unknown")
      do s=1,nspin
         !write(iounit,'(A5,I2)') 'Spin=',s
         do m1=1,norb
            do m2=1,norb
               write(iounit,'(ES23.16,2X)',advance='no') hartree(m1,m2,s)
            enddo
            write(iounit,'(A1)') ' '
         enddo
         write(iounit,'(A1)') ' '
      enddo

      ! add empty line
      write(iounit,'(A1)') ' '

      ! Then exchange
      do s=1,nspin
         !write(iounit,'(A5,I2)') 'Spin=',s
         do m1=1,norb
            do m2=1,norb
               write(iounit,'(ES23.16,2X)',advance='no') exchange(m1,m2,s)
            enddo
            write(iounit,'(A1)') ' '
         enddo
         write(iounit,'(A1)') ' '
      enddo
      close(iounit)
   end subroutine write_hartree_exchange
   
!  ============================================================
!  == Write the chemical potential into a file
!  ============================================================
   subroutine write_mu(filename)
      use params
      character(len=*), intent(in) :: filename
      integer :: iounit=10
      
      open(unit=iounit,file=filename,status="unknown")
      write(iounit,'(ES23.16)') mu
      
      close(iounit)
   end subroutine write_mu

!  ============================================================
!  == Print the filling of a local Green's function
!  ============================================================
   subroutine print_local_filling(myname,gf)
      use params
      use matsfunc_ops
      character(len=*), intent(in) :: myname
      complex(kr),intent(in)       :: gf(:,:,:,:) ! norb,norb,nspin,nomega
      integer                      :: s,m

      write(*,'(A,A,A)') 'Filling of ',myname,':'
      do s=1,nspin
         write(*,'((A6),(I1),A1,4X)',advance='no') ' Spin=',s,':'
         do m=1,norb
            write(*,'((F6.3),2X)',advance='no') filling(gf(m,m,s,:))
         enddo
         write(*,'(A1)') ' '
      enddo
   end subroutine print_local_filling

!  ============================================================
!  == Print the local orbital levels
!  ============================================================
   subroutine print_local_orbital_levels(mu_loc)
      use params
      real(kr),intent(in)  :: mu_loc(:,:,:) ! norb,norb,nspin
      integer              :: s,m1,m2

      write(*,'(A)') 'Local orbital levels:'
      do s=1,nspin
         write(*,'((A5),(I1),A1,4X)') 'Spin=',s,':'
         do m1=1,norb
            do m2=1,norb
               write(*,'((F7.3),2X)',advance='no') mu_loc(m1,m2,s)
            enddo
            write(*,'(A1)') ' '
         enddo
         !write(*,'(A1)') ' '
      enddo
      write(*,'(A)') ' '
   end subroutine print_local_orbital_levels

!  ============================================================
!  == Check offdiagonal elements of local matsubara function of size norb x norb
!  ============================================================
   subroutine check_offdiagonal(myname,f)
      use params
      character(len=*), intent(in) :: myname
      complex(kr),intent(in) :: f(:,:,:,:)  ! norb,norb,nspin,nomega
      integer                :: w,s,m1,m2,a1,a2
      integer                :: mw,ms,mm1,mm2      ! indices for overall maximum
!      integer                :: amw,ams,amm1,amm2  ! indices for atomic submatrices, i.e. interatomic terms neglected
      real(kr)               :: maxitem,maxitem_atom
      maxitem=0.0_kr
      mw = -1
      ms = -1
      mm1 = -1
      mm2 = -1

!      maxitem_atom=0.0_kr
!      amw = -1
!      ams = -1
!      amm1 = -1
!      amm2 = -1

      do w=1,nomega
         do s=1,nspin
            do m1=1,norb
               do m2=1,norb
                  ! check the overall maximum
                  if (m1/=m2 .and. maxitem<abs(f(m1,m2,s,w)) ) then
                     maxitem = abs(f(m1,m2,s,w))
                     mw = w
                     ms = s
                     mm1 = m1
                     mm2 = m2
                  endif

!                  ! check for offdiagonals on a single DMFT atom
!                  a1 = int( (m1-1) / (norb/noEquivAtoms) )
!                  a2 = int( (m2-1) / (norb/noEquivAtoms) )
!                  if (m1/=m2 .and. a1==a2 .and. maxitem_atom<abs(f(m1,m2,s,w)) ) then
!                     maxitem_atom = abs(f(m1,m2,s,w))
!                     amw = w
!                     ams = s
!                     amm1 = m1
!                     amm2 = m2
!                  endif
                
               enddo    
            enddo    
         enddo    
      enddo    

      ! Print info about general maximum      
      if ( maxitem>0.001_kr ) then
         write(*,'(A)') '============================================================='
         write(*,'(A,A,A)') 'WARNING: Offdiagonal elements in ',myname,' are not zero!!!'
         write(*,'(A17,I3,A4,I3,A5,I3,A5,I3,A,2F7.3,A)') 'Max.element at n=', &
                        &     mw,', s=',ms,', m1=',mm1,', m2=',mm2,': (',f(mm1,mm2,ms,mw), ')'
         write(*,'(A)') '============================================================='
      endif

!      ! Print info maximum on single atom
!      if ( maxitem_atom>0.001_kr ) then
!         write(*,'(A)') '============================================================================='
!         write(*,'(A)') 'WARNING: Offdiagonal elements in ',myname,' ON A SINGLE ATOM are not zero!!!'
!         write(*,'(A21,I3,A6,I3,A4,I3,A5,I3,A5,I3,A,2F7.3,A)') 'Max.element for atom=',               &
!                        & int((amm1-1)/(norb/noEquivAtoms)),' at n=', &
!                        &     amw,', s=',ams,', m1=',amm1,', m2=',amm2,': (',f(amm1,amm2,ams,amw), ')'
!         write(*,'(A)') 'WARNING: THIS USUALLY MEANS SOMETHING IS GOING REALLY WRONG !!!'
!         write(*,'(A)') 'f(:,:,s,n):'
!         do m1=1,norb
!            do m2=1,norb
!               write(*,'(ES10.3,1X,ES10.3,3X)',advance='no') real(f(m1,m2,ams,amw)),aimag(f(m1,m2,ams,amw))
!            enddo
!            write(*,'(A1)') ' '
!         enddo
!         write(*,'(A)') '============================================================================='
!      endif

   end subroutine check_offdiagonal

!  ============================================================
!  == Check offdiagonal elements of local matsubara function on
!  == a single DMFT atom of size norbPerAtom x norbPerAtom
!  ============================================================
   subroutine check_offdiagonal_dmftatom(myname,f)
      use params
      character(len=*), intent(in) :: myname
      complex(kr),intent(in) :: f(:,:,:,:)  ! norb,norb,nspin,nomega
      integer                :: w,s,m1,m2,a1,a2
      integer                :: mw,ms,mm1,mm2      ! indices for overall maximum
      real(kr)               :: maxitem
      maxitem=0.0_kr
      mw = -1
      ms = -1
      mm1 = -1
      mm2 = -1

      do w=1,nomega
         do s=1,nspin
            do m1=1,norbPerAtom
               do m2=1,norbPerAtom
                  ! check the maximum
                  if (m1/=m2 .and. maxitem<abs(f( dmftOrbsIndex(m1), dmftOrbsIndex(m2), s,w)) ) then
                     maxitem = abs(f( dmftOrbsIndex(m1), dmftOrbsIndex(m2), s,w))
                     mw = w
                     ms = s
                     mm1 = m1
                     mm2 = m2
                  endif
               enddo
            enddo
         enddo
      enddo

      ! Print info about general maximum
      if ( maxitem>0.001_kr ) then
         write(*,'(A)') '============================================================='
         write(*,'(A,A,A)') 'WARNING: Offdiagonal elements in ',myname,' ON THE DMFT ATOM are not zero!!!'
         write(*,'(A17,I3,A4,I3,A5,I3,A5,I3,A,2F7.3,A)') 'Max.element at n=', &
                        &     mw,', s=',ms,', m1=',mm1,', m2=',mm2,': (',f(dmftOrbsIndex(mm1),dmftOrbsIndex(mm2),ms,mw), ')'
         write(*,'(A,I2,A,I2,A)') '(Between the global orbital ',dmftOrbsIndex(mm1),      &
                               &  ' and ',dmftOrbsIndex(mm2),' )'
         write(*,'(A)') '============================================================='
      endif

   end subroutine check_offdiagonal_dmftatom


!  ============================================================
!  == Print a general matrix to screen
!  ============================================================
   subroutine print_matrix(matrix,N)
      use params
      complex(kr),intent(in) :: matrix(:,:)     ! N,N
      integer(ki),intent(in) :: N
      integer(ki)         :: m1,m2
      write(*,'(A)') 'Real Part:'
      do m1=1,N
         do m2=1,N
            write(*,'((F7.4),(4X))',advance='no') real(matrix(m1,m2))
         enddo
         write(*,'(A1)') ' '
      enddo
!      write(*,'(A)') ' '

      write(*,'(A)') 'Imaginary Part:'
      do m1=1,N
         do m2=1,N
            write(*,'((F7.4),(4X))',advance='no') aimag(matrix(m1,m2))
         enddo
         write(*,'(A1)') ' '
      enddo
      write(*,'(A)') ' '

   end subroutine print_matrix

!  ============================================================
!  == Write the Uijkl matrix for the solver
!  ============================================================
   subroutine write_Uijkl_matrix_solver(Uijkl,filename)
      use params
      complex(kr),intent(in) :: Uijkl(:,:,:,:)        ! norbPerAtom**2,norbPerAtom**2,nspin,nspin
      character(len=*), intent(in) :: filename

      integer(ki),parameter :: iounit = 10
      integer(ki)           :: s1,s2,i,j,k,l, nUelem, n, m1,m2
      real(kr)              :: smallestUval = 0.00001

      ! First count nonzero elements
      nUelem = 0
      do i=1,norbPerAtom**2
         do j=1,norbPerAtom**2
            do s1=1,nspin
               do s2=1,nspin

                  if ( abs(Uijkl(i,j,s1,s2)) > smallestUval ) nUelem=nUelem+1

               enddo ! s2
            enddo ! s1
         enddo ! j
      enddo ! i
      write(*,'(A,I4,A)') 'Found ',nUelem,' nonzero Uijkl elements!'

      open(unit=iounit,file=filename,status="unknown")
      write(iounit,'(I4)') nUelem

      n = 0
      do i=0,norbPerAtom-1
         do j=0,norbPerAtom-1
            do k=0,norbPerAtom-1
               do l=0,norbPerAtom-1
                  do s1=0,nspin-1
                     do s2=0,nspin-1

                        m1 = i*norbPerAtom+j
                        m2 = l*norbPerAtom+k   ! different order for ALPS Cthyb: l<->k

                        if ( abs(Uijkl(m1+1,m2+1,s1+1,s2+1)) > smallestUval ) then
                           write(iounit,'(5(I4,2X),2(F14.10,3X))') n, 2*i+s1, 2*j+s2, 2*k+s2, 2*l+s1,  &
                                      &  real( Uijkl(m1+1,m2+1,s1+1,s2+1) ) ,                             &
                                      & aimag( Uijkl(m1+1,m2+1,s1+1,s2+1) )                               
             
                           n = n+1
                        endif
  
                     enddo ! s2
                  enddo ! s1
               enddo ! l
            enddo ! k
         enddo ! j
      enddo ! i

      close(iounit)

   end subroutine write_Uijkl_matrix_solver

!  ============================================================
!  == Print the U-matrix to screen
!  ============================================================
   subroutine print_umatrix(uloc)
      use params
      use matsfunc_ops
      complex(kr),intent(in) :: uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu
      integer(ki)            :: m1,m2,s1,s2
      real(kr)               :: coeffs(2)

      write(*,'(A)') 'Bare interaction matrix elements w=infty ( U_(ik)(jl) notation):'
      do m1=1,size(uloc(:,1,1,1))
         do m2=1,size(uloc(1,:,1,1))
            if ( size(uloc(1,1,1,:)) > nfit ) then
               coeffs = get_highfreq_coeff_bosonic(uloc(m1,m2,1,:))
               write(*,'((F7.4),(4X))',advance='no') coeffs(1)
            else
               write(*,'((F7.4),(4X))',advance='no') real(uloc(m1,m2,1,1)) 
            endif
         enddo
         write(*,'(A1)') ' '
      enddo
      write(*,'(A)') ' '

      write(*,'(A)') 'Screened interaction matrix elements w=0 ( U_(ik)(jl) notation):'
      do m1=1,size(uloc(:,1,1,1))
         do m2=1,size(uloc(1,:,1,1))
            write(*,'((F7.4),(4X))',advance='no') real(uloc(m1,m2,1,1))
         enddo
         write(*,'(A1)') ' '
      enddo
      write(*,'(A)') ' '


   end subroutine print_umatrix
   
!  ============================================================
!  == Print the high-freq. coefficients to screen
!  ============================================================
   subroutine print_highfreq_coeffs(myname,f)
      use params
      use matsfunc_ops
      character(len=*),intent(in) :: myname
      complex(kr),intent(in)      :: f(:,:,:,:)     ! norb,norb,nspin,nomega
      real(kr)                    :: coeffs(5)
      integer                     :: s,m,i

      write(*,'(A,A,A)') ' High-frequency coeffs. for ',myname,':'
      do m=1,size(f(:,1,1,1))
         do s=1,nspin
            coeffs = get_highfreq_coeff(f(m,m,s,:),0)
            write(*,'(A5,I2,A7,I1,A2)',advance='no') ' orb=',m,', spin=',s,': '
            do i=1,5
               write(*,'(A1,I1,A2,ES10.3,3X)',advance='no') 'c',i-1,'= ',coeffs(i)
            enddo
            write(*,'(A1)') ' '
         enddo
      enddo
      
   end subroutine print_highfreq_coeffs

!  ============================================================
!  == Print the maximal difference between two local Matsubara functions at a
!  frequency
!  ============================================================
   subroutine print_matsf_diff_atw(f1,f2,n,name1,name2)
      use params
      !use matsfunc_ops
      complex(kr),intent(in)      :: f1(:,:,:,:),f2(:,:,:,:) ! norb,norb,nspin,nomega
      integer(ki),intent(in)      :: n                       ! at which omega
      character(len=*),intent(in) :: name1,name2
      integer                     :: w,s,m1,m2
      complex(kr)  :: diff,fval
      
      diff = 0.0_kr
      fval = 1.0_kr
      
      !do w=1,nomega
         w = n
         do s=1,nspin

            do m1=1,norbPerAtom
              m2 = dmftOrbsIndex(m1)
              if ( abs(f1(m2,m2,s,w)-f2(m2,m2,s,w)) > abs(diff) ) then
                 diff =   f2(m2,m2,s,w) -f1(m2,m2,s,w)  
                  
                 fval = f2(m2,m2,s,w)
              endif
            enddo !m1

         enddo !s
      !enddo !w
      write(*,'(A,A,A,A,I4,A,SP,F8.4,A,F6.1,A,F8.4,A,F6.1,A)') name1,' <-> ',name2,' diff at Freq n= ',n,':',  &
            & real(diff),' (',100*real(diff)/real(fval),'%), ',aimag(diff),'i (',100*aimag(diff)/aimag(fval),'%), '
   
   end subroutine print_matsf_diff_atw
   
!  ============================================================
!  == Print the maximal difference between two local Matsubara functions
!  ============================================================
   subroutine print_matsf_diff(f1,f2,name1,name2)
      use params
      !use matsfunc_ops
      complex(kr),intent(in)      :: f1(:,:,:,:),f2(:,:,:,:) ! norb,norb,nspin,nomega
      character(len=*),intent(in) :: name1,name2
      integer                     :: w,s,m1,m2
      real(kr)  :: diff, diffDiag, diffDiagAtom
      
      diff = 0.0_kr
      diffDiag = 0.0_kr
      diffDiagAtom = 0.0_kr
      do w=1,size(f1(1,1,1,:))
         do s=1,nspin
            do m1=1,size(f1(:,1,1,1))
               do m2=1,size(f1(:,1,1,1))

                  if (abs(f1(m1,m2,s,w)-f2(m1,m2,s,w))>diff) then
                     diff = abs(f1(m1,m2,s,w)-f2(m1,m2,s,w))
                  endif

                  if (m1==m2 .and. abs(f1(m1,m2,s,w)-f2(m1,m2,s,w))>diffDiag) then
                     diffDiag = abs(f1(m1,m2,s,w)-f2(m1,m2,s,w))
                  endif

               enddo !m2
            enddo !m1

            do m1=1,norbPerAtom
              if (abs(f1(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s,w)-f2(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s,w))>diffDiagAtom) then
                 diffDiagAtom = abs(  f1(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s,w)          &
                            &        -f2(dmftOrbsIndex(m1),dmftOrbsIndex(m1),s,w)  )
              endif
            enddo !m1

         enddo !s
      enddo !w
      write(*,'(A,A,A,A,F8.4,A,F8.4,A,F8.4)') name1,' <-> ',name2,' diagAtomDiff: ',diffDiagAtom,      &
                                               &                  ', diagDiff: ',   diffDiag,          &
                                               &                  ', totDiff: ',    diff
   
   end subroutine print_matsf_diff


!  ============================================================
!  == Print the bare interaction parameters U,J
!  ============================================================
   subroutine print_bare_interaction()
      use params
      !use matsfunc_ops

      write(*,'(A)') 'The Slater-averaged bare interaction values:'
      write(*,'(A9,F8.3)') ' U_bare = ',Uinput
      write(*,'(A9,F8.3)') ' J_bare = ',Jinput

      if (usePolK==1) then
         write(*,'(A)') ' (calculated from V_q input data)'
      else
         write(*,'(A)') ' (as taken from the DMFT input file)'
      endif

   end subroutine print_bare_interaction

!  ============================================================
!  == Read data from previous calculation
!  ============================================================
   subroutine read_old_data(simp,pimp,s_gw_dc,p_gw_dc,dft_hartree,dft_exchange)
      use params
      complex(kr),intent(inout) :: simp(:,:,:,:)         ! norb,norb,nspin,nomega
      complex(kr),intent(inout) :: pimp(:,:,:,:)         ! norb**2,norb**2,nspin,nnu
      complex(kr),intent(inout) :: s_gw_dc(:,:,:,:)      ! norb,norb,nspin,nomega
      complex(kr),intent(inout) :: p_gw_dc(:,:,:,:)      ! norb**2,norb**2,nspin,nnu
      real(kr),intent(inout)    :: dft_hartree(:,:,:)    ! norb,norb,nspin
      real(kr),intent(inout)    :: dft_exchange(:,:,:)   ! norb,norb,nspin
      real(kr) :: tmp(norb)
      integer  :: iounit = 10,s,m1,m2

      simp         = (0.0_kr, 0.0_kr)
      s_gw_dc      = (0.0_kr, 0.0_kr)
      pimp         = (0.0_kr)
      p_gw_dc      = (0.0_kr)
      dft_hartree  = (0.0_kr)
      dft_exchange = (0.0_kr)
      mu = 0.0_kr
 
      write(*,'(A)') '===Read input from old calculation:==='
      write(*,'(A)') ' Read impurity Selfenergy from s_imp.dat...'
      call read_loc_matsfunc_matrix("s_imp.dat",simp)

      if ( useSigK/=1 ) then
          write(*,'(A)') ' Read Selfenergy doublecounting from s_gw_dc.dat...'
          call read_loc_matsfunc_matrix("s_gw_dc.dat",s_gw_dc)
      endif

      if ( usePolK/=0 ) then
         write(*,'(A)') ' Read impurity Polarization from p_imp.dat...'
         call read_loc_matsfunc_matrix("p_imp.dat",pimp)
         write(*,'(A)') ' Read Polarization doublecounting from p_gw_dc.dat...'
         call read_loc_matsfunc_matrix("p_gw_dc.dat",p_gw_dc)
      endif

      write(*,'(A)') ' Read DFT Hartree and Exchange from dft_hartree_exchange.dat...'
      call read_hartree_exchange("dft_hartree_exchange.dat",dft_hartree,dft_exchange)

      write(*,'(A)') ' Read chemical potential from mu.dat...'
      open(unit=iounit,file="mu.dat",status="old")
      read(iounit,*) mu
      close(iounit)
   end subroutine read_old_data

!  ============================================================
!  == Read data from previous DCA calculation
!  ============================================================
   subroutine read_old_dca(s_gw,s_gw_dc)
      use params
      use fourier_interpolation
      complex(kr),intent(inout)    :: s_gw(:,:,:,:,:)     ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(inout)    :: s_gw_dc(:,:,:,:)     ! norb,norb,nspin,nomega

      real(kr) :: tmp(norb)
      complex(kr),allocatable :: simpK(:,:,:)          ! npx*npy,nspin,nomega
      complex(kr),allocatable :: s_gw_coarse(:,:,:,:,:)     ! norb,norb,nspin,npx*npy,nomega
      integer  :: iounit = 10,s,m,k,ikx,iky,ip,i, kdim

      if (npx*npy .gt. 0) then
         allocate(simpK(npx*npy,nspin,nomega))
         allocate( s_gw_coarse(1,1,nspin,npx*npy,nomega ) )

         s_gw = (0.0_kr, 0.0_kr)
 
         write(*,'(A)') ' Read patched DCA Selfenergy from s_impK_DCA.dat...'
         call read_loc_matsfuncvec("s_impK_DCA.dat",simpK)

         do k=1,npx*npy
            s_gw_coarse(1,1,:,k,:)  = simpK(k,:,:)
         enddo
         call write_loc_matsfuncvec("s_impK_DCA_init.dat",simpK)
   
         ! copy selfenergy to s_gw, this is still patchy
         do ikx=0,nkx-1
            do iky=0,nky-1
               k = ikx*nky + iky +1 
               ip = kPatchIndex(k)
               s_gw(1,1,:,k,:)     = simpK(ip,:,:)
               ! norb,norb,nspin,nkpts,nomega
            enddo
         enddo      
         !call write_kdep_static_f("sgw_DCA_patch_init.dat",s_gw(:,:,:,:,1),0)
         call write_kdep_w0_f("sgw_DCA_patch_init.dat",s_gw)
         
         ! or interpolate the DCA selfenergy ??
         if (interpDCA==1) then
            if (npx .ne. npy) stop 'ERROR: Interpolation of DCA Selfenergy only possible for npx==npy !'
            do i=1,nomega
               call interpolate_kmatrix_cubic(s_gw(:,:,:,:,i), s_gw_coarse(:,:,:,:,i),nkx,nky,nkz,npx, npy,1)
            enddo
            call write_sgw_input(s_gw)
         endif
         ! this shouldnt be necessary since we read s_gw_dc from old calculation
!         call get_locpart_mats(s_gw_dc,s_gw) ! set DC to local part of DCA Selfenergy 

         !call write_kdep_static_f("sgw_DCA_interpolated_init.dat",s_gw(:,:,:,:,1),0)
         call write_kdep_w0_f("sgw_DCA_interpolated_init.dat",s_gw)

         deallocate(simpK)
         deallocate( s_gw_coarse )
      endif
   end subroutine read_old_dca

!  ============================================================
!  == Read the unitary transformation for the impurity model
!  == For now we assume trhe transformation to be real!
!  ============================================================
   subroutine read_utrafo()
      use params

      integer(ki) :: s,m1,m2
      integer(ki),parameter :: iounit = 10
      real(kr)              :: readindata(norbPerAtom) ! temparray for readin

      ! Utrafo(norbPerAtom,norbPerAtom,nspin)
      ! UtrafoT(norbPerAtom,norbPerAtom,nspin)

      ! Initialize with identity
      Utrafo = (0.0_kr,0.0_kr)
      UtrafoT = (0.0_kr,0.0_kr)
      do s=1,nspin
         do m1=1,norbPerAtom
            Utrafo(m1,m1,s)  = 1.0_kr
            UtrafoT(m1,m1,s) = 1.0_kr
         enddo
      enddo

      if (diagBath==1) then
         write(*,'(A)') "Reading transformation matrix from Utrafo.dat (assumed to be real...)"

         open(unit=iounit,file="Utrafo.dat",status="old")

         do s=1,nspin
            do m1=1,norbPerAtom
               read(iounit,*) readindata
               Utrafo(m1,:,s) = readindata
            enddo ! m1
 
            call get_inverse( UtrafoT(:,:,s),  Utrafo(:,:,s), norbPerAtom )

         enddo ! s

         close(iounit)
      endif

   end subroutine read_utrafo

!  ============================================================
!  == Read Hartree and exchange terms from a file
!  ============================================================
   subroutine read_hartree_exchange(filename,hartree,exchange)
      use params
      character(len=*), intent(in)    :: filename
      real(kr),intent(inout)          :: hartree(:,:,:)   ! norb,norb,nspin
      real(kr),intent(inout)          :: exchange(:,:,:)  ! norb,norb,nspin
      integer          :: iounit=14, s,m1,m2
      character(len=1)     :: tmp
      real(kr),allocatable :: readdata(:) ! norb

      hartree  = (0.0_kr)      
      exchange  = (0.0_kr)      
      allocate( readdata(norb) )

      ! first Hartree
      open(unit=iounit,file=filename,status="old")
      do s=1,nspin
         do m1=1,norb
            read(iounit,*) readdata
            do m2=1,norb
               hartree(m1,m2,s) = readdata( m2 )
            enddo
         enddo
      enddo

      ! skip empty line
     ! read(iounit,*) tmp
     ! This is not needed, FORTRAN somehow reads the first data point and skips
     ! empty lines... strange

      ! Then exchange
      do s=1,nspin
         do m1=1,norb
            read(iounit,*) readdata
            do m2=1,norb
               exchange(m1,m2,s) = readdata(m2)
            enddo
         enddo
      enddo
      close(iounit)

      deallocate( readdata )
   end subroutine read_hartree_exchange

!  ============================================================
!  == Read U(inu_nu) from a file, i.e. put it into Uavg(inu_n), Javg(inu_n)
!  ============================================================
   subroutine read_Uavg_Javg(Uavg,Javg,filename)
      use params
      real(kr),intent(inout)   :: Uavg(:), Javg(:)      ! nnu
      character(len=*),intent(in) :: filename
      real(kr) :: tmp(2)
      integer  :: iounit = 10,n

      Uavg = (0.0_kr)
      Javg = (0.0_kr)
      
      write(*,'(A)') ' Assume only F0 to be freq.dep. and set J constant to Jinput...'
      write(*,'(A,A,A)') 'Reading freq.dep. F0 values from file ',filename,' ...'
      open(unit=iounit,file=filename,status="old")
      do n=1,nnu
         read(iounit,*) tmp
         Uavg(n) = tmp(2)
         Javg(n) = Jinput
      enddo
      close(iounit)

      write(*,'(A,F9.4,A,F9.4,A,F9.4)') ' Reading: F0(0)=',Uavg(1),' ... F0(nnu/2)=', Uavg(1+nnu/2),' ... F0(nnu)=', Uavg(nnu)
      write(*,'(A,F9.4)') ' F0 bare (from input file) = ',Uinput

      if ( abs(Uavg(nnu)-Uinput)>1.0_kr ) then
        write(*,'(A)')              '!!! ============================================================================== !!!'
        write(*,'(A)')              '!!! WARNING: Uinput and U(inu_n->inf) are very different !'
        write(*,'(A,F9.4,A5,F9.4)') '!!! ',Uavg(nnu),' <-> ',Uinput
        write(*,'(A)')              '!!! Please check whether you forgot to set Uinput,Jinput correctly in the DMFT input file!'
        write(*,'(A)')              '!!! ============================================================================== !!!'
      endif
      write(*,'(A)') ' '

   end subroutine read_Uavg_Javg

!  ============================================================
!  == Obtain Uavg(inu_n), Javg(inu_n), depends useUw and usePolK
!  ============================================================
   subroutine obtain_Uavg_Javg(Uavg,Javg,uloc)
        use params
        use matsfunc_ops
        real(kr),intent(inout) :: Uavg(:)        ! nnu
        real(kr),intent(inout) :: Javg(:)        ! nnu
        complex(kr),intent(in)    :: uloc(:,:,:,:)  ! norb**2,norb**2,nspin,nnu

        Uavg = (0.0_kr)
        Javg = (0.0_kr)

        ! Do the full program
        if ( usePolK==1 .and. useUw==1 ) then
            stop "obtain_Uavg_Javg not implemented for usePolK==1 .and. useUw==1"
!            call convert_uloc_to_2particle_basis(Uavg,Javg,uloc,Prod2PartOverlap)

        ! Use a freq.dep. U but no selfconsistency on U
        else if ( usePolK/=1 .and. useUw==1 ) then
            call read_Uavg_Javg(Uavg,Javg,"Unninu2.dat")

        ! We have no freq.dep. U, so set Uavg,Javg equal to the values we read from the input file
        ! Uavg,Javg only have one element
        else
            if (nnu/=1) then
               write(*,'(A)') 'ERROR: nnu /= 1, but usePolK/=1 AND useUw==1, this should not happen !!!'
               stop
            endif

            Uavg(1) = Uinput
            Javg(1) = Jinput
        endif

        write(*,'(A)') 'Obtained the following interaction values:'
        write(*,'(A,F7.3,A,F7.3)') ' Uavg(0)  = ',Uavg(1), ', Javg(0)  = ',Javg(1)
        write(*,'(A,F7.3,A,F7.3)') ' Uinput    = ',Uinput,   ', Jinput    = ',Jinput
        if (useUw==1) then
            write(*,'(A)') ' U(inu_n) is frequency dependent!'
        else
            write(*,'(A)') ' U is constant! No frequency dependent U calculation!'
        endif
        write(*,'(A)') ' '

   end subroutine obtain_Uavg_Javg

!  ============================================================
!  ============================================================
    subroutine add_to_diag(v_xc, nk, value)
        use params
!        use matsfunc_ops
        complex(kr),intent(inout) :: v_xc(:,:,:,:)   ! norb,norb,nspin,nkpts
        integer,intent(in)        :: nk
        real(kr),intent(in)       :: value

        integer   :: m,s,k

        do k=1,nk
            do s=1,nspin
                do m=1,norb
                    v_xc(m,m,s,k) = v_xc(m,m,s,k) + value
                enddo
            enddo
        enddo

    end subroutine add_to_diag

!  ============================================================
!  == Write DCA kpatch-index map into file
!  ============================================================
   subroutine write_kPatchIndex(filename)
      use params
      character(len=*), intent(in) :: filename
      integer :: iounit=10, ikx,iky,i, effnkpts
      
      open(unit=iounit,file=filename,status="unknown")
      do ikx=0,nkx-1
         do iky=0,nky-1
            write(iounit,'(I4,3X,I4,3X,I4)') ikx,iky,kPatchIndex(ikx*nky+iky+1)
         enddo
         write(iounit,'(A)') ' '
      enddo

      effnkpts = 0
      write(*,'(A)') ' Number of k-Points in each k-patch:'
      do ikx=0,npx-1
         do iky=0,npy-1
            write(*,'(I4,3X,I4,3X,I4)') ikx,iky,nKpatch(ikx*npy+iky+1)
            effnkpts = effnkpts + nKpatch(ikx*npy+iky+1)
         enddo
         write(iounit,'(A)') ' '
      enddo
      write(*,'(A,I6)') ' effnkpts=',effnkpts
      
      
      close(iounit)
   end subroutine write_kPatchIndex

!  ============================================================
!  == Write the full nonlocal Selfenergy into a file that can be used as input
!  (use hartree units, etc)
!  ============================================================
   subroutine write_sgw_input(s_gw)
      use params
      complex(kr),intent(in) :: s_gw(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega
      integer                :: m1,m2,s,w,i,ik,ikx,iky,ikz,iounit=10

      write(*,'(A)') ' Write the GW Selfenergy into a file gw_sigma_out.dat'
      open(unit=iounit,file="gw_sigma_out.dat",status="unknown")
      write(iounit,'(A)') "# Sxc_wann interpolated to Matsubara freq."
      ik = 1
      do ikx=0,nkx-1
         do iky=0,nky-1
            do ikz=0,nkz-1
               write(iounit,'(A,I6,F11.7,2X,F11.7,2X,F11.7,A)') "# ",ik,ikx*1.0/nkx,iky*1.0/nky,ikz*1.0/nkz," 0.0   0.0 "

                do w=1,size(s_gw(1,1,1,1,:))
                   do s=1,gw_nspin
                      do m2=1,size(s_gw(:,1,1,1,1))
                         do m1=1,size(s_gw(1,:,1,1,1))
                            write(iounit,'(I5,3X)', advance='no') w
                            write(iounit,'(ES23.16,3X)', advance='no') real( s_gw(m1,m2,s,ik,w) )/HartreeToEV
                            write(iounit,'(ES23.16,3X)', advance='no') aimag(s_gw(m1,m2,s,ik,w) )/HartreeToEV
                         enddo !m2 loop
                      enddo !m1 loop
                   enddo !s loop
                   write(iounit,'(A)') ' '
                enddo ! w loop
                write(iounit,'(A)') ' '
   
               ik = ik+1
            enddo !ikz
         enddo !iky
      enddo !ikx
      close(iounit)

   end subroutine write_sgw_input


end module io_ops
!  ============================================================
!  End of module
!  ============================================================

