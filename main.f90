!  ============================================================
!  == Main program: DMFT selfconsistency + GW 
!  ============================================================

program main
   use constants
   use params
   use io_ops
   use matsfunc_ops
   use gw
   use observables
   use solver
   implicit none
   !include 'mpif.h'
   complex(kr),allocatable :: h_dft(:,:,:,:)        ! norb,norb,nspin,nkpts
   complex(kr),allocatable :: h_dft_sig0(:,:,:,:)   ! norb,norb,nspin,nkpts
   complex(kr),allocatable :: v_xc(:,:,:,:)         ! norb,norb,nspin,nkpts
   complex(kr),allocatable :: gloc(:,:,:,:),      & ! norb,norb,nspin,nomega
                            & gimp(:,:,:,:),      & ! norb,norb,nspin,nomega
                            & simp(:,:,:,:),      & ! norb,norb,nspin,nomega
                            & gimp_new(:,:,:,:),  & ! norb,norb,nspin,nomega
                            & simp_new(:,:,:,:),  & ! norb,norb,nspin,nomega
                            & simp_old(:,:,:,:),  & ! norb,norb,nspin,nomega
                            & gbath(:,:,:,:),     & ! norb,norb,nspin,nomega
                            & hybrid(:,:,:,:)       ! norb,norb,nspin,nomega
   complex(kr),allocatable :: S_gw(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega
   complex(kr),allocatable :: s_gw_dc(:,:,:,:),   & ! norb,norb,nspin,nomega
                            & s_gw_loc(:,:,:,:)     ! norb,norb,nspin,nomega

   real(kr),allocatable    :: dft_hartree(:,:,:),    & ! norb,norb,nspin
                            & dft_hartree_dc(:,:,:), & ! norb,norb,nspin
                            & dft_exchange(:,:,:),   & ! norb,norb,nspin
                            & dft_exchange_dc(:,:,:)   ! norb,norb,nspin

   !!! Stuff needed for GW 
   complex(kr),allocatable :: gf_tau(:,:,:,:,:)   ! norb,norb,nspin,nkpts,ntau

   complex(kr),allocatable :: P_gw(:,:,:,:,:),  & ! norb**2,norb**2,nspin,nkpts,nnu
                            & p_gw_dc(:,:,:,:)    ! norb**2,norb**2,nspin,nnu
   complex(kr),allocatable :: W_gw(:,:,:,:,:)     ! norb**2,norb**2,nspin,nkpts,nnu
   complex(kr),allocatable :: V_gw(:,:,:,:,:)     ! norb**2,norb**2,nspin,nkpts,nnu
   complex(kr),allocatable :: wloc(:,:,:,:),    & ! norb**2,norb**2,nspin,nnu
                            & wimp(:,:,:,:),    & ! norb**2,norb**2,nspin,nnu
                            & pimp(:,:,:,:),    & ! norb**2,norb**2,nspin,nnu
                            & uloc(:,:,:,:)       ! norb**2,norb**2,nspin,nnu

   real(kr),allocatable   :: mu_loc(:,:,:)         ! norb,norb,nspin

   integer(ki)                       :: iiter,m1,m2,s,n
!   integer(kind = 4)                 :: ntasks_mpi,rank_mpi,ierr_mpi,len_mpi
!   character(MPI_MAX_PROCESSOR_NAME) :: hostname

!  ============================================================
!  == Initialize MPI
!  ============================================================
   !call MPI_INIT(ierr_mpi)
   !call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks_mpi, ierr_mpi)
   !call MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr_mpi)
   !call MPI_GET_PROCESSOR_NAME(hostname, len_mpi, ierr_mpi)
   !write(*,'(A,I,I,A)') 'Starting MPI task ',rank_mpi,' out of ',ntasks_mpi,', running on ',hostname

!  ============================================================
!  == Read basic parameters
!  ============================================================

   call read_params()

!  ============================================================
!  == Allocation and read-in of DFT and GW results
!  ============================================================
   allocate( h_dft     (norb,norb,nspin,nkpts) );        h_dft        = (0.0_kr,0.0_kr)
   allocate( h_dft_sig0(norb,norb,nspin,nkpts) );        h_dft_sig0   = (0.0_kr,0.0_kr)
   allocate( gloc      (norb,norb,nspin,nomega) );       gloc         = (0.0_kr,0.0_kr)
   allocate( gimp      (norb,norb,nspin,nomega) );       gimp         = (0.0_kr,0.0_kr)
   allocate( simp      (norb,norb,nspin,nomega) );       simp         = (0.0_kr,0.0_kr)
   allocate( gimp_new  (norb,norb,nspin,nomega) );       gimp_new     = (0.0_kr,0.0_kr)
   allocate( simp_new  (norb,norb,nspin,nomega) );       simp_new     = (0.0_kr,0.0_kr)
   allocate( simp_old  (norb,norb,nspin,nomega) );       simp_old     = (0.0_kr,0.0_kr)
   allocate( gbath     (norb,norb,nspin,nomega) );       gbath        = (0.0_kr,0.0_kr)
   allocate( hybrid    (norb,norb,nspin,nomega) );       hybrid       = (0.0_kr,0.0_kr)
   allocate( s_gw_dc   (norb,norb,nspin,nomega) );       s_gw_dc      = (0.0_kr,0.0_kr)
   allocate( s_gw_loc  (norb,norb,nspin,nomega) );       s_gw_loc     = (0.0_kr,0.0_kr)
   allocate( dft_hartree(norb,norb,nspin) );             dft_hartree  = (0.0_kr)
   allocate( dft_exchange(norb,norb,nspin) );            dft_exchange = (0.0_kr)
   allocate( dft_hartree_dc(norb,norb,nspin) );          dft_hartree_dc  = (0.0_kr)
   allocate( dft_exchange_dc(norb,norb,nspin) );         dft_exchange_dc = (0.0_kr)

   allocate( wloc   (norb**2,norb**2,nspin,nnu) );               wloc   = (0.0_kr)
   allocate( wimp   (norb**2,norb**2,nspin,nnu) );               wimp   = (0.0_kr)
   allocate( pimp   (norb**2,norb**2,nspin,nnu) );               pimp   = (0.0_kr)
   allocate( p_gw_dc(norb**2,norb**2,nspin,nnu) );               p_gw_dc   = (0.0_kr)

   if ( useSigK/=0 ) then
      allocate( v_xc(norb,norb,nspin,nkpts) );        v_xc         = (0.0_kr,0.0_kr)
      allocate( S_gw(norb,norb,nspin,nkpts,nomega) ); S_gw         = (0.0_kr,0.0_kr)
   else
      allocate( v_xc(norb,norb,nspin,1) );                    v_xc         = (0.0_kr,0.0_kr)
      allocate( S_gw(norb,norb,nspin,1,1) );                  S_gw         = (0.0_kr,0.0_kr)
   endif
   if ( usePolK/=0 .or. useSigK==2) then
      allocate( P_gw(norb**2,norb**2,nspin,nkpts,nnu) );        P_gw         = (0.0_kr,0.0_kr)
      allocate( W_gw(norb**2,norb**2,nspin,nkpts,nnu) );        W_gw         = (0.0_kr,0.0_kr)
      allocate( V_gw(norb**2,norb**2,nspin,nkpts,nnu) );        V_gw         = (0.0_kr,0.0_kr)
      allocate( gf_tau(norb,norb,nspin,nkpts,ntau+1)) ;         gf_tau   = (0.0_kr)
   else
      allocate( V_gw(norb**2,norb**2,nspin,1,1) );        V_gw         = (0.0_kr,0.0_kr)
      allocate( P_gw(norb**2,norb**2,nspin,1,1) );        P_gw         = (0.0_kr,0.0_kr)
      allocate( W_gw(norb**2,norb**2,nspin,1,1) );        W_gw         = (0.0_kr,0.0_kr)
      allocate( gf_tau(norb,norb,nspin,1,1)) ;         gf_tau   = (0.0_kr)
   endif
   if ( useUw/=0 .or. usePolK/=0 ) then
      allocate( uloc(norb**2,norb**2,nspin,nnu) );        uloc         = (0.0_kr,0.0_kr)
   else
      allocate( uloc(norb**2,norb**2,nspin,1) );        uloc         = (0.0_kr,0.0_kr)
      ! we only need to store the value at w=0
   endif

   allocate( mu_loc(norb,norb,nspin) );                           mu_loc     = (0.0_kr)
   allocate( gw_kvecs(gw_nkpts,3) );                         gw_kvecs   = (0.0_kr)
   allocate( Utrafo(norbPerAtom,norbPerAtom,nspin) )             ; Utrafo=(0.0_kr,0.0_kr)
   allocate( UtrafoT(norbPerAtom,norbPerAtom,nspin) )            ; UtrafoT=(0.0_kr,0.0_kr)
   
   write(*,'(A,F8.3,A)') ' Allocated approx. ~',                     &
      &  ( 2*sizeof(h_dft) + sizeof(v_xc) + 10*sizeof(gloc) + sizeof(S_gw)         &
      &   +3*sizeof(P_gw) ) &
      &    /(1024.0**3),' Gbyte of memory...'

   ! Read the DFT and GW results
   call read_gw_input(h_dft,v_xc,S_gw,V_gw,P_gw)
   call print_vxc_local(v_xc)
   call get_locpart_mats(s_gw_loc,S_gw)
   call write_loc_matsfunc_matrix("s_gw_loc.dat",s_gw_loc)
!   call make_nonloc_qp(S_gw)

   call write_kdep_static_f("h_interpolated.dat",h_dft,1)
   if ( useSigK/=0 )then
      h_dft_sig0 = h_dft - v_xc + REAL(S_gw(:,:,:,:,1))
      call write_kdep_static_f("h_sig0_interpolated.dat",h_dft_sig0,1)
      h_dft_sig0 = h_dft - v_xc
      call write_kdep_static_f("h_-vxc_interpolated.dat",h_dft_sig0,1)

      call write_kdep_static_f("vxc_interpolated.dat",v_xc,0)
      !call write_kdep_static_f("sgw_interpolated.dat",S_gw(:,:,:,:,1),0)
      call write_kdep_w0_f("sgw_interpolated.dat",S_gw)
   endif

!  ============================================================
!  == Prepare everything for SC cycle =========================
!  ============================================================
   write(*,'(A)') 'WARNING: MU_LOC IS ASSUMED TO BE REAL!!!! CHANGE THIS !!!!!!'

   ! First calculate the DFT Green's function since we need it for
   ! the DFT Hartree and Fock terms and the Doublecounting
   ! s_gw_dc, dft_hartree_dc, dft_exchange_dc are still zero here
   write(*,'(A)') 'Calculate initial DFT local Greens function...'
   call get_gloc(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.true.)
   write(*,'(A)') 'Check chemical potential for correct filling...'
   call adjust_chem_potential(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,   &
                         &    dft_hartree_dc,dft_exchange_dc,.true.)
   call print_local_filling("G_loc DFT",gloc)
   call check_offdiagonal("G_loc DFT",gloc)
   call print_Etot(h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.true.)
   call write_loc_matsfunc_matrix("g_loc_DFT.dat",gloc)
   call write_pade_loc(gloc,"g_loc_DFT_pade.dat",.true.)
   call write_pade_bandstructure("bandstruct_DFT.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.true.)
   call write_fermisurface("fermisurface_DFT.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.true.)
   call write_fermisurface_sinterp("fermisurface_DFT_sinterp.dat",h_dft,v_xc,S_gw,s_gw_dc,simp, &
                                                        & dft_hartree_dc,dft_exchange_dc,.true.)
   call check_sym("G_loc_DFT init",gloc)

   ! Setup the bare interaction
   call setup_bare_interaction(V_gw)
   call write_kdep_fs("V_bare_fs.dat",V_gw,0)
   call get_uloc(uloc,P_gw,pimp,p_gw_dc,W_gw,V_gw,.true.)                      ! just get the bare interactions
   call print_umatrix(uloc)
   call write_loc_matsfunc_bos_matrix("u_loc.dat",uloc)

   ! Calculate the Hartree and Fock exchange terms from DFT
   call get_hartree_fock(dft_hartree,dft_exchange,gloc,uloc)
   call print_dft_hartree_exchange(dft_hartree,dft_exchange)
   call get_hartree_fock_dc(dft_hartree_dc,dft_exchange_dc,dft_hartree,dft_exchange)

   ! if we need it, get the polarization from the noninteracting Green's function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call polarization_gg(P_gw,gf_tau,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree,dft_exchange,.true.)
   call write_kdep_fs("pol_G0G0_fs.dat",P_gw,0)
   call write_locpart_matsfunc_bos("p_G0G0_loc.dat",P_gw)
   call write_loc_matsfunc_bos_matrix("p_G0G0_Gamma.dat",P_gw(:,:,:,1,:))

   call get_W(W_gw,P_gw,pimp,p_gw_dc,V_gw)
   call write_kdep_fs("W_G0G0_fs.dat",W_gw,0)
   call write_locpart_matsfunc_bos("W_G0G0_loc.dat",W_gw)
   
   call sigma_gw(S_gw,W_gw,gf_tau)
   call write_kdep_fs("sigma_G0W0_fs.dat",S_gw,1)
   call write_locpart_matsfunc("S_G0W0_loc.dat",S_gw)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  ============================================================
!  Now we have set up all raw objects where we need the Gloc from DFT
!  Start filling simp and pimp with meaningful values and calculate Doublecounting
!  ============================================================

   if (contcalc==1) then
      call read_old_data(simp,pimp,s_gw_dc,p_gw_dc,dft_hartree,dft_exchange)
      call read_old_dca(S_gw,s_gw_dc)
      call print_dft_hartree_exchange(dft_hartree,dft_exchange)
      call get_hartree_fock_dc(dft_hartree_dc,dft_exchange_dc,dft_hartree,dft_exchange)
      call symmetrize_matsubara("s_imp_old",simp)
   else
      call get_sigma_gw_doublecounting(s_gw_dc,S_gw,gloc,wloc)
      call get_pol_doublecounting(p_gw_dc,P_gw,gloc)
      simp = s_gw_dc
      pimp = p_gw_dc
      call adjust_to_hartree_fock(simp,dft_hartree_dc,dft_exchange,dft_exchange_dc)
      if (contcalc==2) then
         call init_insulating(simp)
      endif

      call write_loc_matsfunc_matrix("s_gw_dc.dat",s_gw_dc)
      call write_loc_matsfunc_bos_matrix("p_gw_dc.dat",p_gw_dc)
      call write_hartree_exchange("dft_hartree_exchange.dat",dft_hartree,dft_exchange)
   endif
!!!!!!!!!!!!!!
!call reperiod_sigma(s_gw,simp)
!!!!!!!!!!!!!!!!!!!
   call check_sym("s_imp_init.dat",simp)

   ! Whatever we read in the input, make sure the DC is correct because SigK might be new
   call get_sigma_gw_doublecounting(s_gw_dc,S_gw,gloc,wloc)
   call get_pol_doublecounting(p_gw_dc,P_gw,gloc)

   call write_loc_matsfunc_matrix("s_imp_init.dat",simp)
   call write_loc_sigmas(s_gw_loc,s_gw_dc,simp)
   call write_loc_matsfunc_bos_matrix("p_imp_init.dat",pimp)

!do iiter=1,niter
   write(*,'(A)') 'Calculate initial interacting local Greens function...'
   call get_gloc(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   write(*,'(A)') 'Check chemical potential for correct filling...'
   call adjust_chem_potential(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,   &
                         &    dft_hartree_dc,dft_exchange_dc,.false.)
   call print_local_filling("G_loc init",gloc)
   call check_offdiagonal("G_loc init",gloc)
   call print_Etot(h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_loc_matsfunc_matrix("g_loc_init.dat",gloc)
   call write_gfsig_nonloc(1,0,0,"x",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_gfsig_nonloc(1,1,0,"xy",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_gfsig_nonloc(2,0,0,"2x",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_gfsig_nonloc(2,2,0,"2x2y",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_pade_loc(gloc,"g_loc_init_pade.dat",.true.)
   call write_pade_bandstructure("bandstruct_init.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
!   call write_poly_bandstructure("bandstruct_init_poly.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_fermisurface("fermisurface_init.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
   call write_fermisurface_sinterp("fermisurface_init_sinterp.dat",h_dft,v_xc,S_gw,s_gw_dc,simp, & 
                                                           & dft_hartree_dc,dft_exchange_dc,.false.)
   call check_sym("G_loc init",gloc)

   write(*,'(A)') 'Calculate initial local interaction W(inu) from G_init...'
   call polarization_gg(P_gw,gf_tau,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree,dft_exchange,.false.)
   call get_pol_doublecounting(p_gw_dc,P_gw,gloc)
   call write_loc_matsfunc_bos_matrix("p_gw_dc.dat",p_gw_dc)
   call write_kdep_fs("pol_GG_fs.dat",P_gw,0)
   call write_locpart_matsfunc_bos("p_GG_loc.dat",P_gw)
   call write_loc_matsfunc_bos_matrix("p_GG_Gamma.dat",P_gw(:,:,:,1,:))

   call get_W(W_gw,P_gw,pimp,p_gw_dc,V_gw)
   call write_kdep_fs("W_fs.dat",W_gw,0)
   call write_locpart_matsfunc_bos("W_loc.dat",W_gw)
   call get_locpart_mats(wloc,W_gw)
   
   call sigma_gw(S_gw,W_gw,gf_tau)
   call get_sigma_gw_doublecounting(s_gw_dc,S_gw,gloc,wloc)
   call write_loc_matsfunc_matrix("s_gw_dc.dat",s_gw_dc)
   call write_kdep_fs("sigma_GW_fs.dat",S_gw,1)
   call write_locpart_matsfunc("S_GW_loc.dat",S_gw)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!simp = s_gw_dc
!enddo

   ! Initialize G_imp with G_loc since G_imp will be mixed
   do iiter=1,norb
      gimp(iiter,iiter,:,:) = gloc(iiter,iiter,:,:)
      simp_old(iiter,iiter,:,:) = simp(iiter,iiter,:,:)
   enddo
   simp_new = simp
   simp_old = simp

!  ======================================================================
!  ======================================================================
!  == Start the SC cycle ================================================
!  ======================================================================
!  ======================================================================

   write(*,'(A,/)') 'Starting the selfconsistency cycle...'

   do iiter=1,niter
      write(*,'(A11,I4)') 'Iteration: ',iiter
   
     !  ============================================================
     !  == Calculate the effective bath ============================
     !  ============================================================
      call neglect_offdiag(gloc)
      write(*,'(A)') ' Calculate Weiss field and local orbital levels...'
      if (bathMethod==0) then
         call get_weiss_field_loc_levels(gbath,mu_loc,gloc,h_dft,v_xc,S_gw,s_gw_dc,simp, &
                                         & dft_hartree_dc,dft_exchange_dc,.false.)
      else
         call get_weiss_field_loc_levels_dyson(gbath,mu_loc,gloc,simp)
      endif
      call check_sym("G_bath",gbath)
      call check_offdiagonal("G_bath",gbath)
      call write_loc_matsfunc_matrix("g_bath.dat",gbath)
      call write_pade_loc(gbath,"g_bath_pade.dat",.true.)
      call print_local_orbital_levels(mu_loc)
      write(*,'(A)') ' Calculate hybridization function...'
      call get_hybrid(hybrid,gbath,mu_loc)
      call check_sym("hybrid",hybrid)

      call check_offdiagonal("Hybridization",hybrid)
      call print_highfreq_coeffs("Hybridization",hybrid)
      call write_loc_matsfunc_matrix("hybrid.dat",hybrid)
      call write_pade_loc(hybrid,"hybrid_pade.dat",.true.)


      !  ============================================================
      !  == Calculate the effective interaction =====================
      !  ============================================================
      write(*,'(A)') ' Obtain the effective interaction...'
      call get_uloc(uloc,P_gw,pimp,p_gw_dc,W_gw,V_gw,.false.)
      call print_umatrix(uloc)
      call write_loc_matsfunc_bos_matrix("u_loc.dat",uloc)

     !  ============================================================
     !  == Solve the impurity model ================================
     !  ============================================================

      write(*,'(A)') ' Write input for impurity solver...'
      call write_solver_input(uloc,hybrid,gbath,mu_loc)
      write(*,'(A)') ' Starting impurity solver...'
      if (npx*npy .gt. 0) then
         call solver_DCA(gimp_new,simp_new,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree,dft_exchange,.false.)
      else
         call run_solver()
         write(*,'(A)') ' Impurity solver finished, reading output...'
         call read_loc_matsfunc_matrix_solver("Gwl.dat",gimp_new)
   
         call read_wimp_pimp("nnw.dat","nnt.dat",pimp,wimp,uloc)
         call symmetrize_matsubara("G_imp",gimp_new)
         call check_sym("g_imp after symmetrization",gimp_new)
         call symmetrize_matsubara("G_bath",gbath)
         call check_sym("g_bath after symmetrization",gbath)  ! Since we use gbath and gimp both have to be symmetrized
         call get_sigma_dyson(simp_new,gimp_new,gbath)
   !      call copy_sigma_afm(simp_new) 
!      call hf_solver(S_gw,h_dft,umatrix,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc)
      endif
      
      simp = (1-mixing)*simp + mixing*simp_new
      gimp = (1-mixing)*gimp + mixing*gimp_new

!      call check_sym("s_imp before reperiod",simp)
      ! reperiodize the Selfenergy in Cluster calculations (lots of hardcoding needed)
!      call reperiod_sigma(S_gw,simp)
!      call check_sym("s_imp after reperiod",simp)

      if (updateHFdc==1) then
         ! ReCalculate the Hartree and Fock exchange terms from DFT
         write(*,*) ' Update the Hartree-Fock Doublecounting from G_imp...'
         call get_hartree_fock(dft_hartree,dft_exchange,gimp,uloc)
         call print_dft_hartree_exchange(dft_hartree,dft_exchange)
         call get_hartree_fock_dc(dft_hartree_dc,dft_exchange_dc,dft_hartree,dft_exchange)
         call write_hartree_exchange("dft_hartree_exchange.dat",dft_hartree,dft_exchange)
      endif

     !  ============================================================
     !  == Update the Gloc and Wloc ================================
     !  ============================================================
      write(*,'(A)') ' Calculate new local Greens Gloc function and interaction Wloc...'
      call get_gloc(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      write(*,'(A)') 'WARNING: At some point we should also recalculate the doublecounting...'
            
      write(*,'(A)') ' Check chemical potential for constant filling...'
      call adjust_chem_potential(gloc,h_dft,v_xc,S_gw,s_gw_dc,simp,   &
                            &    dft_hartree_dc,dft_exchange_dc,.false.)
      call check_offdiagonal("G_loc",gloc)
      call check_sym("new g_loc",gloc)

     !  ============================================================
     !  == Update the polarization, interactons, etc ===============
     !  ============================================================
     call polarization_gg(P_gw,gf_tau,h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree,dft_exchange,.false.)
     call get_pol_doublecounting(p_gw_dc,P_gw,gloc)
     call write_loc_matsfunc_bos_matrix("p_gw_dc.dat",p_gw_dc)
     call write_kdep_fs("pol_GG_fs.dat",P_gw,0)
     call write_locpart_matsfunc_bos("p_GG_loc.dat",P_gw)
     call write_loc_matsfunc_bos_matrix("p_GG_Gamma.dat",P_gw(:,:,:,1,:))
   
     call get_W(W_gw,P_gw,pimp,p_gw_dc,V_gw)
     call write_kdep_fs("W_fs.dat",W_gw,0)
     call write_locpart_matsfunc_bos("W_loc.dat",W_gw)
     call get_locpart_mats(wloc,W_gw)

     call sigma_gw(S_gw,W_gw,gf_tau)
     call get_sigma_gw_doublecounting(s_gw_dc,S_gw,gloc,wloc)
     call write_loc_matsfunc_matrix("s_gw_dc.dat",s_gw_dc)
     call write_kdep_fs("sigma_gw_fs.dat",S_gw,1)
     call write_locpart_matsfunc("S_GW_loc.dat",S_gw)

!  === Calculate and print information about convergence =============

      call print_local_filling("G_bath",gbath)
      call print_local_filling("G_imp",gimp)
      call print_local_filling("G_loc",gloc)      

      call print_matsf_diff(gloc,gimp,'G_loc','G_imp')
      call print_matsf_diff(simp_old,simp_new,'S_imp old','S_imp new')
      call print_matsf_diff_atw(simp_old,simp_new,1,'S_imp old','S_imp new')
      call print_matsf_diff_atw(simp_old,simp_new,nomega,'S_imp old','S_imp new')
      call print_matsf_diff(wloc,wimp,'W_loc','W_imp')
      simp_old = simp

      ! Write all data into files
      call print_Etot(h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)

      call write_loc_matsfunc_matrix("s_imp.dat",simp)
      call write_loc_matsfunc_matrix("g_loc.dat",gloc)
      call write_loc_matsfunc_matrix("g_imp.dat",gimp)
      call write_loc_matsfunc_matrix("g_bath.dat",gbath)
      call write_gfsig_nonloc(1,0,0,"x",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_gfsig_nonloc(1,1,0,"xy",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_gfsig_nonloc(2,0,0,"2x",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_gfsig_nonloc(2,2,0,"2x2y",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_loc_sigmas(S_gw_loc,s_gw_dc,simp)
      call write_loc_matsfunc_matrix("s_gw_dc.dat",s_gw_dc)
!      call write_loc_matsfunc_bos_real("p_imp.dat",pimp)
!      call write_loc_matsfunc_bos_real("w_loc.dat",wloc)
!      call write_loc_matsfunc_bos_real("w_imp.dat",wimp)
!      call write_loc_matsfunc_bos_real("u_loc.dat",uloc)
      call write_loc_matsfunc_bos_matrix("p_imp.dat",pimp)
      call write_loc_matsfunc_bos_matrix("w_imp.dat",wimp)
      call write_mu("mu.dat")
      call write_pade_loc(gloc,"g_loc_pade.dat",.true.)
      call write_pade_loc(gimp,"g_imp_pade.dat",.true.)
      call write_pade_loc(simp,"s_imp_pade.dat",.true.)
      call write_pade_loc(hybrid,"hybrid_pade.dat",.true.)
!      call write_pade_loc(wloc,"w_loc_pade.dat",.false.)
!      call write_pade_loc(wimp,"w_imp_pade.dat",.false.)
!      call write_pade_loc(uloc,"u_loc_pade.dat",.false.)
!      call write_pade_loc(pimp,"p_imp_pade.dat",.false.)
      call write_pade_bandstructure("bandstruct.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
!      call write_poly_bandstructure("bandstruct_poly.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_fermisurface("fermisurface.dat",h_dft,v_xc,S_gw,s_gw_dc,simp,dft_hartree_dc,dft_exchange_dc,.false.)
      call write_fermisurface_sinterp("fermisurface_sinterp.dat",h_dft,v_xc,S_gw,s_gw_dc,simp, & 
                                                     & dft_hartree_dc,dft_exchange_dc,.false.)
      write(*,'(A)') ' '
   enddo
   write(*,'(A)') 'Maximum number of DMFT iterations reached, finish program...'

!  ============================================================
!  == Deallocate and cleanup ==================================
!  ============================================================

   deallocate( h_dft )
   deallocate( h_dft_sig0 )
   deallocate( gloc )
   deallocate( gimp )
   deallocate( simp )
   deallocate( gimp_new )
   deallocate( simp_new )
   deallocate( simp_old )
   deallocate( gbath )
   deallocate( hybrid )
   deallocate( dft_hartree )
   deallocate( dft_exchange )
   deallocate( dft_hartree_dc )
   deallocate( dft_exchange_dc )
   deallocate( s_gw_dc )
   deallocate( s_gw_loc )

   deallocate( v_xc )
   deallocate( S_gw ) 
      deallocate( gf_tau )

   deallocate( P_gw )
   deallocate( W_gw )
   deallocate( V_gw )
   deallocate( uloc )
   deallocate( wloc   )
   deallocate( wimp   )
   deallocate( pimp   )
   deallocate( p_gw_dc)
!   deallocate( Uavg )
!   deallocate( Javg )

   deallocate( mu_loc )
!   deallocate( umatrix )
   deallocate( gw_kvecs )
   deallocate( Utrafo )
   deallocate( UtrafoT )

   ! This has been allocated when reading the input in io_ops.f90
   deallocate( dmftOrbsIndex )
!   deallocate( degOrbs%noDegenOrbs )
!   deallocate( degOrbs%degenOrbsIndex )
   
end program main
