!  ============================================================
!  == Basic parameters accessible to anyone
!  ============================================================

module params
   use constants
   implicit none
   ! this is used to store information about the DMFT degenerate orbitals
!   type degenarateOrbitals
!      integer             :: noGroups
!      integer,allocatable :: noDegenOrbs(:)         ! noGroups
!      integer,allocatable :: degenOrbsIndex(:,:)    ! noGroups, noDegenOrbs(group)
!   end type degenarateOrbitals
   integer(ki)          :: orbSym   ! 0: no sym, 1:determine yourself, 2:  provide orbsym.dat

   ! these are variables given by the GW/DFT input from gw_input.dat
   integer(ki)          :: nomega, nnu, norb, nspin, gw_nspin
   integer(ki)          :: gw_nkpts,nqx,nqy,nqz ! nqpts, nprodstates ! product basis stuff
   real(kr),allocatable :: gw_kvecs(:,:)              ! nk from GW, (kx,ky,kz)

   ! these are variables given by the DMFT input from dmft_input.dat
   integer(ki)              :: useSigK, usePolK, useUw
   integer(ki)              :: nkx,nky,nkz, nkpts, ntau, niter, nfit=15
   integer(ki)              :: DCtype, spinorder,contcalc,readUmatrix,avgHartreeFock, neglectOffdiagBath, diagBath
   integer(ki)              :: updateHFdc
   integer(ki)              :: norb_dmft, noEquivAtoms, norbPerAtom, bathMethod
   integer(ki),allocatable  :: dmftOrbsIndex(:)
!   type(degenarateOrbitals) :: degOrbs
   real(kr)                 :: beta, nel, mu, charge_error, mixing
   real(kr)                 :: Uinput, Jinput, Unninput, U0scaleFac
   real(kr)                 :: hfield(2)        
   complex(kr),allocatable  :: Utrafo(:,:,:),UtrafoT(:,:,:)         ! Basis transformation (norbPerAtom,norbPerAtom.nspin)
   !!
   ! DCA !
   integer(ki)              :: interpDCA, npx, npy  !,nkpx,nkpy
   integer(ki),allocatable  :: kPatchIndex(:), nKpatch(:)
end module params
