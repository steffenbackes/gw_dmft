!  ============================================================
!  == Mathematical Operations 
!  ============================================================

module mathfunc
   use constants
   use params
   implicit none

 contains
!  ============================================================
!  == Return fermionic Matsubara frequencies: start at n=0
!  ============================================================
   function wn(n)
      real(kr)               :: wn
      integer(ki),intent(in) :: n
      wn = (2*n+1)*pi/beta
   end function wn

!  ============================================================
!  == Return bosonic Matsubara frequencies: start at n=0
!  ============================================================
   function vn(n)
      real(kr)               :: vn
      integer(ki),intent(in) :: n
      vn = 2*n*pi/beta
   end function vn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! matrix diagonalization !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_eigenvectors(evals,evectors,matrix,N)
      real(kr),intent(inout) :: evals(:)
      complex(kr),intent(inout) :: evectors(:,:)
      complex(kr),intent(in)    :: matrix(:,:)
      integer(ki),intent(in)    :: N 
      !   Lapack stuff
      COMPLEX*16,allocatable       :: matrix_zheev(:,:)  ! matrix for  CZHEEV
      INTEGER                      :: N_zheev, LDA_zheev, LWORK_zheev,info_lpck
      DOUBLE PRECISION,allocatable :: eigvals(:), rwork(:)
      COMPLEX*16,allocatable       :: work(:)
      external ZHEEV

      N_zheev = N
      LDA_zheev = N
      LWORK_zheev = N*2
      allocate( matrix_zheev(LDA_zheev,LDA_zheev) )
      allocate( eigvals(LDA_zheev) )
      allocate( rwork(3*LDA_zheev-2) )
      allocate( work(LWORK_zheev) )

      ! copy matrix
      matrix_zheev = matrix

      ! diagonalize
      call ZHEEV( 'V', 'U', N_zheev, matrix_zheev, LDA_zheev, eigvals, work, &
               &  LWORK_zheev, rwork, info_lpck  )

      if ( info_lpck /= 0) stop 'ERROR: Could not diagonalize matrix !!!'

!      write(*,*) 'Matrix:'
!      write(*,*)  matrix
!      write(*,*) 'Eigenvalues: '
!      write(*,*) eigvals
!      write(*,*) 'Eigenvectors'
!      write(*,*) matrix_zheev

      ! copy to result
      evals = eigvals
      evectors = matrix_zheev

      deallocate( matrix_zheev )
      deallocate( eigvals )
      deallocate( rwork )
      deallocate( work )
   end subroutine get_eigenvectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! matrix diagonalization single prec !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_eigenvectors_sp(evals,evectors,matrix,N)
      real(kr),intent(inout) :: evals(:)
      complex(kr),intent(inout) :: evectors(:,:)
      complex(kr),intent(in)    :: matrix(:,:)
      integer(ki),intent(in)    :: N 
      !   Lapack stuff
      COMPLEX,allocatable       :: matrix_zheev(:,:)  ! matrix for  CZHEEV
      INTEGER                   :: N_zheev, LDA_zheev, LWORK_zheev,info_lpck
      real,allocatable          :: eigvals(:), rwork(:)
      COMPLEX,allocatable       :: work(:)
      external CHEEV

      N_zheev = N
      LDA_zheev = N
      LWORK_zheev = N*2
      allocate( matrix_zheev(LDA_zheev,LDA_zheev) )
      allocate( eigvals(LDA_zheev) )
      allocate( rwork(3*LDA_zheev-2) )
      allocate( work(LWORK_zheev) )

      ! copy matrix
      matrix_zheev = matrix

      ! diagonalize
      call CHEEV( 'V', 'U', N_zheev, matrix_zheev, LDA_zheev, eigvals, work, &
               &  LWORK_zheev, rwork, info_lpck  )

      if ( info_lpck /= 0) stop 'ERROR: Could not diagonalize matrix !!!'

!      write(*,*) 'Matrix:'
!      write(*,*)  matrix
!      write(*,*) 'Eigenvalues: '
!      write(*,*) eigvals
!      write(*,*) 'Eigenvectors'
!      write(*,*) matrix_zheev

      ! copy to result
      evals = eigvals
      evectors = matrix_zheev

      deallocate( matrix_zheev )
      deallocate( eigvals )
      deallocate( rwork )
      deallocate( work )
   end subroutine get_eigenvectors_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! complex exact matrix inversion !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_inverse(invMat,matrix,N)
      complex(kr),intent(inout) :: invMat(:,:)  ! N,N
      complex(kr),intent(in)    :: matrix(:,:)  ! N,N
      integer(ki),intent(in)    :: N 

  !   Lapack stuff
      complex(8),allocatable :: mat_tmp(:,:)  ! N,N

      complex(8), dimension(N) :: work  ! work array for LAPACK
      integer, dimension(N)    :: ipiv  ! pivot indices
      integer :: info_lpck
      external ZGETRF
      external ZGETRI

      allocate( mat_tmp(N,N) )

      ! copy matrix
      mat_tmp = matrix

      ! Then invert with LAPACK
      call ZGETRF(N, N, mat_tmp, N, ipiv, info_lpck)
      if (info_lpck /= 0) then
         write(*,'(A,I3)') 'ERROR: Matrix for inversion is numerically singular! Return value', &
                                       & info_lpck
         stop 
      end if

      call ZGETRI(N, mat_tmp, N, ipiv, work, N, info_lpck)
      if (info_lpck /= 0) then
         stop 'Matrix inversion failed!'
      end if

      ! copy to output
      invMat = mat_tmp

      deallocate( mat_tmp )

   end subroutine get_inverse

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function isnan_c(a) ! complex version
      logical                :: isnan_c
      complex(kr),intent(in) :: a

      if (a.ne.a) then
         isnan_c = .true.
      else
         isnan_c = .false.
      end if
      return
   end function isnan_c

   function isnan_r(a) ! real version
      logical                :: isnan_r
      real(kr),intent(in) :: a

      if (a.ne.a) then
         isnan_r = .true.
      else
         isnan_r = .false.
      end if
      return
   end function isnan_r

   function isnan_i(a) ! integer version
      logical                :: isnan_i
      integer(ki),intent(in) :: a

      if (a.ne.a) then
         isnan_i = .true.
      else
         isnan_i = .false.
      end if
      return
   end function isnan_i

   function isinf_r(a) ! real version
      logical             :: isinf_r
      real(kr),intent(in) :: a

      if ((a*0.0_kr) .ne. 0.0_kr) then
         isinf_r = .true.
      else
         isinf_r = .false.
      end if
      return
   end function isinf_r

!  ============================================================
!  == Calculate the matrix maximum norm
!  ============================================================
   function get_matrixnorm_inf(mat,n)
      use params
      real(kr)                    :: get_matrixnorm_inf
      complex(kr),intent(in)      :: mat(:,:)  ! n,n
      integer(ki),intent(in)      :: n                ! dimension
      integer                     :: m1,m2
      real(kr)  :: diff
      
      diff = 0.0_kr
      do m1=1,n
         do m2=1,n
            if ( abs( mat(m1,m2) ) > diff ) then
               diff = abs( mat(m1,m2) )
            endif
         enddo !m2
      enddo !m1

      get_matrixnorm_inf = diff

      end function get_matrixnorm_inf

!  ============================================================
!  == Check orbital symmetry of GF
!  ============================================================
   subroutine check_sym(filename,gf)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       ::gf(:,:,:,:) ! norb,norb,nspin,nomega
      integer                     :: m1,m2,n,s
      
      outer: do n=1,nomega
         do s=1,2
            do m1=1,norb
               do m2=m1+1,norb
                  if ( abs( gf(m1,m2,s,n) - gf(m2,m1,s,n)  ) > 1e-5 ) then
                     write(*,'(A,A,A,I2,A,I2,A,I2,A,I2,A)') "WARNING: ",filename," not symmetric at: n=",  &
                                   &  n,", s=",s,", (m1,m2)=",m1,",",m2," !!" 
                     write(*,*) 'gf(m1,m2,s,n)=',gf(m1,m2,s,n)
                     write(*,*) 'gf(m2,m1,s,n)=',gf(m2,m1,s,n)
                     exit outer
                  endif
               enddo
            enddo !m2
         enddo !m1
      enddo outer

      end subroutine check_sym

!  ============================================================
!  == Check orbital symmetry of matrix
!  ============================================================
   subroutine check_sym_orb(filename,gf)
      use params
      character(len=*), intent(in) :: filename
      complex(kr),intent(in)       ::gf(:,:) ! norb,norb
      integer                     :: m1,m2,n,s
      
      outer: do m1=1,norb
         do m2=m1+1,norb
            if ( abs( gf(m1,m2) - gf(m2,m1)  ) > 1e-5 ) then
               write(*,'(A,A,A,A,I2,A,I2,A)') "WARNING: ",filename," not symmetric at:",  &
                             & " (m1,m2)=",m1,",",m2," !!" 
               write(*,*) 'gf(m1,m2)',gf(m1,m2)
               write(*,*) 'gf(m2,m1)',gf(m2,m1)
               exit outer
            endif
         enddo
      enddo outer 

      end subroutine check_sym_orb

!  ============================================================
!  == Return index of patch for given ikx,iky (indexing from 0)
!  ============================================================
!   function patchindex(ikx,iky)
!      use params
!      integer(ki),intent(in) :: ikx,iky
!      integer                :: ix,iy,xshift,yshift,patchindex
!
!      xshift = floor(nkpx/2.0)
!      yshift = floor(nkpy/2.0)
!!      npx = Int(nkx/nkpx)     ! #patches along x and y
!!      npy = Int(nky/nkpy) 
!
!      ix = mod( (ikx+xshift)/nkpx ,npx)
!      iy = mod( (iky+yshift)/nkpy ,npy)
!
!      patchindex = ix*npy+iy +1
!    end function patchindex

!  ============================================================
!  == Fill the index storing which kpoint belongs to which k-patch
!  ============================================================
   subroutine fillKpatchIndex(kPIndex,nKp)
      use params
      integer(ki),intent(inout) :: kPIndex(nkpts)
      integer(ki),intent(inout) :: nKp(npx*npy)

      integer(ki) :: ikpx,ikpy
      integer     :: ikx,iky
      real(kr)    :: kx,ky,kpx,kpy

      kPIndex = (0) 
      do ikpx=0,npx-1
         do ikpy=0,npy-1
            ! if quadratic grid
            kpx = ikpx*2*pi/npx 
            kpy = ikpy*2*pi/npy

            ! angle for 8 of 32 !
            if (npx .ne. npy) then
               kpx = ikpx*2*pi/npx + ikpy*2*pi/npy
               kpy = ikpy*2*pi/npy

               if (kpx+0.000001>=2*pi) then
                  kpx = kpx - 2*pi
               endif
            endif

            do ikx=0,nkx-1
               do iky=0,nky-1
                  kx = ikx*2*pi/nkx
                  ky = iky*2*pi/nky

                  ! Rotated grid uses abs norm
                  if (npx .ne. npy) then
                     if ( abs(kx-kpx)+abs(ky-kpy)<=2*pi/npy .or. abs(kx-kpx-2*pi)+abs(ky-kpy)<=2*pi/npy  .or.  &
                   &      abs(kx-kpx)+abs(ky-kpy-2*pi)<=2*pi/npy .or. abs(kx-kpx-2*pi)+abs(ky-kpy-2*pi)<=2*pi/npy ) then
   
                        kPIndex(ikx*nky+iky+1) = ikpx*npy+ikpy+1
                     endif
                  else
                     if ( max(abs(kx-kpx     ),abs(ky-kpy)     )<=pi/npy   .or.    &
                   &      max(abs(kx-kpx-2*pi),abs(ky-kpy)     )<=pi/npy   .or.    &
                   &      max(abs(kx-kpx     ),abs(ky-kpy-2*pi))<=pi/npy   .or.    &
                   &      max(abs(kx-kpx-2*pi),abs(ky-kpy-2*pi))<=pi/npy ) then
   
                        kPIndex(ikx*nky+iky+1) = ikpx*npy+ikpy+1
                     endif
                  endif

               enddo
            enddo

         enddo
      enddo

      nKp = (0) 
      do ikx=0,nkx-1
         do iky=0,nky-1
            nKp( kPIndex(ikx*nky+iky+1) ) = nKpatch( kPIndex(ikx*nky+iky+1) ) +1
         enddo
      enddo

    end subroutine fillKpatchIndex

!!!!!!!!!!!!!!!!!!!!!
! Print progress
!!!!!!!!!!!!!!!!!!!!
   subroutine progress(i,N)
      integer(ki),intent(in) :: i
      integer(ki),intent(in) :: N

      if (i .lt. N) then
          if ( ( modulo( i, max(N/20,1))==0 ) ) then        
             write(*,'(A,I3,A)',advance='no') " ",int(i*100.0_kr/N),"%"
         endif

         ! This creates to much garbage in the job.out
!         write(*,'(A,I3,A)',advance='no') " ",int(i*100.0_kr/N)," %"
!         write(*,'(6(A))',advance='no') char(8),char(8),char(8), char(8),char(8),char(8)
      else
         write(*,'(A)') " 100%, done."
      endif
      call flush(6)
   end subroutine progress

end module mathfunc
