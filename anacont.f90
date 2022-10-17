!  ============================================================
!  == Analytic continuation procedures, Pade. polynomial fit etc.
!  ============================================================

module anacont
   use constants
   implicit none
   complex(kr),allocatable :: pade_array(:,:)       ! nfreq,nfreq
   logical                 :: fermionic
   integer                 :: nfreq,maxn

   private :: pade_array,fermionic,nfreq,maxn

 contains
 !  ============================================================
 !  == Initialize the pade array with a matsubara function
 !  ============================================================
   subroutine init_pade(f,nfreq_in,fermionic_in)
      use params
      use matsfunc_ops
      complex(kr),intent(in) :: f(:)          ! nfreq
      integer,intent(in)     :: nfreq_in
      logical,intent(in)     :: fermionic_in

      integer                :: i,j
      real(kr)               :: wj,wi

      fermionic = fermionic_in
      nfreq = nfreq_in
      ! Check if we have used this already, if yes start with a fresh one
      if (allocated(pade_array)) deallocate( pade_array )
      allocate( pade_array(nfreq,nfreq) )
      pade_array = (0.0_kr,0.0_kr)

      ! Set the upper row of the pade matrix
      pade_array(1,:) = f(1:nfreq)

      ! Set remaining part of upper triangular matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      maxn = 1
      do i=1,nfreq-1
         do j=i,nfreq-1
            if (fermionic) then
               wj = wn(j)
               wi = wn(i-1)
            else
               wj = vn(j)
               wi = vn(i-1)
            endif
 
             if ( abs(pade_array(i,j+1)) < 1.0E-10 ) exit 
            
             pade_array(i+1,j+1) = ( pade_array(i,i) - pade_array(i,j+1) )    &
                       &             /( pade_array(i,j+1)*ci*( wj-wi ) )
         enddo
         if ( abs(pade_array(i+1,i+1)) > 1.0E-10)  maxn = i+1
      enddo
      ! after maxn all entries are zero, therefore we do not need to sum over
      ! them later in the get_pade() function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do j=1,nfreq-1
!         do i=1,j
!            if (fermionic) then
!               wj = wn(j)
!               wi = wn(i-1)
!            else
!               wj = vn(j)
!               wi = vn(i-1)
!            endif
!            
!            if ( abs(pade_array(i,i)-pade_array(i,j+1))< 1E-10 ) then
!               pade_array(i+1,j+1) = (0.0_kr,0.0_kr)
!            else
!               pade_array(i+1,j+1) = ( pade_array(i,i) - pade_array(i,j+1) )    &
!                       &             /( pade_array(i,j+1)*ci*( wj-wi ) )
!            endif
!         enddo
!      enddo
   end subroutine init_pade

!  ============================================================
!  == Obtain the Pade result after init_pade has been run
!  ============================================================
  function get_pade(w)
     use params
     use matsfunc_ops
     complex(kr)         :: get_pade
     real(kr),intent(in) :: w

     real(kr)            :: wi
     integer             :: n

     if (fermionic) then
        wi = wn(maxn-2)
     else
       wi = vn(maxn-2)
     endif
     get_pade = 1.0_kr + pade_array(maxn,maxn)*( w+0.001*ci - ci*wi )

     do n=maxn-2,1,-1
        if (fermionic) then
           wi = wn(n-1)
        else
           wi = vn(n-1)
        endif
        get_pade = 1.0_kr + pade_array(n+1,n+1)*( w+0.001*ci - ci*wi )/get_pade
     enddo
     get_pade = pade_array(1,1)/get_pade

  end function get_pade

end module anacont
