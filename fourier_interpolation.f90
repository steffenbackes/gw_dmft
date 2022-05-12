!  ============================================================
!  == This module interpolates some data from a coarse k-grid to an arbitrary grid
!  ============================================================

module fourier_interpolation
   use constants
   implicit none

 contains
!  ============================================================
!  == Interpolation in k-space / Wannier-interpolation
!  ============================================================
   subroutine interpolate_kmatrix(kgrid,qgrid)
      use params
      complex(kr),intent(inout) :: kgrid(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: qgrid(:,:,:,:)       ! norb,norb,nspin,gw_nkpts
      
      integer                 :: kx,ky,kz, ik, qx,qy,qz, iq, x,y,z ,ix, nnqx,nnqy,nnqz
      complex(kr),allocatable :: realgrid(:,:,:,:)    ! norb,norb,nspin,gw_nkpts
      real(kr)                :: qlength = 1.0_kr
      
      nnqx = nqx+2  ! make space for preceeding 0 and appened datapoint
      nnqy = nqy+2
      nnqz = nqz+2
      allocate( realgrid(norb,norb,nspin, nnqx*nnqy*nnqz) )
      
      ! This function interpolates a matrix in k-space from a coarse
      ! grid to an arbitrary grid
      !
      ! qgrid is the coarse input grid on gw_nkpts points with the
      ! dimensions gw_nkx x gw_nky x gw_nkz taken from the params module
      !
      ! the interpolated result will be saved in kgrid on nkpts points with
      ! the dimensions nkx,nky,nkz
      !
      ! for this a discrete Fourier transform to real space is done on 
      ! the same grid as the coarse one. Then in each dimension a 0 is appended 
      ! to the front of the array and the initial first data point is copied to the
      ! end of the array. This is done for all dimensions. By this, we artifically
      ! can do a Fourier transform on a much larger real space interval where
      ! the function is zero except for small x,y,z where we have the original data.
      ! Nyquist's theorem then allows us to use a finer k-mesh, since 
      ! deltaK ~ 1/len(real_space_interval)
      ! This data is then transformed back and evaluated at an arbitrary finer
      ! k-mesh defined by nkx,nky,nkz
      
      ! Attention: Any modifications here also have to be copied to interpolate_kmatrix_w !!!
 
      ! gw_kvecs(:,:)         ! ik, (kx,ky,kz) are coming from params

      write(*,*) '!!! ============================================================================== !!!'
      write(*,*) 'WARNING: Assuming that the maximum length of each GW k-vector is normalized to 1.0 !!!'
      write(*,*) 'WARNING: For Wien2K it should be changed to 1.0, will fail otherwise               !!!'
      write(*,*) '!!! ============================================================================== !!!'

      realgrid = (0.0_kr,0.0_kr)
      do x=0,nqx-1
        do y=0,nqy-1
           do z=0,nqz-1
              ix = (x+1)*nnqy*nnqz + (y+1)*nnqz + z+1 + 1
           
              ! now do the sum over q-points
              iq = 1
              do qx=0,nqx-1
                 do qy=0,nqy-1
                    do qz=0,nqz-1
                       
                       realgrid(:,:,:,ix) = realgrid(:,:,:,ix) + &
                            &  qgrid(:,:,:,iq)*EXP( 2*pi*ci*( (x-nqx/2)*gw_kvecs(iq,1)/qlength      &
                                                          &   +(y-nqy/2)*gw_kvecs(iq,2)/qlength      &
                                                          &   +(z-nqz/2)*gw_kvecs(iq,3)/qlength ) ) /gw_nkpts
                       iq = iq+1
                    enddo
                 enddo
              enddo
              
           enddo ! z loop
        enddo ! y loop
      enddo ! x loop
      
      !copy boundary data from the front to the back of the array
      ! z plane
      do x=0,nnqx-1
        do y=0,nnqy-1
           realgrid(:,:,:,x*nnqy*nnqz + y*nnqz + (nnqz-1) +1) = realgrid(:,:,:,x*nnqy*nnqz + y*nnqz + 1 +1)
        enddo
      enddo
      ! y plane
      do x=0,nnqx-1
        do z=0,nnqz-1
           realgrid(:,:,:,x*nnqy*nnqz + (nnqy-1)*nnqz + z +1) = realgrid(:,:,:,x*nnqy*nnqz + 1*nnqz + z +1)
        enddo
      enddo
      ! x plane
      do z=0,nnqz-1
        do y=0,nnqy-1
           realgrid(:,:,:,(nnqx-1)*nnqy*nnqz + y*nnqz + z +1) = realgrid(:,:,:,1*nnqy*nnqz + y*nnqz + z +1)
        enddo
      enddo
      
     ! write(*,*) realgrid(1,1,1,:)
      
      ! now we have the real space data
      ! go back to finer k-space
      kgrid = (0.0_kr,0.0_kr)
      ik = 1
      do kx=0,nkx-1
        do ky=0,nky-1
           do kz=0,nkz-1
           
              ! now do the sum over realspace-points
              ix = 1
              do x=0,nnqx-1
                 do y=0,nnqy-1
                    do z=0,nnqz-1
                       kgrid(:,:,:,ik) = kgrid(:,:,:,ik) + &
                         &  realgrid(:,:,:,ix)*EXP(-2*pi*ci*( (x-nnqx/2)*kx*1.0/nkx      &
                                                          &   +(y-nnqy/2)*ky*1.0/nky      &
                                                          &   +(z-nnqz/2)*kz*1.0/nkz ) ) 
                       ix = ix+1
                    enddo
                 enddo
              enddo
              
              ik = ik+1
           enddo
        enddo
      enddo      
      
      deallocate( realgrid )
      
      !write(*,*) ' '
      !write(*,*) gw_kvecs(1,:)
      !write(*,*) qgrid(:,:,1,1)
      !write(*,*) kgrid(:,:,1,1)
      !write(*,*) gw_kvecs(gw_nkpts,:)
      !write(*,*) ' '      
   end subroutine interpolate_kmatrix
   
!  ============================================================
!  == Interpolation in k-space / Wannier-interpolation for Sigma_GW (copy of interpolate_kmatrix)
!  ============================================================
   subroutine interpolate_kmatrix_w(kgrid,qgrid)
      use params
      complex(kr),intent(inout) :: kgrid(:,:,:,:,:)       ! norb,norb,nspin,nkpts,nomega
      complex(kr),intent(in)    :: qgrid(:,:,:,:,:)       ! norb,norb,nspin,gw_nkpts,nomega
      
      integer                 :: kx,ky,kz, ik, qx,qy,qz, iq, x,y,z ,ix, nnqx,nnqy,nnqz
      complex(kr),allocatable :: realgrid(:,:,:,:,:)    ! norb,norb,nspin,gw_nkpts,nomega
      real(kr)                :: qlength = 1.0_kr

      write(*,*) '!!! ============================================================================== !!!'
      write(*,*) 'WARNING: Assuming that the maximum length of each GW k-vector is normalized to 1.0 !!!'
      write(*,*) 'WARNING: For Wien2K it should be changed to 1.0, will fail otherwise               !!!'
      write(*,*) '!!! ============================================================================== !!!'


      nnqx = nqx+2  ! make space for preceeding 0 and appened datapoint
      nnqy = nqy+2
      nnqz = nqz+2
      allocate( realgrid(norb,norb,nspin,nnqx*nnqy*nnqz,nomega) )
 
      ! This function is an exact copy of interpolate_kmatrix only with another dimension
      ! for the grid for frequencies! It should be only used for gw_sigma !
      ! If you change something in interpolate_kmatrix, copy all changes to here also!     
      
      realgrid = (0.0_kr,0.0_kr)
      do x=0,nqx-1
        do y=0,nqy-1
           do z=0,nqz-1
              ix = (x+1)*nnqy*nnqz + (y+1)*nnqz + z+1 + 1
           
              ! now do the sum over q-points
              iq = 1
              do qx=0,nqx-1
                 do qy=0,nqy-1
                    do qz=0,nqz-1
                       realgrid(:,:,:,ix,:) = realgrid(:,:,:,ix,:) + &
                            &  qgrid(:,:,:,iq,:)*EXP( 2*pi*ci*( (x-nqx/2)*gw_kvecs(iq,1)/qlength      &
                                                            &   +(y-nqy/2)*gw_kvecs(iq,2)/qlength      &
                                                            &   +(z-nqz/2)*gw_kvecs(iq,3)/qlength ) ) /gw_nkpts
                       iq = iq+1
                    enddo
                 enddo
              enddo
              
           enddo
        enddo
      enddo
      
      !copy boundary data from the front to the back of the array
      ! z plane
      do x=0,nnqx-1
        do y=0,nnqy-1
           realgrid(:,:,:,x*nnqy*nnqz + y*nnqz + (nnqz-1) +1,:) = realgrid(:,:,:,x*nnqy*nnqz + y*nnqz + 1 +1,:)
        enddo
      enddo
      ! y plane
      do x=0,nnqx-1
        do z=0,nnqz-1
           realgrid(:,:,:,x*nnqy*nnqz + (nnqy-1)*nnqz + z +1,:) = realgrid(:,:,:,x*nnqy*nnqz + 1*nnqz + z +1,:)
        enddo
      enddo
      ! x plane
      do z=0,nnqz-1
        do y=0,nnqy-1
           realgrid(:,:,:,(nnqx-1)*nnqy*nnqz + y*nnqz + z +1,:) = realgrid(:,:,:,1*nnqy*nnqz + y*nnqz + z +1,:)
        enddo
      enddo
      
      ! now we have the real space data
      ! go back to finer k-space
      kgrid = (0.0_kr,0.0_kr)
      ik = 1
      do kx=0,nkx-1
        do ky=0,nky-1
           do kz=0,nkz-1
           
              ! now do the sum over realspace-points
              ix = 1
              do x=0,nnqx-1
                 do y=0,nnqy-1
                    do z=0,nnqz-1
                       kgrid(:,:,:,ik,:) = kgrid(:,:,:,ik,:) + &
                         &  realgrid(:,:,:,ix,:)*EXP(-2*pi*ci*( (x-nnqx/2)*kx*1.0/nkx      &
                                                            &   +(y-nnqy/2)*ky*1.0/nky      &
                                                            &   +(z-nnqz/2)*kz*1.0/nkz ) ) 
                       ix = ix+1
                    enddo
                 enddo
              enddo
              
              ik = ik+1
           enddo
        enddo
      enddo      
      
      deallocate( realgrid )
   end subroutine interpolate_kmatrix_w
   
!  ============================================================
!  == Cubic interpolation for kmatrix 
!  ============================================================
   subroutine interpolate_kmatrix_cubic(kgrid,qgrid,inkx,inky,inkz,inqx,inqy,inqz)
      use params
      complex(kr),intent(inout) :: kgrid(:,:,:,:)       ! norb,norb,nspin,nkpts
      complex(kr),intent(in)    :: qgrid(:,:,:,:)       ! norb,norb,nspin,gw_nkpts
      integer(ki),intent(in)    :: inkx,inky,inkz,inqx,inqy,inqz
      
      integer                 :: m1,m2,s,kx,ky,kz,nx,ny,nz, i,j,k
      complex(kr)             :: tmat(4,4), Cmat(4,4), fmat(4,4)
      real(kr)                :: xvec(4),yvec(4),zvec(4),x,y,z
      
      kgrid = (0.0_kr,0.0_kr)
      
      ! first initialize the Cmatrix
      Cmat = transpose(reshape([0.0_kr, 2.0_kr, 0.0_kr, 0.0_kr, &
                            &  -1.0_kr, 0.0_kr, 1.0_kr, 0.0_kr, &
                            &   2.0_kr,-5.0_kr, 4.0_kr,-1.0_kr, &
                            &  -1.0_kr, 3.0_kr,-3.0_kr, 1.0_kr ],[4,4]) )

      do kx=0,inkx-1
         ! get the coordinates for x
         nx = floor(kx*inqx*1.0_kr/inkx)
         x = kx*inqx*1.0_kr/inkx - nx
         xvec = [1.0_kr, x, x*x, x**3]

         do ky=0,inky-1
            ! get the coordinates for y
            ny = floor(ky*inqy*1.0_kr/inky)
            y = ky*inqy*1.0_kr/inky - ny
            yvec = [1.0_kr, y, y*y, y**3]

            do kz=0,inkz-1
               ! get the coordinates for z
               nz = floor(kz*inqz*1.0_kr/inkz)
               z = kz*inqz*1.0_kr/inkz - nz
               zvec = [1.0_kr, z, z*z, z**3]
                
               do s=1,nspin
                  do m1=1,norb
                     do m2=1,norb
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! now do the interpolation at given s,m1,m2
                        tmat = (0.0_kr,0.0_kr)
                         
                        do i=-1,2
                           ! fill fmatrix for given i
                           do j=-1,2
                              do k=-1,2
                                 fmat(k+2,j+2) = qgrid(m1,m2,s, modulo(nx+i+inqx,inqx)*inqy*inqz &
                                                            & + modulo(ny+j+inqy,inqy)*inqz     &
                                                            & + modulo(nz+k+inqz,inqz)   +1    )
                              enddo
                           enddo
                           tmat(i+2,:) = 0.5*matmul( zvec, matmul(Cmat,fmat) )
                        enddo
                        ! now we have the full tmatrix
                        ! so we can get the interpolated value directly by
                        ! x * C * t * C * y

                        kgrid(m1,m2,s, kx*inky*inkz + ky*inkz + kz +1) =                      &
                                      & 0.25*dot_product( xvec , matmul( Cmat, matmul( tmat,      &
                                           matmul( transpose(Cmat), yvec )  ) ) )  
!                        if (kx==2 .and. ky==0 .and. kz==0) then
!                          write(*,*) ,'m1,m2,s',m1,m2,s
!                          write(*,*) qgrid(m1,m2,s, 0*inqy*inqz + 0*inqz + 0 +1)
!                          write(*,*) kgrid(m1,m2,s, 0*inky*inkz + ky*inkz + kz +1)
!                          write(*,*) kgrid(m1,m2,s, 1*inky*inkz + ky*inkz + kz +1)
!                          write(*,*) kgrid(m1,m2,s, 2*inky*inkz + ky*inkz + kz +1)
 !                         write(*,*) qgrid(m1,m2,s, 1*inqy*inqz + 0*inqz + 0 +1)
 !                                                 stop
  !                     endif

                        
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! sanity check: Imaginary part on the diagonal should
                        ! never be positive!
                        if ( m1==m2 .and. aimag(kgrid(m1,m2,s, kx*inky*inkz + ky*inkz + kz +1)) > 0.0_kr ) then
                           kgrid(m1,m2,s, kx*inky*inkz + ky*inkz + kz +1) =     &
                           &      real(kgrid(m1,m2,s, kx*inky*inkz + ky*inkz + kz +1)) - 0.003*ci
                        endif

                     enddo ! m2
                  enddo ! m1
               enddo ! s       

            enddo ! kz
         enddo ! ky
      enddo ! kx 

   end subroutine interpolate_kmatrix_cubic
      
 
 end module fourier_interpolation
