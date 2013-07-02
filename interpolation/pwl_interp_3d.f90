subroutine pwl_interp_3d ( nwd, nxd, nyd, wd, xd, yd, zd, ni, wi, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_3D: piecewise linear interpolant to data defined on a 3D grid.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    1 July 2013
!
!  Author:
!
!    David Pugh
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NWD, NXD, NYD, the number of W, X and Y data 
!    values.
!
!    Input, real ( kind = 8 ) WD(NWD), XD(NXD), YD(NYD), the sorted W, X and Y 
!    data.
!
!    Input, real ( kind = 8 ) ZD(NWD,NXD,NYD), the Z data.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) WI(NI), XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nwd
  integer ( kind = 4 ) nxd
  integer ( kind = 4 ) nyd

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) wd(nwd)
  real ( kind = 8 ) wi(ni)
  real ( kind = 8 ) xd(nxd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nyd)
  real ( kind = 8 ) yi(ni)
  real ( kind = 8 ) tmpzi(ni,nwd)
  real ( kind = 8 ) zd(nwd,nxd,nyd)
  real ( kind = 8 ) zi(ni)

  do j = 1, nwd 

    call pwl_interp_2d ( nxd, nyd, xd, yd, zd(j,:,:), ni, xi, yi, tmpzi(:,j) )
    
  end do
  
  do i = 1, ni
  
    call pwl_interp_1d ( nwd, wd, tmpzi(i,:), 1, wi(i), zi(i) )

  end do

  return
end

!  do i = 1, ni
!
!    do j = 1, nwd 
!
!      do k = 1, nxd 
!
!        call pwl_interp_1d ( nyd, yd, zd(j,k,:), 1, yi(i), tmpzi1(k) )
!
!      end do
!
!      call pwl_interp_1d ( nxd, xd, tmpzi1, 1, xi(i), tmpzi2(j) )
!
!    end do
!
!    call pwl_interp_1d ( nwd, wd, tmpzi2, wi(i), zi(i) )
!
!  end do
!
!  return
!end
