subroutine pwl_basis_1d ( nd, xd, k, ni, xi, bk )

!*****************************************************************************80
!
!! PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, integer ( kind = 4 ) K, the index of the desired basis function,
!    1 <= K <= ND.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) BK(NI), the basis function at the 
!    interpolation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) bk(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)

  bk(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    bk(1:ni) = 1.0D+00
    return
  end if

  do i = 1, ni

    if ( k == 1 .and. xi(i) <= xd(k) ) then

      t = ( xi(i) - xd(k) ) / ( xd(k+1) - xd(k) )
      bk(i) = 1.0D+00 - t

    else if ( k == nd .and. xd(k) <= xi(i) ) then

      t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
      bk(i) = t

    else if ( xd(k-1) < xi(i) .and. xi(i) <= xd(k) ) then

      t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
      bk(i) = t

    else if ( xd(k) <= xi(i) .and. xi(i) < xd(k+1) ) then

      t = ( xi(i) - xd(k) ) / ( xd(k+1) - xd(k) )
      bk(i) = 1.0D+00 - t

    end if

  end do

  return
end
subroutine pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_INTERP_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
  return
end
