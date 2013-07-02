subroutine gamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_VALUES returns some values of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as:
!
!      Gamma(Z) = integral ( 0 <= T < +oo) T^(Z-1) exp(-T) dT
!
!    It satisfies the recursion:
!
!      Gamma(X+1) = X * Gamma(X)
!
!    Gamma is undefined for nonpositive integral X.
!    Gamma(0.5) = sqrt(PI)
!    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.3544907701811032D+01, &
    -0.1005871979644108D+03, &
     0.9943258511915060D+02, &
     0.9513507698668732D+01, &
     0.4590843711998803D+01, &
     0.2218159543757688D+01, &
     0.1772453850905516D+01, &
     0.1489192248812817D+01, &
     0.1164229713725303D+01, &
     0.1000000000000000D+01, &
     0.9513507698668732D+00, &
     0.9181687423997606D+00, &
     0.8974706963062772D+00, &
     0.8872638175030753D+00, &
     0.8862269254527580D+00, &
     0.8935153492876903D+00, &
     0.9086387328532904D+00, &
     0.9313837709802427D+00, &
     0.9617658319073874D+00, &
     0.1000000000000000D+01, &
     0.2000000000000000D+01, &
     0.6000000000000000D+01, &
     0.3628800000000000D+06, &
     0.1216451004088320D+18, &
     0.8841761993739702D+31 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.50D+00, &
    -0.01D+00, &
     0.01D+00, &
     0.10D+00, &
     0.20D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine i4int_to_r8int ( imin, imax, i, rmin, rmax, r )

!*****************************************************************************80
!
!! I4INT_TO_R8INT maps an I4INT to an R8INT.
!
!  Discussion:
!
!    The formula used is:
!
!      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IMIN, IMAX, the range.
!
!    Input, integer ( kind = 4 ) I, the integer to be converted.
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Output, real ( kind = 8 ) R, the corresponding value in [RMIN,RMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin

  if ( imax == imin ) then

    r = 0.5D+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i,        kind = 8 ) * rmin   &
        + real (        i - imin, kind = 8 ) * rmax ) &
        / real ( imax     - imin, kind = 8 )

  end if

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_permute ( n, p, a )

!*****************************************************************************80
!
!! I4VEC_PERMUTE permutes an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_temp
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine legendre_zeros ( n, x )

!*****************************************************************************80
!
!! LEGENDRE_ZEROS computes the zeros of the Legendre polynomial of degree N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2011
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    0 < N.
!
!    Output, real ( kind = 8 ) X(N), the locations of the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtemp

  e1 = real ( n * ( n + 1 ), kind = 8 )

  m = ( n + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * n + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( n, kind = 8 ) ) &
      / real ( 8 * n * n, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, n
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( n, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 ) &
      / ( 1.0D+00 + x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    x(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

  end do

  if ( mod ( n, 2 ) == 1 ) then
    x(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( n + 1 ) / 2
  ncopy = n - nmove

  do i = 1, nmove
    iback = n + 1 - i
    x(iback) = x(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, n - nmove
    x(i) = - x(n+1-i)
  end do

  return
end
subroutine perm_check ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_uniform ( n, base, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input, integer ( kind = 4 ) BASE, is 0 for a 0-based permutation and 1 for
!    a 1-based permutation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  do i = 1, n
    p(i) = ( i - 1 ) + base
  end do

  do i = 1, n
    j = i4_uniform_ab ( i, n, seed )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

  return
end
function r8_abs ( x )

!*****************************************************************************80
!
!! R8_ABS returns the absolute value of an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    FORTRAN90 supplies the ABS function, which should be used instead
!    of this function!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose absolute value is desired.
!
!    Output, real ( kind = 8 ) R8_ABS, the absolute value of X.
!
  implicit none

  real ( kind = 8 ) r8_abs
  real ( kind = 8 ) x

  if ( 0.0D+00 <= x ) then
    r8_abs = + x
  else
    r8_abs = - x
  end if

  return
end
function r8_acos ( c )

!*****************************************************************************80
!
!! R8_ACOS computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) R8_ACOS, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) r8_acos

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  r8_acos = acos ( c2 )

  return
end
function r8_acosh ( x )

!*****************************************************************************80
!
!! R8_ACOSH evaluates the arc-hyperbolic cosine of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_ACOSH, the arc-hyperbolic cosine of X.
!
  implicit none

  real ( kind = 8 ), parameter :: dln2 = 0.69314718055994530941723212145818D+00
  real ( kind = 8 ) r8_acosh
  real ( kind = 8 ) r8_tiny
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00

  if ( xmax == 0.0D+00 ) then
    xmax = 1.0D+00 / sqrt ( r8_tiny ( ) )
  end if

  if ( x < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ACOSH - Fatal error!'
    write ( *, '(a)' ) '  X < 1.0'
    stop
  else if ( x < xmax ) then
    value = log ( x + sqrt ( x * x - 1.0D+00 ) )
  else
    value = dln2 + log ( x )
  end if

  r8_acosh = value

  return
end
function r8_add ( x, y )

!*****************************************************************************80
!
!! R8_ADD returns the sum of two R8's.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    FORTRAN90 supplies the + operator, which should generally be used instead
!    of this function!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the numbers to be added.
!
!    Output, real ( kind = 8 ) R8_ADD, the sum.
!
  implicit none

  real ( kind = 8 ) r8_add
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  r8_add = x + y

  return
end
function r8_aint ( x )

!****************************************************************************80
!
!! R8_AINT truncates an R8 argument to an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2011
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_AINT, the truncated version of X.
!
  implicit none

  real ( kind = 8 ) r8_aint
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = - int ( abs ( x ) )
  else
    value =   int ( abs ( x ) )
  end if

  r8_aint = value

  return
end
function r8_asin ( s )

!*****************************************************************************80
!
!! R8_ASIN computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S, the argument.
!
!    Output, real ( kind = 8 ) R8_ASIN, an angle whose sine is S.
!
  implicit none

  real ( kind = 8 ) r8_asin
  real ( kind = 8 ) s
  real ( kind = 8 ) s2

  s2 = s
  s2 = max ( s2, -1.0D+00 )
  s2 = min ( s2, +1.0D+00 )

  r8_asin = asin ( s2 )

  return
end
function r8_atan ( y, x )

!*****************************************************************************80
!
!! R8_ATAN computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) R8_ATAN, an angle between 0 and 2 * PI, whose
!    tangent is (Y/X), and which lies in the appropriate quadrant so that
!    the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 8 ) abs_x
  real ( kind = 8 ) abs_y
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_atan
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_0
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  r8_atan = theta

  return
end
function r8_cas ( x )

!*****************************************************************************80
!
!! R8_CAS returns the "casine" of an R8.
!
!  Discussion:
!
!    The "casine", used in the discrete Hartley transform, is abbreviated
!    CAS(X), and defined by:
!
!      CAS(X) = cos ( X ) + sin( X )
!             = sqrt ( 2 ) * sin ( X + pi/4 )
!             = sqrt ( 2 ) * cos ( X - pi/4 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose casine is desired.
!
!    Output, real ( kind = 8 ) R8_CAS, the casine of X, which will be between
!    plus or minus the square root of 2.
!
  implicit none

  real ( kind = 8 ) r8_cas
  real ( kind = 8 ) x

  r8_cas = cos ( x ) + sin ( x )

  return
end
function r8_ceiling ( r )

!*****************************************************************************80
!
!! R8_CEILING rounds an R8 "up" (towards +oo) to an integral R8.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1.0
!   -1.0  -1.0
!   -0.9   0.0
!    0.0   0.0
!    5.0   5.0
!    5.1   6.0
!    5.9   6.0
!    6.0   6.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded up.
!
!    Output, real ( kind = 8 ) R8_CEILING, the rounded value.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) r8_ceiling
  integer ( kind = 4 ) value

  value = real ( int ( r ), kind = 8 )
  if ( value < r ) then
    value = value + 1.0D+00
  end if

  r8_ceiling = value

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, real ( kind = 8 ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 8 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 8 ) ) / real ( i, kind = 8 )
    end do

  end if

  r8_choose = value

  return
end
function r8_chop ( place, x )

!*****************************************************************************80
!
!! R8_CHOP chops an R8 to a given number of binary places.
!
!  Example:
!
!    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
!
!    The following values would be returned for the 'chopped' value of
!    3.875:
!
!    PLACE  Value
!
!       1      2
!       2      3     = 2 + 1
!       3      3.5   = 2 + 1 + 1/2
!       4      3.75  = 2 + 1 + 1/2 + 1/4
!       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLACE, the number of binary places to preserve.
!    PLACE = 0 means return the integer part of X.
!    PLACE = 1 means return the value of X, correct to 1/2.
!    PLACE = 2 means return the value of X, correct to 1/4.
!    PLACE = -1 means return the value of X, correct to 2.
!
!    Input, real ( kind = 8 ) X, the number to be chopped.
!
!    Output, real ( kind = 8 ) R8_CHOP, the chopped number.
!
  implicit none

  real ( kind = 8 ) fac
  integer ( kind = 4 ) place
  real ( kind = 8 ) r8_chop
  real ( kind = 8 ) r8_log_2
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) s
  integer ( kind = 4 ) temp
  real ( kind = 8 ) x

  s = r8_sign ( x )
  temp = int ( r8_log_2 ( abs ( x ) ) )
  fac = 2.0D+00**( temp - place + 1 )
  r8_chop = s * real ( int ( abs ( x ) / fac ), kind = 8 ) * fac

  return
end
function r8_csc ( theta )

!*****************************************************************************80
!
!! R8_CSC returns the cosecant of X.
!
!  Discussion:
!
!    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
!
!    The cosecant is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) THETA, the angle, in radians, whose
!    cosecant is desired.  It must be the case that SIN ( THETA ) is not zero.
!
!    Output, real ( kind = 8 ) R8_CSC, the cosecant of THETA.
!
  implicit none

  real ( kind = 8 ) r8_csc
  real ( kind = 8 ) theta
  real ( kind = 8 ) value

  value = sin ( theta )

  if ( value == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSC - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Cosecant undefined for THETA = ', theta
    stop
  end if

  r8_csc = 1.0D+00 / value

  return
end
function r8_csqrt ( x )

!*****************************************************************************80
!
!! R8_CSQRT returns the complex square root of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Output, complex ( kind = 8 ) R8_CSQRT, the square root of X:
!
  implicit none

  real ( kind = 8 ) argument
  real ( kind = 8 ) magnitude
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  complex ( kind = 8 ) r8_csqrt
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    magnitude = x
    argument = 0.0D+00
  else if ( 0.0D+00 == x ) then
    magnitude = 0.0D+00
    argument = 0.0D+00
  else if ( x < 0.0D+00 ) then
    magnitude = -x
    argument = pi
  end if

  magnitude = sqrt ( magnitude )
  argument = argument / 2.0D+00

  r8_csqrt = magnitude * cmplx ( cos ( argument ), sin ( argument ), kind = 8 )

  return
end
function r8_cube_root ( x )

!*****************************************************************************80
!
!! R8_CUBE_ROOT returns the cube root of an R8.
!
!  Discussion:
!
!    This routine is designed to avoid the possible problems that can occur
!    when formulas like 0.0^(1/3) or (-1.0)^(1/3) are to be evaluated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose cube root is desired.
!
!    Output, real ( kind = 8 ) R8_CUBE_ROOT, the cube root of X.
!
  implicit none

  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    value = x ** ( 1.0D+00 / 3.0D+00 )
  else if ( x == 0.0D+00 ) then
    value = 0.0D+00
  else
    value = -( abs ( x ) ) ** ( 1.0D+00 / 3.0D+00 )
  end if

  r8_cube_root = value

  return
end
function r8_degrees ( radians )

!*****************************************************************************80
!
!! R8_DEGREES converts an angle from radian to degree measure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIANS, the angle measurement in radians.
!
!    Output, real ( kind = 8 ) R8_DEGREES, the angle measurement in degrees.
!
  implicit none

  real ( kind = 8 ) r8_degrees
  real ( kind = 8 ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ) radians

  r8_degrees = radians * 180.0D+00 / r8_pi

  return
end
function r8_diff ( x, y, n )

!*****************************************************************************80
!
!! R8_DIFF computes the difference of two R8's to a specified accuracy.
!
!  Discussion:
!
!    The user controls how many binary digits of accuracy
!    are to be used.
!
!    N determines the accuracy of the value of the result.  If N = 10,
!    for example, only 11 binary places will be used in the arithmetic.
!    In general, only N+1 binary places will be used.
!
!    N may be zero.  However, a negative value of N should
!    not be used, since this will cause both X and Y to look like 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the two values whose difference is desired.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to use.
!
!    Output, real ( kind = 8 ) R8_DIFF, the value of X-Y.
!
  implicit none

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  integer ( kind = 4 ) n
  real ( kind = 8 ) pow2
  real ( kind = 8 ) r8_diff
  real ( kind = 8 ) size
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x == y ) then
    r8_diff = 0.0D+00
    return
  end if

  pow2 = 2.0D+00**n
!
!  Compute the magnitude of X and Y, and take the larger of the
!  two.  At least one of the two values is not zero!
!
  size = max ( abs ( x ), abs ( y ) )
!
!  Make normalized copies of X and Y.  One of the two values will
!  actually be equal to 1.
!
  cx = x / size
  cy = y / size
!
!  Here's where rounding comes in.  We know that the larger of the
!  the two values equals 1.  We multiply both values by 2^N,
!  where N+1 is the number of binary digits of accuracy we want
!  to use, truncate the values, and divide back by 2^N.
!
  cx = real ( int ( cx * pow2 + sign ( 0.5D+00, cx ) ), kind = 8 ) / pow2
  cy = real ( int ( cy * pow2 + sign ( 0.5D+00, cy ) ), kind = 8 ) / pow2
!
!  Take the difference now.
!
  r8_diff = cx - cy
!
!  Undo the scaling.
!
  r8_diff = r8_diff * size

  return
end
subroutine r8_digit ( x, idigit, digit )

!*****************************************************************************80
!
!! R8_DIGIT returns a particular decimal digit of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose NDIG-th decimal digit
!    is desired.  If X is zero, all digits will be returned as 0.
!
!    Input, integer ( kind = 4 ) IDIGIT, the position of the desired decimal
!    digit.  A value of 1 means the leading digit, a value of 2 the second digit
!    and so on.
!
!    Output, integer ( kind = 4 ) DIGIT, the value of the IDIGIT-th decimal
!    digit of X.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ival
  real ( kind = 8 ) x
  real ( kind = 8 ) xcopy

  if ( x == 0.0D+00 ) then
    digit = 0
    return
  end if

  if ( idigit <= 0 ) then
    digit = 0
    return
  end if
!
!  Set XCOPY = X, and then force XCOPY to lie between 1 and 10.
!
  xcopy = abs ( x )

  do while ( xcopy < 1.0D+00 )
    xcopy = xcopy * 10.0D+00
  end do

  do while ( 10.0D+00 <= xcopy )
    xcopy = xcopy / 10.0D+00
  end do

  do i = 1, idigit
    ival = int ( xcopy )
    xcopy = ( xcopy - ival ) * 10.0D+00
  end do

  digit = ival

  return
end
function r8_divide_i4 ( i, j )

!*****************************************************************************80
!
!! R8_DIVIDE_I4 returns an I4 fraction as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the numerator and denominator.
!
!    Output, real ( kind = 8 ) R8_DIVIDE_I4, the value of (I/J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_divide_i4

  r8_divide_i4 = real ( i, kind = 8 ) / real ( j, kind = 8 )

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the round-off unit.
!
  implicit none

  real ( kind = 8 ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end
function r8_epsilon_compute ( )

!*****************************************************************************80
!
!! R8_EPSILON_COMPUTE computes the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON_COMPUTE, the computed round-off unit.
!
  implicit none

  real ( kind = 8 ) one
  real ( kind = 8 ) r8_add
  real ( kind = 8 ) r8_epsilon_compute
  real ( kind = 8 ) temp
  real ( kind = 8 ) test
  real ( kind = 8 ) value

  one = real ( 1, kind = 8 )

  value = one
  temp = value / 2.0D+00
  test = r8_add ( one, temp )

  do while ( one < test )
    value = temp
    temp = value / 2.0D+00
    test = r8_add ( one, temp )
  end do

  r8_epsilon_compute = value

  return
end
function r8_exp ( x )

!*****************************************************************************80
!
!! R8_EXP computes the exponential of an R8, avoiding overflow and underflow.
!
!  Discussion:
!
!    My experience with the G95 compiler has included many unpleasant
!    floating point exceptions when very small arguments are given to
!    the exponential function.
!
!    This routine is designed to avoid such problems.
!
!    Ideally, the rule would be:
!
!                    X <= log ( TINY ) => R8_EXP ( X ) = 0
!    log ( HUGE ) <= X                 => R8_EXP ( X ) = HUGE
!
!    However, the G95 math library seems to produce infinity for
!    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
!    included a fudge factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the exponential function.
!
!    Output, real ( kind = 8 ) R8_EXP, the value of exp ( X ).
!
  implicit none

  real ( kind = 8 ), parameter :: log_max = 709.711D+00
  real ( kind = 8 ), parameter :: log_min = -708.467D+00
  real ( kind = 8 ) r8_exp
  real ( kind = 8 ) x

  if ( x <= log_min ) then
    r8_exp = 0.0D+00
  else if ( x < log_max ) then
    r8_exp = exp ( x )
  else
    r8_exp = huge ( x )
  end if

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1, the value is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = 8 )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

  return
end
function r8_floor ( r )

!*****************************************************************************80
!
!! R8_FLOOR rounds an R8 "down" (towards -oo) to the nearest integral R8.
!
!  Example:
!
!    R     Value
!
!   -1.1  -2.0
!   -1.0  -1.0
!   -0.9  -1.0
!    0.0   0.0
!    5.0   5.0
!    5.1   5.0
!    5.9   5.0
!    6.0   6.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded down.
!
!    Output, real ( kind = 8 ) R8_FLOOR, the rounded value.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) r8_floor
  real ( kind = 8 ) value

  value = real ( int ( r ), kind = 8 )
  if ( r < value ) then
    value = value - 1.0D+00
  end if

  r8_floor = value

  return
end
function r8_fraction ( i, j )

!*****************************************************************************80
!
!! R8_FRACTION uses real arithmetic on an integer ratio.
!
!  Discussion:
!
!    Given integer variables I and J, both FORTRAN and C will evaluate
!    an expression such as "I/J" using what is called "integer division",
!    with the result being an integer.  It is often convenient to express
!    the parts of a fraction as integers but expect the result to be computed
!    using real arithmetic.  This function carries out that operation.
!
!  Example:
!
!       I     J   I/J  R8_FRACTION
!
!       1     2     0  0.5
!       7     4     1  1.75
!       8     4     2  2.00
!       9     4     2  2.25
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the arguments.
!
!    Output, real ( kind = 8 ) R8_FRACTION, the value of the ratio.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_fraction

  r8_fraction = real ( i, kind = 8 ) / real ( j, kind = 8 )

  return
end
function r8_fractional ( x )

!*****************************************************************************80
!
!! R8_FRACTIONAL returns the fractional part of an R8.
!
!  Discussion:
!
!    If we regard a real number as
!
!      R = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R      FRACTION
!
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_FRACTIONAL, the fractional part of X.
!
  implicit none

  real ( kind = 8 ) r8_fractional
  real ( kind = 8 ) x

  r8_fractional = abs ( x ) - real ( int ( abs ( x ) ), kind = 8 )

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.79D+308
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < r8_epsilon ( ) ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ) :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ) :: d2 = 4.227843350984671393993777D-01
  real ( kind = 8 ) :: d4 = 1.791759469228055000094023D+00
  real ( kind = 8 ), parameter :: frtbig = 2.25D+76
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 2.55D+305
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.79D+308
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq

  y = x

  if ( 0.0D+00 < y .and. y <= xbig ) then

    if ( y <= epsilon ( y ) ) then

      res = - log ( y )
!
!  EPS < X <= 1.5.
!
    else if ( y <= 1.5D+00 ) then

      if ( y < 0.6796875D+00 ) then
        corr = -log ( y )
        xm1 = y
      else
        corr = 0.0D+00
        xm1 = ( y - 0.5D+00 ) - 0.5D+00
      end if

      if ( y <= 0.5D+00 .or. 0.6796875D+00 <= y ) then

        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm1 + p1(i)
          xden = xden * xm1 + q1(i)
        end do

        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

      else

        xm2 = ( y - 0.5D+00 ) - 0.5D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

      end if
!
!  1.5 < X <= 4.0.
!
    else if ( y <= 4.0D+00 ) then

      xm2 = y - 2.0D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0 < X <= 12.0.
!
    else if ( y <= 12.0D+00 ) then

      xm4 = y - 4.0D+00
      xden = -1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
      end do

      res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
    else

      res = 0.0D+00

      if ( y <= frtbig ) then

        res = c(7)
        ysq = y * y

        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if

      res = res / y
      corr = log ( y )
      res = res + sqrtpi - 0.5D+00 * corr
      res = res + y * ( corr - 1.0D+00 )

    end if
!
!  Return for bad arguments.
!
  else

    res = xinf

  end if
!
!  Final adjustments and return.
!
  r8_gamma_log = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_hypot ( x, y )

!*****************************************************************************80
!
!! R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_hypot
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( abs ( x ) < abs ( y ) ) then
    a = abs ( y )
    b = abs ( x )
  else
    a = abs ( x )
    b = abs ( y )
  end if
!
!  A contains the larger value.
!
  if ( a == 0.0D+00 ) then
    c = 0.0D+00
  else
    c = a * sqrt ( 1.0D+00 + ( b / a )**2 )
  end if

  r8_hypot = c

  return
end
function r8_in_01 ( a )

!*****************************************************************************80
!
!! R8_IN_01 is TRUE if an R8 is in the range [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the value.
!
!    Output, logical R8_IN_01, is TRUE if 0 <= A <= 1.
!
  implicit none

  real ( kind = 8 ) a
  logical r8_in_01
  logical value

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    value = .false.
  else
    value = .true.
  end if

  r8_in_01 = value

  return
end
function r8_insignificant ( r, s )

!*****************************************************************************80
!
!! R8_INSIGNIFICANT determines if an R8 is insignificant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be compared against.
!
!    Input, real ( kind = 8 ) S, the number to be compared.
!
!    Output, logical R8_INSIGNIFICANT, is TRUE if S is insignificant
!    compared to R.
!
  implicit none

  real ( kind = 8 ) r
  logical r8_insignificant
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  logical value

  value = .true. 

  t = r + s
  tol = epsilon ( r ) * abs ( r )

  if ( tol < abs ( r - t ) ) then 
    value = .false.
  end if
  
  r8_insignificant = value

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if an R8 represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  real ( kind = 8 ) r
  logical r8_is_int
  logical value

  if ( real ( i4_huge, kind = 8 ) < r ) then
    value = .false.
  else if ( r < - real ( i4_huge, kind = 8 ) ) then
    value = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    value = .true.
  else
    value = .false.
  end if

  r8_is_int = value

  return
end
function r8_log_2 ( x )

!*****************************************************************************80
!
!! R8_LOG_2 returns the logarithm base 2 of an R8.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2^R8_LOG_2.
!
  implicit none

  real ( kind = 8 ) r8_log_2
  real ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_2 = - huge ( x )
  else
    r8_log_2 = log ( abs ( x ) ) / log ( 2.0D+00 )
  end if

  return
end
function r8_log_10 ( x )

!*****************************************************************************80
!
!! R8_LOG_10 returns the logarithm base 10 of an R8.
!
!  Discussion:
!
!    value = Log10 ( |X| )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_10, the logarithm base 10 of the absolute
!    value of X.  It should be true that |X| = 10**R8_LOG_10.
!
  implicit none

  real ( kind = 8 ) r8_log_10
  real ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_10 = - huge ( x )
  else
    r8_log_10 = log10 ( abs ( x ) )
  end if

  return
end
function r8_log_b ( x, b )

!*****************************************************************************80
!
!! R8_LOG_B returns the logarithm base B of an R8.
!
!  Discussion:
!
!    value = log ( |X| ) / log ( |B| )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base B logarithm is desired.
!    X should not be 0.
!
!    Input, real ( kind = 8 ) B, the base, which should not be 0, 1 or -1.
!
!    Output, real ( kind = 8 ) R8_LOG_B, the logarithm base B of the absolute
!    value of X.  It should be true that |X| = |B|**R8_LOG_B.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) r8_log_b
  real ( kind = 8 ) x

  if ( b == 0.0D+00 .or. b == 1.0D+00 .or. b == - 1.0D+00 ) then
    r8_log_b = - huge ( x )
  else if ( abs ( x ) == 0.0D+00 ) then
    r8_log_b = - huge ( x )
  else
    r8_log_b = log ( abs ( x ) ) / log ( abs ( b ) )
  end if

  return
end
subroutine r8_mant ( x, s, r, l )

!*****************************************************************************80
!
!! R8_MANT computes the "mantissa" or "fraction part" of an R8.
!
!  Discussion:
!
!    X = S * R * 2^L
!
!    S is +1 or -1,
!    R is an real value between 1.0 and 2.0,
!    L is an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, integer ( kind = 4 ) S, the "sign" of the number.
!    S will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, real ( kind = 8 ) R, the mantissa of X.  R will be greater
!    than or equal to 1, and strictly less than 2.  The one
!    exception occurs if X is zero, in which case R will also
!    be zero.
!
!    Output, integer ( kind = 4 ) L, the integer part of the logarithm
!    (base 2) of X.
!
  implicit none

  integer ( kind = 4 ) l
  real ( kind = 8 ) r
  integer ( kind = 4 ) s
  real ( kind = 8 ) x
!
!  Determine the sign.
!
  if ( x < 0.0D+00 ) then
    s = -1
  else
    s = + 1
  end if
!
!  Set R to the absolute value of X, and L to zero.
!  Then force R to lie between 1 and 2.
!
  if ( x < 0.0D+00 ) then
    r = - x
  else
    r = + x
  end if

  l = 0
!
!  Time to bail out if X is zero.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  do while ( 2.0D+00 <= r )
    r = r / 2.0D+00
    l = l + 1
  end do

  do while ( r < 1.0D+00 )
    r = r * 2.0D+00
    l = l - 1
  end do

  return
end
function r8_mod ( x, y )

!*****************************************************************************80
!
!! R8_MOD returns the remainder of R8 division.
!
!  Discussion:
!
!    If
!      REM = R8_MOD ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM has the same sign as X, and abs ( REM ) < Y.
!
!  Example:
!
!        X         Y     R8_MOD  R8_MOD Factorization
!
!      107        50       7      107 =  2 *  50 + 7
!      107       -50       7      107 = -2 * -50 + 7
!     -107        50      -7     -107 = -2 *  50 - 7
!     -107       -50      -7     -107 =  2 * -50 - 7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MOD, the remainder when X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_mod
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MOD - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MOD ( X, Y ) called with Y = ', y
    stop
  end if

  r8_mod = x - real ( int ( x / y ), kind = 8 ) * y

  if ( x < 0.0D+00 .and. 0.0D+00 < r8_mod ) then
    r8_mod = r8_mod - abs ( y )
  else if ( 0.0D+00 < x .and. r8_mod < 0.0D+00 ) then
    r8_mod = r8_mod + abs ( y )
  end if

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of R8 division.
!
!  Discussion:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        X         Y     MOD R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder
!    when X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 8 ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
function r8_nint ( x )

!*****************************************************************************80
!
!! R8_NINT returns the nearest integer to an R8.
!
!  Example:
!
!        X        R8_NINT
!
!      1.3         1
!      1.4         1
!      1.5         1 or 2
!      1.6         2
!      0.0         0
!     -0.7        -1
!     -1.1        -1
!     -1.6        -2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, integer ( kind = 4 ) R8_NINT, the nearest integer to X.
!
  implicit none

  integer ( kind = 4 ) r8_nint
  integer ( kind = 4 ) s
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    s = - 1
  else
    s = + 1
  end if

  r8_nint = s * int ( abs ( x ) + 0.5D+00 )

  return
end
function r8_normal ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( - 2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal = a + b * x

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( - 2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI, the value of pi.
!
  implicit none

  real ( kind = 8 ) r8_pi

  r8_pi = 3.141592653589793D+00

  return
end
function r8_pi_sqrt ( )

!*****************************************************************************80
!
!! R8_PI_SQRT returns the square root of pi as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI_SQRT, the square root of pi.
!
  implicit none

  real ( kind = 8 ) r8_pi_sqrt

  r8_pi_sqrt = 1.7724538509055160273D+00

  return
end
function r8_power ( r, p )

!*****************************************************************************80
!
!! R8_POWER computes the P-th power of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 8 ) R8_POWER, the value of the P-th power of R.
!
  implicit none

  integer ( kind = 4 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_power
  real ( kind = 8 ) value
!
!  Special case.  R^0 = 1.
!
  if ( p == 0 ) then

    value = 1.0D+00
!
!  Special case.  Positive powers of 0 are 0.
!  For negative powers of 0, we go ahead and compute R^P,
!  relying on the software to complain.
!
  else if ( r == 0.0D+00 ) then

    if ( 0 < p ) then
      value = 0.0D+00
    else
      value = r**p
    end if

  else if ( 1 <= p ) then
    value = r**p
  else
    value = 1.0D+00 / r**(-p)
  end if

  r8_power = value

  return
end
subroutine r8_power_fast ( r, p, rp, mults )

!*****************************************************************************80
!
!! R8_POWER_FAST computes an integer power of an R8.
!
!  Discussion:
!
!    Obviously, R^P can be computed using P-1 multiplications.
!
!    However, R^P can also be computed using at most 2*LOG2(P) multiplications.
!    To do the calculation this way, let N = LOG2(P).
!    Compute A, A^2, A^4, ..., A^N by N-1 successive squarings.
!    Start the value of R^P at A, and each time that there is a 1 in
!    the binary expansion of P, multiply by the current result of the squarings.
!
!    This algorithm is not optimal.  For small exponents, and for special
!    cases, the result can be computed even more quickly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 8 ) RP, the value of R^P.
!
!    Output, integer ( kind = 4 ) MULTS, the number of multiplications
!    and divisions.
!
  implicit none

  integer ( kind = 4 ) mults
  integer ( kind = 4 ) p
  integer ( kind = 4 ) p_mag
  integer ( kind = 4 ) p_sign
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) rp

  mults = 0
!
!  Special bases.
!
  if ( r == 1.0D+00 ) then
    rp = 1.0D+00
    return
  end if

  if ( r == -1.0D+00 ) then

    if ( mod ( p, 2 ) == 1 ) then
      rp = -1.0D+00
    else
      rp = 1.0D+00
    end if

    return

  end if

  if ( r == 0.0D+00 ) then

    if ( p <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_POWER_FAST - Fatal error!'
      write ( *, '(a)' ) '  Base R is zero, and exponent is negative.'
      write ( *, '(a,i8)' ) '  Exponent P = ', p
      stop
    end if

    rp = 0.0D+00
    return

  end if
!
!  Special powers.
!
  if ( p == -1 ) then
    rp = 1.0D+00 / r
    mults = mults + 1
    return
  else if ( p == 0 ) then
    rp = 1.0D+00
    return
  else if ( p == 1 ) then
    rp = r
    return
  end if
!
!  Some work to do.
!
  p_mag = abs ( p )
  p_sign = sign ( 1, p )

  rp = 1.0D+00
  r2 = r

  do while ( 0 < p_mag )

    if ( mod ( p_mag, 2 ) == 1 ) then
      rp = rp * r2
      mults = mults + 1
    end if

    p_mag = p_mag / 2
    r2 = r2 * r2
    mults = mults + 1

  end do

  if ( p_sign == -1 ) then
    rp = 1.0D+00 / rp
    mults = mults + 1
  end if

  return
end
function r8_pythag ( a, b )

!*****************************************************************************80
!
!! R8_PYTHAG computes sqrt ( A * A + B * B ), avoiding overflow and underflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the values for which sqrt ( A * A + B * B )
!    is desired.
!
!    Output, real ( kind = 8 ) R8_PYTHAG, the value of sqrt ( A * A + B * B ).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_abs
  real ( kind = 8 ) b
  real ( kind = 8 ) b_abs
  real ( kind = 8 ) r8_pythag

  a_abs = abs ( a )
  b_abs = abs ( b )

  if ( b_abs < a_abs ) then
    r8_pythag = a_abs * sqrt ( 1.0D+00 + ( b_abs / a_abs ) * ( b_abs / a_abs ) )
  else if ( b_abs == 0.0D+00 ) then
    r8_pythag = 0.0D+00
  else if ( a_abs <= b_abs ) then
    r8_pythag = b_abs * sqrt ( 1.0D+00 + ( a_abs / b_abs ) * ( a_abs / b_abs ) )
  end if

  return
end
function r8_radians ( degrees )

!*****************************************************************************80
!
!! R8_RADIANS converts an angle from degree to radian measure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DEGREES, the angle measurement in degrees.
!
!    Output, real ( kind = 8 ) R8_RADIANS, the angle measurement in radians.
!
  implicit none

  real ( kind = 8 ) degrees
  real ( kind = 8 ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ) r8_radians

  r8_radians = degrees * r8_pi / 180.0D+00

  return
end
function r8_round ( x )

!*****************************************************************************80
!
!! R8_ROUND sets an R8 to the nearest integral value.
!
!  Example:
!
!        X        R8_ROUND
!
!      1.3         1.0
!      1.4         1.0
!      1.5         1.0 or 2.0
!      1.6         2.0
!      0.0         0.0
!     -0.7        -1.0
!     -1.1        -1.0
!     -1.6        -2.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, real ( kind = 8 ) R8_ROUND, the rounded value.
!
  implicit none

  real ( kind = 8 ) r8_round
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = - real ( int ( - x + 0.5D+00 ), kind = 8 )
  else
    value =   real ( int ( + x + 0.5D+00 ), kind = 8 )
  end if

  r8_round = value

  return
end
function r8_round_i4 ( x )

!*****************************************************************************80
!
!! R8_ROUND_I4 sets an R8 to the nearest integral value, returning an I4
!
!  Example:
!
!        X        R8_ROUND_I4
!
!      1.3         1
!      1.4         1
!      1.5         1 or 2
!      1.6         2
!      0.0         0
!     -0.7        -1
!     -1.1        -1
!     -1.6        -2
!
!  Discussion:
!
!    In FORTRAN90, we rely on the fact that, for positive X, int ( X )
!    is the "floor" function, returning the largest integer less than
!    or equal to X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, integer ( kind = 4 ) R8_ROUND_I4, the rounded value.
!
  implicit none

  integer ( kind = 4 ) r8_round_i4
  integer ( kind = 4 ) value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = - int ( - x + 0.5D+00 )
  else
    value =   int ( + x + 0.5D+00 )
  end if

  r8_round_i4 = value

  return
end
subroutine r8_round2 ( nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUND2 rounds an R8 in base 2.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 2^L
!
!    where S is plus or minus 1, L is an integer, and J is a binary
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.5 and strictly less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 2^L
!
!    where S and L are unchanged, and K is a binary mantissa which
!    agrees with J in the first NPLACE binary digits and is zero
!    thereafter.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
!    or 0.75.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of binary digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  integer ( kind = 4 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) xmant
  real ( kind = 8 ) xround
  real ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign S.
!
  if ( 0.0D+00 < x ) then
    s = 1
    xtemp = x
  else
    s = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the
!  logarithm L.
!
  l = 0

  do while ( 2.0D+00 <= xtemp )
    xtemp = xtemp / 2.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 2.0D+00
    l = l - 1
  end do
!
!  4: Strip out the digits of the mantissa as XMANT, and decrease L.
!
  xmant = 0.0D+00
  iplace = 0

  do

    xmant = 2.0D+00 * xmant

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + 1.0D+00
      xtemp = xtemp - 1.0D+00
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = s * xmant * 2.0D+00**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * 2.0D+00

  end do

  return
end
subroutine r8_roundb ( base, nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUNDB rounds an R8 in a given base.
!
!  Discussion:
!
!    The code does not seem to do a good job of rounding when
!    the base is negative.
!
!    Assume that the input quantity X has the form
!
!      X = S * J * BASE^L
!
!    where S is plus or minus 1, L is an integer, and J is a
!    mantissa base BASE which is either exactly zero, or greater
!    than or equal to (1/BASE) and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * BASE^L
!
!    where S and L are unchanged, and K is a mantissa base BASE
!    which agrees with J in the first NPLACE digits and is zero
!    thereafter.
!
!    Note that because of rounding, for most bases, most numbers
!    with a fractional quantities cannot be stored exactly in the
!    computer, and hence will have trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0,
!    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0,
!    BASE/BASE^2, (BASE+1)/BASE^2, ...,
!    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base of the arithmetic.
!    BASE must not be zero.  Theoretically, BASE may be negative.
!
!    Input, integer ( kind = 4 ) NPLACE, the number of digits base BASE to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) js
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real ( kind = 8 ) x
  real ( kind = 8 ) xmant
  real ( kind = 8 ) xround
  real ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  0: Error checks.
!
  if ( base == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ROUNDB - Fatal error!'
    write ( *, '(a)' ) '  The base BASE cannot be zero.'
    stop
  end if
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0D+00 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
!  logarithm L.
!
  l = 0

  do while ( abs ( base ) <= abs ( xtemp ) )

    xtemp = xtemp / real ( base, kind = 8 )

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l + 1

  end do

  do while ( abs ( xtemp ) < 1.0D+00 )

    xtemp = xtemp * base

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l - 1

  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0D+00
  iplace = 0
  js = is

  do

    xmant = base * xmant

    if ( xmant < 0.0D+00 ) then
      js = -js
      xmant = -xmant
    end if

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = js * xmant * ( real ( base, kind = 8 ) )**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * base

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

  end do

  return
end
subroutine r8_roundx ( nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUNDX rounds an R8 in base 10.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 10^L
!
!    where S is plus or minus 1, L is an integer, and J is a decimal
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.1 and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 10^L
!
!    where S and L are unchanged, and K is a decimal mantissa which
!    agrees with J in the first NPLACE decimal digits and is zero
!    thereafter.
!
!    Note that because of rounding, most decimal fraction quantities
!    cannot be stored exactly in the computer, and hence will have
!    trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
!    0.2, ..., or 0.9.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
!    0.03, ..., 0.98, 0.99.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of decimal digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real ( kind = 8 ) x
  real ( kind = 8 ) xmant
  real ( kind = 8 ) xround
  real ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0D+00 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 10, and compute the
!  logarithm L.
!
  l = 0

  do while ( 10.0D+00 <= x )
    xtemp = xtemp / 10.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 10.0D+00
    l = l - 1
  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0D+00
  iplace = 0

  do

    xmant = 10.0D+00 * xmant

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = is * xmant * ( 10.0D+00**l )
      exit
    end if

    l = l - 1
    xtemp = xtemp * 10.0D+00

  end do

  return
end
function r8_sech ( x )

!*****************************************************************************80
!
!! R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_SECH, the value of the function.
!
  implicit none

  real ( kind = 8 ) :: log_huge = 80.0D+00
  real ( kind = 8 ) r8_sech
  real ( kind = 8 ) x

  if ( log_huge < abs ( x ) ) then
    r8_sech = 0.0D+00
  else
    r8_sech = 1.0D+00 / cosh ( x )
  end if

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8_sign_char ( x )

!*****************************************************************************80
!
!! R8_SIGN_CHAR returns a character indicating the sign of an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, character R8_SIGN_CHAR, the sign of X, '-', '0' or '+'.
!
  implicit none

  character r8_sign_char
  character value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = '-'
  else if ( x == 0.0D+00 ) then
    value = '0'
  else
    value = '+'
  end if

  r8_sign_char = value

  return
end
function r8_sign_match ( r1, r2 )

!*****************************************************************************80
!
!! R8_SIGN_MATCH is TRUE if two R8's are of the same sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( 0 <= r1 * r2 ) then ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the values to check.
!
!    Output, logical R8_SIGN_MATCH, is TRUE if ( R1 <= 0 and R2 <= 0 )
!    or ( 0 <= R1 and 0 <= R2 ).
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical r8_sign_match

  r8_sign_match = ( r1 <= 0.0D+00 .and. r2 <= 0.0D+00 ) .or. &
                  ( 0.0D+00 <= r1 .and. 0.0D+00 <= r2 )

  return
end
function r8_sign_match_strict ( r1, r2 )

!*****************************************************************************80
!
!! R8_SIGN_MATCH_STRICT is TRUE if two R8's are of the same strict sign.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the values to check.
!
!    Output, logical R8_SIGN_MATCH_STRICT, is TRUE if the signs match.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical r8_sign_match_strict

  r8_sign_match_strict = &
    (           r1 <  0.0D+00 .and. r2 <  0.0D+00 ) .or. &
    (           r1 == 0.0D+00 .and. r2 == 0.0D+00 ) .or. &
    ( 0.0D+00 < r1            .and.       0.0D+00 < r2 )

  return
end
function r8_sign_opposite ( r1, r2 )

!*****************************************************************************80
!
!! R8_SIGN_OPPOSITE is TRUE if two R8's are not of the same sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( r1 * r2 <= 0.0 ) then ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the values to check.
!
!    Output, logical R8_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
!    or ( R2 <= 0 and 0 <= R1 ).
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical r8_sign_opposite

  r8_sign_opposite = ( r1 <= 0.0D+00 .and. 0.0D+00 <= r2 ) .or. &
                     ( r2 <= 0.0D+00 .and. 0.0D+00 <= r1 )

  return
end
function r8_sign_opposite_strict ( r1, r2 )

!*****************************************************************************80
!
!! R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( r1 * r2 < 0.0 ) then ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the values to check.
!
!    Output, logical R8_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
!    or ( R2 < 0 and 0 < R1 ).
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical r8_sign_opposite_strict

  r8_sign_opposite_strict = ( r1 < 0.0D+00 .and. 0.0D+00 < r2 ) .or. &
                            ( r2 < 0.0D+00 .and. 0.0D+00 < r1 )

  return
end
function r8_sqrt_i4 ( i )

!*****************************************************************************80
!
!! R8_SQRT_I4 returns the square root of an I4 as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose square root is desired.
!
!    Output, real ( kind = 8 ) R8_SQRT_I4, the value of sqrt(I).
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_sqrt_i4

  r8_sqrt_i4 = sqrt ( real ( i, kind = 8 ) )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8_swap3 ( x, y, z )

!*****************************************************************************80
!
!! R8_SWAP3 swaps three R8's.
!
!  Example:
!
!    Input:
!
!      X = 1, Y = 2, Z = 3
!
!    Output:
!
!      X = 2, Y = 3, Z = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  w = x
  x = y
  y = z
  z = w

  return
end
function r8_tiny ( )

!*****************************************************************************80
!
!! R8_TINY returns a very small but positive R8.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine TINY ( X ) that
!    is more suitable for this purpose, returning the smallest positive
!    but normalized real number.
!
!    This routine does NOT try to provide an accurate value for TINY.
!    Instead, it simply returns a "reasonable" value, that is, a rather
!    small, but representable, real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_TINY, a "tiny" value.
!
  implicit none

  real ( kind = 8 ) r8_tiny

  r8_tiny = 1.0D-30

  return
end
subroutine r8_to_r8_discrete ( r, rmin, rmax, nr, rd )

!*****************************************************************************80
!
!! R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
!
!  Formula:
!
!    if ( R < RMIN ) then
!      RD = RMIN
!    else if ( RMAX < R ) then
!      RD = RMAX
!    else
!      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
!      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
!
!    In the special case where NR = 1, when
!
!      XD = 0.5 * ( RMAX + RMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, real ( kind = 8 ) RMAX, RMIN, the maximum and minimum
!    values for RD.
!
!    Input, integer ( kind = 4 ) NR, the number of allowed values for XD.
!    NR should be at least 1.
!
!    Output, real ( kind = 8 ) RD, the corresponding discrete value.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) nr
  real ( kind = 8 ) r
  real ( kind = 8 ) rd
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin
!
!  Check for errors.
!
  if ( nr < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TO_R8_DISCRETE - Fatal error!'
    write ( *, '(a,i8)' ) '  NR = ', nr
    write ( *, '(a)' ) '  but NR must be at least 1.'
    stop
  end if

  if ( nr == 1 ) then
    rd = 0.5D+00 * ( rmin + rmax )
    return
  end if

  if ( rmax == rmin ) then
    rd = rmax
    return
  end if

  f = nint ( real ( nr, kind = 8 ) * ( rmax - r ) / ( rmax - rmin ) )
  f = max ( f, 0 )
  f = min ( f, nr )

  rd = ( real (      f, kind = 8 ) * rmin   &
       + real ( nr - f, kind = 8 ) * rmax ) &
       / real ( nr,     kind = 8 )

  return
end
subroutine r8_to_dhms ( r, d, h, m, s )

!*****************************************************************************80
!
!! R8_TO_DHMS converts decimal days into days, hours, minutes, seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, a decimal number representing a time
!    period measured in days.
!
!    Output, integer ( kind = 4 ) D, H, M, S, the equivalent number of days,
!    hours, minutes and seconds.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) r_copy
  integer ( kind = 4 ) s

  r_copy = abs ( r )

  d = int ( r_copy )

  r_copy = r_copy - d
  r_copy = 24.0D+00 * r_copy
  h = int ( r_copy )

  r_copy = r_copy - h
  r_copy = 60.0D+00 * r_copy
  m = int ( r_copy )

  r_copy = r_copy - m
  r_copy = 60.0D+00 * r_copy
  s = int ( r_copy )

  if ( r < 0.0D+00 ) then
    d = -d
    h = -h
    m = -m
    s = -s
  end if

  return
end
subroutine r8_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

!*****************************************************************************80
!
!! R8_TO_I4 maps X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
!
!  Formula:
!
!    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
!    IX := min ( IX, max ( IXMIN, IXMAX ) )
!    IX := max ( IX, min ( IXMIN, IXMAX ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be converted.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the range.  XMAX and
!    XMIN must not be equal.  It is not necessary that XMIN be less than XMAX.
!
!    Input, integer ( kind = 4 ) IXMIN, IXMAX, the allowed range of the output
!    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
!    It is not necessary that IXMIN be less than IXMAX.
!
!    Output, integer ( kind = 4 ) IX, the value in the range [IXMIN,IXMAX] that
!    corresponds to X.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) ixmin
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  if ( xmax == xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  XMAX = XMIN, making a zero divisor.'
    write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
    write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
    stop
  end if

  temp = &
      ( ( xmax - x        ) * real ( ixmin, kind = 8 )  &
      + (        x - xmin ) * real ( ixmax, kind = 8 ) ) &
      / ( xmax     - xmin )

  if ( 0.0D+00 <= temp ) then
    temp = temp + 0.5D+00
  else
    temp = temp - 0.5D+00
  end if

  ix = int ( temp )

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM_AB returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_AB, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_ab = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8_unswap3 ( x, y, z )

!*****************************************************************************80
!
!! R8_UNSWAP3 unswaps three R8's.
!
!  Example:
!
!    Input:
!
!      X = 2, Y = 3, Z = 1
!
!    Output:
!
!      X = 1, Y = 2, Z = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  w = z
  z = y
  y = x
  x = w

  return
end
function r8_walsh_1d ( x, digit )

!*****************************************************************************80
!
!! R8_WALSH_1D evaluates the Walsh function.
!
!  Discussion:
!
!    Consider the binary representation of X, and number the digits
!    in descending order, from leading to lowest, with the units digit
!    being numbered 0.
!
!    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Walsh function.
!
!    Input, integer ( kind = 4 ) DIGIT, the index of the Walsh function.
!
!    Output, real ( kind = 8 ) R8_WALSH_1D, the value of the Walsh function.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_walsh_1d
  real ( kind = 8 ) x
  real ( kind = 8 ) x_copy
!
!  Hide the effect of the sign of X.
!
  x_copy = abs ( x )
!
!  If DIGIT is positive, divide by 2 DIGIT times.
!  If DIGIT is negative, multiply by 2 (-DIGIT) times.
!
  x_copy = x_copy / 2.0D+00**digit
!
!  Make it an integer.
!  Because it's positive, and we're using INT, we don't change the
!  units digit.
!
  n = int ( x_copy )
!
!  Is the units digit odd or even?
!
  if ( mod ( n, 2 ) == 0 ) then
    r8_walsh_1d = 0.0D+00
  else
    r8_walsh_1d = 1.0D+00
  end if

  return
end
function r8_wrap ( r, rlo, rhi )

!*****************************************************************************80
!
!! R8_WRAP forces an R8 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Example:
!
!    RLO = 4.0, RHI = 8.0
!
!     R  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, a value.
!
!    Input, real ( kind = 8 ) RLO, RHI, the desired bounds.
!
!    Output, real ( kind = 8 ) R8_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_wrap
  real ( kind = 8 ) rhi
  real ( kind = 8 ) rhi2
  real ( kind = 8 ) rlo
  real ( kind = 8 ) rlo2
  real ( kind = 8 ) rwide
  real ( kind = 8 ) value
!
!  Guarantee RLO2 < RHI2.
!
  rlo2 = min ( rlo, rhi )
  rhi2 = max ( rlo, rhi )
!
!  Find the width.
!
  rwide = rhi2 - rlo2
!
!  Add enough copies of (RHI2-RLO2) to R so that the
!  result ends up in the interval RLO2 - RHI2.
!
  if ( rwide == 0.0D+00 ) then
    value = rlo
  else if ( r < rlo2 ) then
    n = int ( ( rlo2 - r ) / rwide ) + 1
    value = r + n * rwide
    if ( value == rhi ) then
      value = rlo
    end if
  else
    n = int ( ( r - rlo2 ) / rwide )
    value = r - n * rwide
    if ( value == rlo ) then
      value = rhi
    end if
  end if

  r8_wrap = value

  return
end
subroutine r82_cheby ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R82_CHEBY sets up the Chebyshev abscissas in an R8 interval.
!
!  Discussion:
!
!    The routine sets up a vector of X values spaced between the values
!    XLO and XHI in a similar way to the spacing of the Chebyshev
!    points of the same order in the interval [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to compute.
!
!    Input, real ( kind = 8 ) ALO, AHI, the range.
!
!    Output, real ( kind = 8 ) A(N), the computed X values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else if ( 1 < n ) then

    do i = 1, n

      arg = real ( 2 * i - 1, kind = 8 ) * pi &
          / real ( 2 * n, kind = 8 )

      a(i) = 0.5D+00 * ( ( 1.0D+00 + cos ( arg ) ) * alo &
                       + ( 1.0D+00 - cos ( arg ) ) * ahi )

    end do

  end if

  return
end
function r82_dist_l2 ( a1, a2 )

!*****************************************************************************80
!
!! R82_DIST_L2 returns the L2 distance between a pair of R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The vector L2 norm is defined as:
!
!      sqrt ( sum ( 1 <= I <= N ) A(I) * A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), the vectors.
!
!    Output, real ( kind = 8 ) R82_DIST_L2, the L2 norm of the distance
!    between A1 and A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  real ( kind = 8 ) r82_dist_l2

  r82_dist_l2 = sqrt ( sum ( ( a1(1:dim_num) - a2(1:dim_num) )**2 ) )

  return
end
function r82_eq ( a1, a2 )

!*****************************************************************************80
!
!! R82_EQ == ( A1 == A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 == A2  <=>  A1(1) == A2(1) and A1(2) == A2(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), two R82 vectors to be compared.
!
!    Output, logical R82_EQ, is TRUE if and only if A1 == A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  logical r82_eq

  if ( all ( a1(1:dim_num) == a2(1:dim_num) ) ) then
    r82_eq = .true.
  else
    r82_eq = .false.
  end if

  return
end
function r82_ge ( a1, a2 )

!*****************************************************************************80
!
!! R82_GE == ( A1 >= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 >= A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) >= A2(2) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R92_GE, is TRUE if and only if A1 >= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_ge

  r82_ge = .true.

  do i = 1, dim_num

    if ( a2(i) < a1(i) ) then
      r82_ge = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r82_ge = .false.
      exit
    end if

  end do

  return
end
function r82_gt ( a1, a2 )

!*****************************************************************************80
!
!! R82_GT == ( A1 > A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R2, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) > A2(2) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_gt

  r82_gt = .false.

  do i = 1, dim_num

    if ( a2(i) < a1(i) ) then
      r82_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r82_gt = .false.
      exit
    end if

  end do

  return
end
function r82_le ( a1, a2 )

!*****************************************************************************80
!
!! R82_LE == ( A1 <= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 <= A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) <= A2(2) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_LE, is TRUE if and only if A1 <= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_le

  r82_le = .true.

  do i = 1, dim_num

    if ( a1(i) < a2(i) ) then
      r82_le = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r82_le = .false.
      exit
    end if

  end do

  return
end
function r82_lt ( a1, a2 )

!*****************************************************************************80
!
!! R82_LT == ( A1 < A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) < A2(2) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_lt

  r82_lt = .false.

  do i = 1, dim_num

    if ( a1(i) < a2(i) ) then
      r82_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r82_lt = .false.
      exit
    end if

  end do

  return
end
function r82_ne ( a1, a2 )

!*****************************************************************************80
!
!! R82_NE == ( A1 /= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 /= A2  <=>  A1(1) /= A2(1) or A1(2) /= A2(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_NE, is TRUE if and only if A1 /= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a1(dim_num)
  real ( kind = 8 ) a2(dim_num)
  logical r82_ne

  if ( any ( a1(1:dim_num) /= a2(1:dim_num) ) ) then
    r82_ne = .true.
  else
    r82_ne = .false.
  end if

  return
end
function r82_norm ( a )

!*****************************************************************************80
!
!! R82_NORM returns the Euclidean norm of an R82.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), the vector.
!
!    Output, real ( kind = 8 ) R82_NORM, the norm.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) r82_norm

  r82_norm = sqrt ( a(1) * a(1) + a(2) * a(2) )

  return
end
subroutine r82_normalize ( a )

!*****************************************************************************80
!
!! R82_NORMALIZE Euclidean normalizes an R82.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(2), the components of the vector.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) norm

  norm = sqrt ( a(1) * a(1) + a(2) * a(2) )

  if ( norm /= 0.0D+00 ) then
    a(1:2) = a(1:2) / norm
  end if

  return
end
subroutine r82_print ( a, title )

!*****************************************************************************80
!
!! R82_PRINT prints an R82.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  real ( kind = 8 ) a(2)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', a(1), ',', a(2), ')'
  else
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1 )' ) '(', a(1), ',', a(2), ')'

  end if

  return
end
subroutine r82_swap ( x, y )

!*****************************************************************************80
!
!! R82_SWAP swaps two R82 values.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(2), Y(2).  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) y(dim_num)
  real ( kind = 8 ) z(dim_num)

  z(1:dim_num) = x(1:dim_num)
  x(1:dim_num) = y(1:dim_num)
  y(1:dim_num) = z(1:dim_num)

  return
end
subroutine r82_uniform_ab ( b, c, seed, a )

!*****************************************************************************80
!
!! R82_UNIFORM_AB returns a random R82 value in a given range.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) B, C, the minimum and maximum values.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(2), the randomly chosen value.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed

  do i = 1, dim_num
    a(i) = r8_uniform_ab ( b, c, seed )
  end do

  return
end
subroutine r82poly2_print ( a, b, c, d, e, f )

!*****************************************************************************80
!
!! R82POLY2_PRINT prints a second order polynomial in two variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, E, F, the coefficients.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f

  write ( *, &
    '( 2x, f8.4, '' * x^2 + '', f8.4, '' * y^2 + '', f8.4, '' * xy  + '' )' ) &
    a, b, c

  write ( *, &
    '( 2x, f8.4, '' * x + '', f8.4, '' * y + '', f8.4, '' = 0 '' )' ) d, e, f

  return
end
subroutine r82poly2_type ( a, b, c, d, e, f, type )

!*****************************************************************************80
!
!! R82POLY2_TYPE analyzes a second order polynomial in two variables.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A x^2 + B y^2 + C xy + Dx + Ey + F = 0
!
!    The possible types of the solution set are:
!
!     1: a hyperbola;
!        9x^2 -  4y^2       -36x - 24y -  36 = 0
!     2: a parabola;
!        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
!     3: an ellipse;
!        9x^2 + 16y^2       +36x - 32y -  92 = 0;
!     4: an imaginary ellipse (no real solutions);
!         x^2 +   y^2       - 6x - 10y + 115 = 0;
!     5: a pair of intersecting lines;
!                        xy + 3x -   y -   3 = 0
!     6: one point;
!         x^2 +  2y^2       - 2x + 16y +  33 = 0;
!     7: a pair of distinct parallel lines;
!                 y^2            -  6y +   8 = 0
!     8: a pair of imaginary parallel lines (no real solutions);
!                 y^2            -  6y +  10 = 0
!     9: a pair of coincident lines.
!                 y^2            -  2y +   1 = 0
!    10: a single line;
!                             2x -   y +   1 = 0;
!    11; all space;
!                                          0 = 0;
!    12; no solutions;
!                                          1 = 0;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 30th Edition, 1996, pages 282-284.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, E, F, the coefficients.
!
!    Output, integer ( kind = 4 ) TYPE, indicates the type of the solution set.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) delta
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) j
  real ( kind = 8 ) k
  integer ( kind = 4 ) type
!
!  Handle the degenerate case.
!
  if ( a == 0.0D+00 .and. &
       b == 0.0D+00 .and. &
       c == 0.0D+00 ) then
    if ( d == 0.0D+00 .and. e == 0.0D+00 ) then
      if ( f == 0.0D+00 ) then
        type = 11
      else
        type = 12
      end if
    else
      type = 10
    end if
    return
  end if

  delta = &
      8.0D+00 * a * b * f &
    + 2.0D+00 * c * e * d &
    - 2.0D+00 * a * e * e &
    - 2.0D+00 * b * d * d &
    - 2.0D+00 * f * c * c

  j = 4.0D+00 * a * b - c * c

  if ( delta /= 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 1
    else if ( j == 0.0D+00 ) then
      type = 2
    else if ( 0.0D+00 < j ) then
      if ( sign ( 1.0D+00, delta ) /= sign ( 1.0D+00, ( a + b ) ) ) then
        type = 3
      else if ( sign ( 1.0D+00, delta ) == sign ( 1.0D+00, ( a + b ) ) ) then
        type = 4
      end if
    end if
  else if ( delta == 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 5
    else if ( 0.0D+00 < j ) then
      type = 6
    else if ( j == 0.0D+00 ) then

      k = 4.0D+00 * ( a + b ) * f - d * d - e * e

      if ( k < 0.0D+00 ) then
        type = 7
      else if ( 0.0D+00 < k ) then
        type = 8
      else if ( k == 0.0D+00 ) then
        type = 9
      end if

    end if
  end if

  return
end
subroutine r82poly2_type_print ( type )

!*****************************************************************************80
!
!! R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TYPE, the type index returned by R82POLY2_TYPE.
!
  implicit none

  integer ( kind = 4 ) type

  if ( type == 1 ) then
    write ( *, '(a)' ) '  The set of solutions forms a hyperbola.'
  else if ( type == 2 ) then
    write ( *, '(a)' ) '  The set of solutions forms a parabola.'
  else if ( type == 3 ) then
    write ( *, '(a)' ) '  The set of solutions forms an ellipse.'
  else if ( type == 4 ) then
    write ( *, '(a)' ) '  The set of solutions forms an imaginary ellipse.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 5 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of intersecting lines.'
  else if ( type == 6 ) then
    write ( *, '(a)' ) '  The set of solutions is a single point.'
  else if ( type == 7 ) then
    write ( *, '(a)' ) &
      '  The set of solutions form a pair of distinct parallel lines.'
  else if ( type == 8 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of imaginary parallel lines.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 9 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of coincident lines.'
  else if ( type == 10 ) then
    write ( *, '(a)' ) '  The set of solutions forms a single line.'
  else if ( type == 11 ) then
    write ( *, '(a)' ) '  The set of solutions is all space.'
  else if ( type == 12 ) then
    write ( *, '(a)' ) '  The set of solutions is empty.'
  else
    write ( *, '(a)' ) '  This type index is unknown.'
  end if

  return
end
subroutine r82vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R82VEC_MAX returns the maximum value in an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array.
!
!    Output, real ( kind = 8 ) AMAX(2); the largest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) amax(2)

  amax(1) = maxval ( a(1,1:n) )
  amax(2) = maxval ( a(2,1:n) )

  return
end
subroutine r82vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R82VEC_MIN returns the minimum value in an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array.
!
!    Output, real ( kind = 8 ) AMIN(2); the smallest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) amin(2)

  amin(1) = minval ( a(1,1:n) )
  amin(2) = minval ( a(2,1:n) )

  return
end
subroutine r82vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R82VEC_ORDER_TYPE finds the order type of an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8 values.
!
!    The dictionary or lexicographic ordering is used.
!
!    (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1,1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( &
         a(1,1) <  a(1,i) .or. &
       ( a(1,1) == a(1,i) .and. a(2,1) < a(2,i) ) &
    ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( &
        a(1,i) <  a(1,1)  .or. &
      ( a(1,i) == a(1,1) .and. a(2,i) < a(2,1) ) &
    ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do

    i = i + 1
    if ( n < i ) then
      exit
    end if

    if ( order == 1 ) then

      if ( &
          a(1,i) <  a(1,i-1) .or. &
        ( a(1,i) == a(1,i-1) .and. a(2,i) < a(2,i-1) ) &
      ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( &
          a(1,i) <  a(1,i-1) .or. &
        ( a(1,i) == a(1,i-1) .and. a(2,i) < a(2,i-1) ) &
      ) then
        order = -1
        exit
      else if ( &
         a(1,i) == a(1,i-1) .and. a(2,i) == a(2,i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( &
          a(1,i-1) <  a(1,i) .or. &
        ( a(1,i-1) == a(1,i) .and. a(2,i-1) < a(2,i) ) &
      ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( &
          a(1,i-1) <  a(1,i) .or. &
        ( a(1,i-1) == a(1,i) .and. a(2,i-1) < a(2,i) ) &
      ) then
        order = -1
        exit
      else if ( a(1,i) == a(1,i-1) .and. a(2,i) == a(2,i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine r82vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82 values.
!
!    The routine reorders the entries of A.  Using A(1:2,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
!             -----------          ----------------------------------
!             LEFT          KEY    RIGHT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three
!    segments.  Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    KEY < A(1:2,I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) key(dim_num)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:dim_num) = a(1:dim_num,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      r = r - 1
      call r8vec_swap ( dim_num, a(1:dim_num,r), a(1:dim_num,l+1) )
    else if ( r8vec_eq ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      m = m + 1
      call r8vec_swap ( dim_num, a(1:dim_num,m), a(1:dim_num,l+1) )
      l = l + 1
    else if ( r8vec_lt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:dim_num,i) = a(1:dim_num,i+m)
  end do

  l = l - m

  do i = 1, dim_num
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r82vec_permute ( n, p, a )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8 values.
!
!    The same logic can be used to permute an array of objects of any
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  real ( kind = 8 ) a_temp(dim_num)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)

    else

      a_temp(1:dim_num) = a(1:dim_num,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:dim_num,iput) = a_temp(1:dim_num)
          exit
        end if

        a(1:dim_num,iput) = a(1:dim_num,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine r82vec_print ( n, a, title )

!*****************************************************************************80
!
!! R82VEC_PRINT prints an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,(5g14.6))' ) i, a(1:dim_num,i)
  end do

  return
end
subroutine r82vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R82VEC_PRINT_PART prints "part" of an R82VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(2,N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(1:2,i), &
      '...more entries...'

  end if

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r82vec_permute ( n, indx, a )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  real ( kind = 8 ) aval(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:dim_num) = a(1:dim_num,indxt)

    else

      indxt = indx(ir)
      aval(1:dim_num) = a(1:dim_num,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r82vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 30
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r82vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
function r83_norm ( x, y, z )

!*****************************************************************************80
!
!! R83_NORM returns the Euclidean norm of an R83.
!
!  Discussion:
!
!    An R83 is a vector of 3 R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, Z, the vector.
!
!    Output, real ( kind = 8 ) R83_NORM, the norm of the vector.
!
  implicit none

  real ( kind = 8 ) r83_norm
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  r83_norm = sqrt ( x * x + y * y + z * z )

  return
end
subroutine r83_normalize ( x, y, z )

!*****************************************************************************80
!
!! R83_NORMALIZE normalizes an R83.
!
!  Discussion:
!
!    An R83 is a vector of 3 R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, the components of the vector.
!
  implicit none

  real ( kind = 8 ) norm
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  norm = sqrt ( x * x + y * y + z * z )

  if ( norm /= 0.0D+00 ) then
    x = x / norm
    y = y / norm
    z = z / norm
  end if

  return
end
subroutine r83_print ( x, y, z, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83.
!
!  Discussion:
!
!    An R83 is a vector of 3 R8's.
!
!    A format is used which suggests a coordinate triple:
!
!  Example:
!
!    Center : ( 1.23, 7.45, -1.45 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, Z, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  character ( len = * ) title
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  if ( 0 < len_trim ( title ) ) then
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', x, ',', y, ',', z, ')'
  else
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      '(', x, ',', y, ',', z, ')'
  end if

  return
end
subroutine r83_swap ( x, y )

!*****************************************************************************80
!
!! R83_SWAP swaps two R83's.
!
!  Discussion:
!
!    An R83 is a vector of 3 R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(3), Y(3).  On output, the values
!    of X and Y have been interchanged.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) y(dim_num)
  real ( kind = 8 ) z(dim_num)

  z(1:dim_num) = x(1:dim_num)
  x(1:dim_num) = y(1:dim_num)
  y(1:dim_num) = z(1:dim_num)

  return
end
subroutine r83vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R83VEC_MAX returns the maximum value in an R83VEC.
!
!  Discussion:
!
!    An R83VEC is an array of R83's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(3,N), the array.
!
!    Output, real ( kind = 8 ) AMAX(3); the largest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) amax(3)
  integer ( kind = 4 ) i

  do i = 1, 3
    amax(i) = maxval ( a(i,1:n) )
  end do

  return
end
subroutine r83vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R83VEC_MIN returns the minimum value in an R83VEC.
!
!  Discussion:
!
!    An R83VEC is an array of R83's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(3,N), the array.
!
!    Output, real ( kind = 8 ) AMIN(3); the smallest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) amin(3)
  integer ( kind = 4 ) i

  do i = 1, 3
    amin(i) = minval ( a(i,1:n) )
  end do

  return
end
subroutine r83vec_normalize ( n, x )

!*****************************************************************************80
!
!! R83VEC_NORMALIZE normalizes each R83 in an R83VEC.
!
!  Discussion:
!
!    An R83VEC is a vector of R83's.
!
!    An R83 is a vector of 3 R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of R83 vectors.
!
!    Input/output, real ( kind = 8 ) X(3,N), the coordinates of N R83 vectors.
!    On output, the nonzero vectors have been scaled to have unit L2 norm.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) i
  real ( kind = 8 ) norm
  real ( kind = 8 ) x(dim_num,n)

  do i = 1, n

    norm = sqrt ( sum ( x(1:dim_num,i)**2 ) )

    if ( norm /= 0.0D+00 ) then
      x(1:dim_num,i) = x(1:dim_num,i) / norm
    end if

  end do

  return
end
subroutine r83vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R83VEC_PRINT_PART prints "part" of an R83VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(3,N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do
    write ( *, '(a)' ) &
      '  ........  ..............  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(1:3,i), &
      '...more entries...'

  end if

  return
end
subroutine r84_normalize ( v )

!*****************************************************************************80
!
!! R84_NORMALIZE normalizes an R84.
!
!  Discussion:
!
!    An R84 is a vector of four R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) V(4), the components of the vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4

  real ( kind = 8 ) norm
  real ( kind = 8 ) v(dim_num)

  norm = sqrt ( sum ( v(1:dim_num)**2 ) )

  if ( norm /= 0.0D+00 ) then
    v(1:dim_num) = v(1:dim_num) / norm
  end if

  return
end
subroutine r8block_expand_linear ( l, m, n, x, lfat, mfat, nfat, xfat )

!*****************************************************************************80
!
!! R8BLOCK_EXPAND_LINEAR linearly interpolates new data into an R8BLOCK.
!
!  Discussion:
!
!    An R8BLOCK is a 3D array of R8 values.
!
!    In this routine, the expansion is specified by giving the number
!    of intermediate values to generate between each pair of original
!    data rows and columns.
!
!    The interpolation is not actually linear.  It uses the functions
!
!      1, x, y, z, xy, xz, yz, xyz.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the input data.
!
!    Input, real ( kind = 8 ) X(L,M,N), the original data.
!
!    Input, integer ( kind = 4 ) LFAT, MFAT, NFAT, the number of data values
!    to interpolate original data values in the first, second and third
!    dimensions.
!
!    Output, real ( kind = 8 ) XFAT(L2,M2,N2), the fattened data, where
!    L2 = (L-1)*(LFAT+1)+1,
!    M2 = (M-1)*(MFAT+1)+1,
!    N2 = (N-1)*(NFAT+1)+1.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lfat
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mfat
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkk
  integer ( kind = 4 ) kp1
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(l,m,n)
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xfat((l-1)*(lfat+1)+1,(m-1)*(mfat+1)+1,(n-1)*(nfat+1)+1)

  do i = 1, l

    if ( i < l ) then
      ihi = lfat
    else
      ihi = 0
    end if

    do j = 1, m

      if ( j < m ) then
        jhi = mfat
      else
        jhi = 0
      end if

      do k = 1, n

        if ( k < n ) then
          khi = nfat
        else
          khi = 0
        end if

        if ( i < l ) then
          ip1 = i + 1
        else
          ip1 = i
        end if

        if ( j < m ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        if ( k < n ) then
          kp1 = k + 1
        else
          kp1 = k
        end if

        x000 = x(i,j,k)
        x001 = x(i,j,kp1)
        x100 = x(ip1,j,k)
        x101 = x(ip1,j,kp1)
        x010 = x(i,jp1,k)
        x011 = x(i,jp1,kp1)
        x110 = x(ip1,jp1,k)
        x111 = x(ip1,jp1,kp1)

        do ii = 0, ihi

          r = real ( ii,      kind = 8 ) &
            / real ( ihi + 1, kind = 8 )

          do jj = 0, jhi

            s = real ( jj,      kind = 8 ) &
              / real ( jhi + 1, kind = 8 )

            do kk = 0, khi

              t = real ( kk,      kind = 8 ) &
                / real ( khi + 1, kind = 8 )

              iii = 1 + ( i - 1 ) * ( lfat + 1 ) + ii
              jjj = 1 + ( j - 1 ) * ( mfat + 1 ) + jj
              kkk = 1 + ( k - 1 ) * ( nfat + 1 ) + kk

              xfat(iii,jjj,kkk) = &
                  x000 * ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) &
                + x001 * ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * (           t ) &
                + x010 * ( 1.0D+00 - r ) * (           s ) * ( 1.0D+00 - t ) &
                + x011 * ( 1.0D+00 - r ) * (           s ) * (           t ) &
                + x100 * (           r ) * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) &
                + x101 * (           r ) * ( 1.0D+00 - s ) * (           t ) &
                + x110 * (           r ) * (           s ) * ( 1.0D+00 - t ) &
                + x111 * (           r ) * (           s ) * (           t )

            end do

          end do

        end do

      end do

    end do

  end do

  return
end
subroutine r8block_print ( l, m, n, a, title )

!*****************************************************************************80
!
!! R8BLOCK_PRINT prints an R8BLOCK.
!
!  Discussion:
!
!    An R8BLOCK is a 3D array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the block.
!
!    Input, real ( kind = 8 ) A(L,M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(l,m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do k = 1, n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K = ', k

    do jlo = 1, m, 5
      jhi = min ( jlo + 4, m )
      write ( *, '(a)' ) ' '
      write ( *, '(10x,5(i8,6x))' ) (j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = 1, l
        write ( *, '(2x,i8,5g14.6)' ) i, a(i,jlo:jhi,k)
      end do
    end do

  end do

  return
end
subroutine r8col_compare ( m, n, a, i, j, value )

!*****************************************************************************80
!
!! R8COL_COMPARE compares columns in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      VALUE = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) VALUE, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) value
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  value = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      value = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      value = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8col_duplicates ( m, n, n_unique, seed, a )

!*****************************************************************************80
!
!! R8COL_DUPLICATES generates an R8COL with some duplicate columns.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    This routine generates a random R8COL with a specified number of
!    duplicate columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) N_UNIQUE, the number of unique columns in A.
!    1 <= N_UNIQUE <= N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) A(M,N), the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp(m)

  if ( n_unique < 1 .or. n < n_unique ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_DUPLICATES - Fatal error!'
    write ( *, '(a)' ) '  1 <= N_UNIQUE <= N is required.'
    stop
  end if

  call r8mat_uniform_01 ( m, n_unique, seed, a )
!
!  Randomly copy unique columns.
!
  do j1 = n_unique + 1, n
    j2 = i4_uniform_ab ( 1, n_unique, seed )
    a(1:m,j1) = a(1:m,j2)
  end do
!
!  Permute the columns.
!
  do j1 = 1, n
    j2 = i4_uniform_ab ( j1, n, seed )
    temp(1:m) = a(1:m,j1)
    a(1:m,j1) = a(1:m,j2)
    a(1:m,j2) = temp(1:m)
  end do

  return
end
subroutine r8col_find ( m, n, a, x, col )

!*****************************************************************************80
!
!! R8COL_FIND seeks a column value in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      COL = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real ( kind = 8 ) X(M), a vector to be matched with a column of A.
!
!    Output, integer ( kind = 4 ) COL, the index of the first column of A
!    which exactly matches every entry of X, or -1 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)

  col = -1

  do j = 1, n

    col = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        col = -1
        exit
      end if
    end do

    if ( col /= -1 ) then
      return
    end if

  end do

  return
end
subroutine r8col_first_index ( m, n, a, tol, first_index )

!*****************************************************************************80
!
!! R8COL_FIRST_INDEX indexes the first occurrence of values in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
!    the first column whose entries are equal to A(1:M,J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!    The length of an "element" of A, and the number of "elements".
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) FIRST_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) first_index(n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) tol

  first_index(1:n) = -1

  do j1 = 1, n

    if ( first_index(j1) == -1 ) then

      first_index(j1) = j1

      do j2 = j1 + 1, n
        if ( maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) ) <= tol ) then
          first_index(j2) = j1
        end if
      end do

    end if

  end do

  return
end
subroutine r8col_insert ( n_max, m, n, a, x, col )

!*****************************************************************************80
!
!! R8COL_INSERT inserts a column into an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      N_MAX = 10,
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      X = ( 3., 4., 18. )
!
!    Output:
!
!      N = 5,
!
!      A = (
!        1.  2.  3.  3.  4.
!        5.  6.  4.  7.  8.
!        9. 10. 18. 11. 12. )
!
!      COL = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum number of columns in A.
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input/output, integer ( kind = 4 ) N, the number of columns.
!    If the new column is inserted into the table, then the output
!    value of N will be increased by 1.
!
!    Input/output, real ( kind = 8 ) A(M,N_MAX), a table of numbers, regarded
!    as an array of columns.  The columns must have been sorted
!    lexicographically.
!
!    Input, real ( kind = 8 ) X(M), a vector of data which will be inserted
!    into the table if it does not already occur.
!
!    Output, integer ( kind = 4 ) COL.
!    I, X was inserted into column I.
!    -I, column I was already equal to X.
!    0, N = N_MAX.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_max

  real ( kind = 8 ) a(m,n_max)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) high
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) n
  real ( kind = 8 ) x(m)
!
!  Refuse to work if N_MAX <= N.
!
  if ( n_max <= n ) then
    col = 0
    return
  end if
!
!  Stick X temporarily in column N+1, just so it's easy to use R8COL_COMPARE.
!
  a(1:m,n+1) = x(1:m)
!
!  Do a binary search.
!
  low = 1
  high = n

  do

    if ( high < low ) then
      col = low
      exit
    end if

    mid = ( low + high ) / 2

    call r8col_compare ( m, n + 1, a, mid, n + 1, isgn )

    if ( isgn == 0 ) then
      col = -mid
      return
    else if ( isgn == -1 ) then
      low = mid + 1
    else if ( isgn == +1 ) then
      high = mid - 1
    end if

  end do
!
!  Shift part of the table up to make room.
!
  do j = n, col, -1
    a(1:m,j+1) = a(1:m,j)
  end do
!
!  Insert the new column.
!
  a(1:m,col) = x(1:m)

  n = n + 1

  return
end
subroutine r8col_max ( m, n, a, amax )

!*****************************************************************************80
!
!! R8COL_MAX returns the maximums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMAX(N), the maximums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax(n)
  integer ( kind = 4 ) j

  do j = 1, n

    amax(j) = maxval ( a(1:m,j) )

  end do

  return
end
subroutine r8col_max_index ( m, n, a, imax )

!*****************************************************************************80
!
!! R8COL_MAX_INDEX returns the indices of column maximums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IMAX(N); IMAX(I) is the row of A in which
!    the maximum for column I occurs.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax(n)
  integer ( kind = 4 ) j

  do j = 1, n

    imax(j) = 1
    amax = a(1,j)
    do i = 2, m
      if ( amax < a(i,j) ) then
        imax(j) = i
        amax = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_max_one ( m, n, a )

!*****************************************************************************80
!
!! R8COL_MAX_ONE rescales an R8COL so each column maximum is 1.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the array to be rescaled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_big
  integer ( kind = 4 ) j

  do j = 1, n

    i_big = 1
    do i = 2, m
      if ( abs ( a(i_big,j) ) < abs ( a(i,j) ) ) then
        i_big = i
      end if
    end do

    if ( a(i_big,j) /= 0.0D+00 ) then
      a(1:m,j) = a(1:m,j) / a(i_big,j)
    end if

  end do

  return
end
subroutine r8col_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! R8COL_MEAN returns the column means of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      1.5  4.0  5.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) MEAN(N), the means, or averages, of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(n)

  do j = 1, n
    mean(j) = sum ( a(1:m,j) )
  end do

  mean(1:n) = mean(1:n) / real ( m, kind = 8  )

  return
end
subroutine r8col_min ( m, n, a, amin )

!*****************************************************************************80
!
!! R8COL_MIN returns the column minimums of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMIN(N), the minimums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amin(n)
  integer ( kind = 4 ) j

  do j = 1, n

    amin(j) = minval ( a(1:m,j) )

  end do

  return
end
subroutine r8col_min_index ( m, n, a, imin )

!*****************************************************************************80
!
!! R8COL_MIN_INDEX returns the indices of column minimums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IMIN(N); IMIN(I) is the row of A in which
!    the minimum for column I occurs.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin(n)
  integer ( kind = 4 ) j

  do j = 1, n

    imin(j) = 1
    amin = a(1,j)
    do i = 2, m
      if ( a(i,j) < amin ) then
        imin(j) = i
        amin = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_normalize_li ( m, n, a )

!*****************************************************************************80
!
!! R8COL_NORMALIZE_LI normalizes an R8COL with the column infinity norm.
!
!  Discussion:
!
!    Each column is scaled so that the entry of maximum norm has the value 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the array to be normalized.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    c = a(1,j)

    do i = 2, m
      if ( abs ( c ) < abs ( a(i,j) ) ) then
        c = a(i,j)
      end if
    end do

    if ( c /= 0.0D+00 ) then
      a(1:m,j) = a(1:m,j) / c
    end if

  end do

  return
end
subroutine r8col_part_quick_a ( m, n, a, l, r )

!*****************************************************************************80
!
!! R8COL_PART_QUICK_A reorders the columns of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The routine reorders the columns of A.  Using A(1:M,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      M = 2, N = 8
!      A = ( 2  8  6  0 10 10  0  5
!            4  8  2  2  6  0  6  8 )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = (  0  0  2  8  6 10 10  5
!             2  6  4  8  2  6  0  8 )
!             ----     -------------
!             LEFT KEY     RIGHT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row dimension of A, and the length of
!    a column.
!
!    Input, integer ( kind = 4 ) N, the column dimension of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three
!    segments.  Let KEY = the input value of A(1:M,1).  Then
!    I <= L                 A(1:M,I) < KEY;
!         L < I < R         A(1:M,I) = KEY;
!                 R <= I    KEY < A(1:M,I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) key(m)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    return
  end if

  if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:m) = a(1:m,1)
  k = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do j = 2, n

    if ( r8vec_gt ( m, a(1:m,l+1), key(1:m) ) ) then
      r = r - 1
      call r8vec_swap ( m, a(1:m,r), a(1:m,l+1) )
    else if ( r8vec_eq ( m, a(1:m,l+1), key(1:m) ) ) then
      k = k + 1
      call r8vec_swap ( m, a(1:m,k), a(1:m,l+1) )
      l = l + 1
    else if ( r8vec_lt ( m, a(1:m,l+1), key(1:m) ) ) then
      l = l + 1
    end if

  end do
!
!  Shift small elements to the left.
!
  do j = 1, l - k
    a(1:m,j) = a(1:m,j+k)
  end do
!
!  Shift KEY elements to center.
!
  do j = l - k + 1, l
    a(1:m,j) = key(1:m)
  end do
!
!  Update L.
!
  l = l - k

  return
end
subroutine r8col_permute ( m, n, p, a )

!*****************************************************************************80
!
!! R8COL_PERMUTE permutes an R8COL in place.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The same logic can be used to permute an array of objects of any
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      M = 2
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of objects.
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, real ( kind = 8 ) A(M,N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_temp(m)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:m) = a(1:m,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:m,iput) = a_temp(1:m)
          exit
        end if

        a(1:m,iput) = a(1:m,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine r8col_reverse ( m, n, a )

!*****************************************************************************80
!
!! R8COL_REVERSE reverses the order of columns in an R8COL.
!
!  Discussion:
!
!    To reverse the columns is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      15 14 13 12 11
!      25 24 23 22 21
!      35 34 33 32 31
!      45 44 43 42 41
!      55 54 53 52 51
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  real ( kind = 8 ) t(m)

  jhi = n / 2

  do j = 1, jhi
    t(1:m)       = a(1:m,j)
    a(1:m,j)     = a(1:m,n+1-j)
    a(1:m,n+1-j) = t(1:m)
  end do

  return
end
subroutine r8col_sort_heap_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) t(m)

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  indx = 0
  isgn = 0
  j1 = 0
  j2 = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, j1, j2, isgn )
!
!  Interchange columns J1 and J2.
!
    if ( 0 < indx ) then

      t(1:m)    = a(1:m,j1)
      a(1:m,j1) = a(1:m,j2)
      a(1:m,j2) = t(1:m)
!
!  Compare columns J1 and J2.
!
    else if ( indx < 0 ) then

      isgn = 0

      do i = 1, m 
 
        if ( a(i,j1) < a(i,j2) ) then
          isgn = -1
          exit
        else if ( a(i,j2) < a(i,j1) ) then
          isgn = +1
          exit
        end if

      end do
!
!  The columns are sorted.
!
    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sort_heap_index_a ( m, n, a, indx )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2)
!    is negative.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(*,INDX(*)) is sorted,
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The I-th element
!    of the sorted array is column INDX(I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) column(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = ( n / 2 ) + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      column(1:m) = a(1:m,indxt)

    else

      indxt = indx(ir)
      column(1:m) = a(1:m,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        call r8vec_compare ( m, a(1:m,indx(j)), a(1:m,indx(j+1)), isgn )

        if ( isgn < 0 ) then
          j = j + 1
        end if

      end if

      call r8vec_compare ( m, column, a(1:m,indx(j)), isgn )

      if ( isgn < 0 ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8col_sort_quick_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_QUICK_A ascending quick sorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row order of A, and the length of
!    a column.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 30
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( m <= 0 ) then
    return
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  end if

  if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r8col_part_quick_a ( m, n_segment, a(1:m,base:base+n_segment-1), &
      l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8col_sorted_tol_undex ( m, n, a, unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_SORTED_TOL_UNDEX indexes tolerably unique entries in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the tolerably unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    Assuming A is already sorted, we examine the entries of A in order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector A could be
!    replaced by a compressed vector XU, containing the unique entries
!    of A in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the unique sort and
!    inverse unique sort vectors and the compressed unique sorted vector.
!
!      I      A      XU  Undx  Xdnu
!    ----+------+------+-----+-----+
!      1 | 11.0 |  11.0    1     1
!      2 | 11.0 |  22.0    5     1
!      3 | 11.0 |  33.0    8     1
!      4 | 11.0 |  55.0    9     1
!      5 | 22.0 |                2
!      6 | 22.0 |                2
!      7 | 22.0 |                2
!      8 | 33.0 |                3
!      9 | 55.0 |                4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(unique_num)
  logical unique
  integer ( kind = 4 ) xdnu(n)
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = i
  xdnu(i) = k
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      i2 = undx(j)
      diff = maxval ( abs ( a(1:m,i) - a(1:m,i2) ) )
      if ( diff <= tol ) then
        unique = .false.
        xdnu(i) = j
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = i
      xdnu(i) = k
    end if

  end do

  return
end
subroutine r8col_sorted_tol_unique ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_TOL_UNIQUE keeps tolerably unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of columns of M-vectors.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol
  logical unique
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    unique = .true.

    do j = 1, unique_num
      diff = maxval ( abs ( a(1:m,j) - a(1:m,i) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      unique_num = unique_num + 1
      a(1:m,unique_num) = a(1:m,i)
    end if

  end do

  return
end
subroutine r8col_sorted_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_TOL_UNIQUE_COUNT: tolerably unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  logical unique
  integer unique_num
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!
  i = 1
  k = 1
  undx(k) = i
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K and set the
!  K-th unique item to I.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      i2 = undx(j)
      diff = maxval ( abs ( a(1:m,i) - a(1:m,i2) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = i
    end if

  end do

  unique_num = k

  return
end
subroutine r8col_sorted_undex ( m, n, a, unique_num, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_SORTED_UNDEX returns unique sorted indexes for a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    Assuming A is already sorted, we examine the entries of A in order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector A could be
!    replaced by a compressed vector XU, containing the unique entries
!    of A in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I      A      XU  Undx  Xdnu
!    ----+------+------+-----+-----+
!      1 | 11.0 |  11.0    1     1
!      2 | 11.0 |  22.0    5     1
!      3 | 11.0 |  33.0    8     1
!      4 | 11.0 |  55.0    9     1
!      5 | 22.0 |                2
!      6 | 22.0 |                2
!      7 | 22.0 |                2
!      8 | 33.0 |                3
!      9 | 55.0 |                4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) AL(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) undx(unique_num)
  integer ( kind = 4 ) xdnu(n)
!
!  Walk through the sorted array.
!
  i = 1
  j = 1
  undx(j) = i
  xdnu(i) = j

  do i = 2, n

    if ( any ( a(1:m,i) /= a(1:m,j) ) ) then
      j = j + 1
      undx(j) = i
    end if

    xdnu(i) = j

  end do

  return
end
subroutine r8col_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE keeps unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of columns of M-vectors.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      j1 = j1 + 1
      a(1:m,j1) = a(1:m,j2)
    end if

  end do

  unique_num = j1

  return
end
subroutine r8col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine r8col_sortr_a ( m, n, a, key )

!*****************************************************************************80
!
!! R8COL_SORTR_A ascending sorts one column of an R8COL, adjusting all columns.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, an unsorted M by N array.
!    On output, rows of the array have been shifted in such
!    a way that column KEY of the array is in nondecreasing order.
!
!    Input, integer ( kind = 4 ) KEY, the column in which the "key" value
!    is stored.  On output, column KEY of the array will be
!    in nondecreasing order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) key

  if ( m <= 0 ) then
    return
  end if

  if ( key < 1 .or. n < key ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SORTR_A - Fatal error!'
    write ( *, '(a)' ) '  The value of KEY is not a legal column index.'
    write ( *, '(a,i8)' ) '  KEY = ', key
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( a(i,key) < a(j,key) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sum ( m, n, a, colsum )

!*****************************************************************************80
!
!! R8COL_SUM sums the columns of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) COLSUM(N), the sums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) colsum(n)
  integer ( kind = 4 ) j

  do j = 1, n
    colsum(j) = sum ( a(1:m,j) )
  end do

  return
end
subroutine r8col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! R8COL_SWAP swaps columns I and J of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  NCOL = ', n
    stop
  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m) = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine r8col_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8COL_TO_R8VEC converts an R8COL to an R8VEC.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, real ( kind = 8 ) X(M*N), a vector containing the N columns of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(m*n)

  k = 1
  do j = 1, n
    x(k:k+m-1) = a(1:m,j)
    k = k + m
  end do

  return
end
subroutine r8col_tol_undex ( m, n, a, unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_TOL_UNDEX indexes tolerably unique entries of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of A, which creates arrays INDX and XDNI.  (If all the entries
!    of A are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of A in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    A   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is A(3).
!    XDNI(2) = 5 means that A(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
!    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = A(I).
!    XU(I)        = A(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(unique_num)
  logical unique
  integer ( kind = 4 ) xdnu(n)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = indx(i)
  xdnu(indx(i)) = k
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )
      if ( diff <= tol ) then
        unique = .false.
        xdnu(indx(i)) = j
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = indx(i)
      xdnu(indx(i)) = k
    end if

  end do

  return
end
subroutine r8col_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  logical unique
  integer ( kind = 4 ) unique_num
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = indx(i)
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = indx(i)
    end if

  end do

  unique_num = k

  return
end
subroutine r8col_tol_unique_index ( m, n, a, tol, unique_index )

!*****************************************************************************80
!
!! R8COL_TOL_UNIQUE_INDEX indexes tolerably unique entries in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
!    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
!    gathered in order, then
!
!      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do j1 = 1, n

    if ( unique_index(j1) == -1 ) then

      unique_num = unique_num + 1
      unique_index(j1) = unique_num

      do j2 = j1 + 1, n
        diff = maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) )
        if ( diff <= tol ) then
          unique_index(j2) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine r8col_undex ( m, n, a, unique_num, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_UNDEX returns unique sorted indexes for an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of A, which creates arrays INDX and XDNI.  (If all the entries
!    of A are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of A in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    A   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is A(3).
!    XDNI(2) = 5 means that A(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
!    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = A(I).
!    XU(I)        = A(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) undx(unique_num)
  integer ( kind = 4 ) xdnu(n)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Walk through the implicitly sorted array.
!
  i = 1
  j = 1
  undx(j) = indx(i)
  xdnu(indx(i)) = j

  do i = 2, n

    diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )

    if ( 0.0D+00 < diff ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine r8col_uniform_abvec ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8COL_UNIFORM_ABVEC fills an R8COL with scaled pseudorandom numbers.
!
!  Discussion:
!
!    An R8COL is an array of R8 values, regarded as a set of column vectors.
!
!    The user specifies a minimum and maximum value for each row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a(i) &
        + ( b(i) - a(i) ) * real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8col_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8COL_UNIQUE_COUNT counts the unique columns in an unsorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  logical unique(n)
  integer ( kind = 4 ) unique_num

  unique_num = 0

  do j1 = 1, n

    unique_num = unique_num + 1
    unique(j1) = .true.

    do j2 = 1, j1 - 1

      if ( unique(j2) ) then
        diff = maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) )
        if ( diff == 0.0D+00 ) then
          unique_num = unique_num - 1
          unique(j1) = .false.
          exit
        end if
      end if

    end do

  end do

  return
end
subroutine r8col_unique_index ( m, n, a, unique_index )

!*****************************************************************************80
!
!! R8COL_UNIQUE_INDEX indexes the unique occurrence of values in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
!    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
!    gathered in order, then
!
!      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!    The length of an "element" of A, and the number of "elements".
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do j1 = 1, n

    if ( unique_index(j1) == -1 ) then

      unique_num = unique_num + 1
      unique_index(j1) = unique_num

      do j2 = j1 + 1, n
        diff = maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) )
        if ( diff == 0.0D+00 ) then
          unique_index(j2) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine r8col_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! R8COL_VARIANCE returns the variances of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array whose variances are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(N), the variances of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(n)

  do j = 1, n

    mean = sum ( a(1:m,j) ) / real ( m, kind = 8  )

    variance(j) = 0.0D+00
    do i = 1, m
      variance(j) = variance(j) + ( a(i,j) - mean )**2
    end do

    if ( 1 < m ) then
      variance(j) = variance(j) / real ( m - 1, kind = 8 )
    else
      variance(j) = 0.0D+00
    end if

  end do

  return
end
subroutine r8int_to_r8int ( rmin, rmax, r, r2min, r2max, r2 )

!*****************************************************************************80
!
!! R8INT_TO_R8INT maps one R8INT to another.
!
!  Discussion:
!
!    The formula used is
!
!      R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the first range.
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, real ( kind = 8 ) R2MAX, R2MIN, the second range.
!
!    Output, real ( kind = 8 ) R2, the corresponding value in
!    the range [R2MIN,R2MAX].
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin
  real ( kind = 8 ) r2
  real ( kind = 8 ) r2max
  real ( kind = 8 ) r2min

  if ( rmax == rmin ) then

    r2 = ( r2max + r2min ) / 2.0D+00

  else

    r2 = ( ( ( rmax - r        ) * r2min   &
           + (        r - rmin ) * r2max ) &
           / ( rmax     - rmin ) )

  end if

  return
end
subroutine r8int_to_i4int ( rmin, rmax, r, imin, imax, i )

!*****************************************************************************80
!
!! R8INT_TO_I4INT maps an R8INT to an integer interval.
!
!  Discussion:
!
!    The formula used is
!
!      I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, integer ( kind = 4 ) IMAX, IMIN, the integer range.
!
!    Output, integer ( kind = 4 ) I, the corresponding value in the
!    range [IMIN,IMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin

  if ( rmax == rmin ) then

    i = ( imax + imin ) / 2

  else

    i = nint ( &
      ( ( rmax - r        ) * real ( imin, kind = 8 )   &
      + (        r - rmin ) * real ( imax, kind = 8 ) ) &
      / ( rmax     - rmin ) )

  end if

  return
end
subroutine r8mat_add ( m, n, alpha, a, beta, b, c )

!*****************************************************************************80
!
!! R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) ALPHA, the multiplier for A.
!
!    Input, real ( kind = 8 ) A(M,N), the first matrix.
!
!    Input, real ( kind = 8 ) BETA, the multiplier for A.
!
!    Input, real ( kind = 8 ) B(M,N), the second matrix.
!
!    Output, real ( kind = 8 ) C(M,N), the sum of alpha*A+beta*B.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(m,n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(m,n)

  c(1:m,1:n) = alpha * a(1:m,1:n) + beta * b(1:m,1:n)

  return
end
function r8mat_amax ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_AMAX returns the maximum absolute value entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 8 ) R8MAT_AMAX, the maximum absolute value 
!    entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_amax

  r8mat_amax = maxval ( abs ( a(1:m,1:n) ) )

  return
end
subroutine r8mat_border_add ( m, n, table, table2 )

!*****************************************************************************80
!
!! R8MAT_BORDER_ADD adds a "border" to an R8MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    in the interior of a 2D grid, and we wish to create a new table
!    with additional positions for the nodes that would be on the
!    border of the 2D grid.
!
!                  0 0 0 0 0 0
!      * * * *     0 * * * * 0
!      * * * * --> 0 * * * * 0
!      * * * *     0 * * * * 0
!                  0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 3 by 4 array
!    is input, and a 5 by 6 array is to be output.
!
!    The old data is shifted to its correct positions in the new array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
!    Output, real ( kind = 8 ) TABLE2(M+2,N+2), the augmented table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) table2(m+2,n+2)

  table2(1,1:n+2) = 0.0D+00
  table2(m+2,1:n+2) = 0.0D+00
  table2(2:m+1,1) = 0.0D+00
  table2(2:m+1,n+2) = 0.0D+00

  table2(2:m+1,2:n+1) = table(1:m,1:n)

  return
end
subroutine r8mat_border_cut ( m, n, table, table2 )

!*****************************************************************************80
!
!! R8MAT_BORDER_CUT cuts the "border" of an R8MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    on a 2D grid, and we wish to create a new table corresponding only
!    to those nodes in the interior of the 2D grid.
!
!      0 0 0 0 0 0
!      0 * * * * 0    * * * *
!      0 * * * * 0 -> * * * *
!      0 * * * * 0    * * * *
!      0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 5 by 6 array
!    is input, and a 3 by 4 array is to be output.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
!    Output, real ( kind = 8 ) TABLE2(M-2,N-2), the new table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) table2(m-2,n-2)

  if ( m <= 2 .or. n <= 2 ) then
    return
  end if

  table2(1:m-2,1:n-2) = table(2:m-1,2:n-1)

  return
end
subroutine r8mat_cholesky_factor ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) sum2

  flag = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end
subroutine r8mat_cholesky_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N Cholesky factor of the
!    system matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)
!
!  Solve L * y = b.
!
  call r8mat_l_solve ( n, a, b, x )
!
!  Solve L' * x = y.
!
  call r8mat_lt_solve ( n, a, x, x )

  return
end
subroutine r8mat_choresky_factor ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHORESKY_FACTOR computes the "Choresky" factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is an upper triangular matrix R such that:
!
!      A = R * R'
!
!    Note that the usual Cholesky factor is a LOWER triangular matrix L
!    such that
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the N by N upper triangular
!    "Choresky" factor.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) sum2

  flag = 0

  c(n:1:-1,n:1:-1) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

  c(n:1:-1,n:1:-1) = c(1:n,1:n)

  return
end
subroutine r8mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! R8MAT_COPY copies an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be copied.
!
!    Output, real ( kind = 8 ) B(M,N), a copy of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET computes the determinant of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Helmut Spaeth.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)
  real ( kind = 8 ) t

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0D+00

  do k = 1, n

    piv = maxloc ( abs ( b(k:n,k) ) )

    m = piv(1) + k - 1

    if ( m /= k ) then
      det = - det
      t      = b(m,k)
      b(m,k) = b(k,k)
      b(k,k) = t
    end if

    det = det * b(k,k)

    if ( b(k,k) /= 0.0D+00 ) then

      b(k+1:n,k) = -b(k+1:n,k) / b(k,k)

      do j = k + 1, n
        if ( m /= k ) then
          t      = b(m,j)
          b(m,j) = b(k,j)
          b(k,j) = t
        end if
        b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
      end do

    end if

  end do

  return
end
function r8mat_det_2d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The formula for the determinant of a 2 by 2 matrix is
!
!      a11 * a22 - a12 * a21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_2D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) r8mat_det_2d

  r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function r8mat_det_3d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The formula for the determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_3D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) r8mat_det_3d

  r8mat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d

  r8mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
function r8mat_det_5d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(5,5), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_5D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(5,5)
  real ( kind = 8 ) b(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) r8mat_det_5d
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  r8mat_det_5d = 0.0D+00

  do k = 1, 5

    do i = 1, 4
      do j = 1, 4

        if ( j < k ) then
          inc = 0
        else
          inc = 1
        end if

        b(i,j) = a(i+1,j+inc)

      end do
    end do

    r8mat_det_5d = r8mat_det_5d + (-1)**( k + 1 ) * a(1,k) * r8mat_det_4d ( b )

  end do

  return
end
subroutine r8mat_diag_add_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be added to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s

  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine r8mat_diag_add_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector to be added to the diagonal of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)

  do i = 1, n
    a(i,i) = a(i,i) + v(i)
  end do

  return
end
subroutine r8mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) V(N), the diagonal entries
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine r8mat_diag_set_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be assigned to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s

  do i = 1, n
    a(i,i) = s
  end do

  return
end
subroutine r8mat_diag_set_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector to be assigned to the
!    diagonal of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)

  do i = 1, n
    a(i,i) = v(i)
  end do

  return
end
subroutine r8mat_expand_linear ( m, n, x, mfat, nfat, xfat )

!*****************************************************************************80
!
!! R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In this routine, the expansion is specified by giving the number
!    of intermediate values to generate between each pair of original
!    data rows and columns.
!
!    The interpolation is not actually linear.  It uses the functions
!
!      1, x, y, and xy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    input data.
!
!    Input, real ( kind = 8 ) X(M,N), the original data.
!
!    Input, integer ( kind = 4 ) MFAT, NFAT, the number of data values
!    to interpolate between each row, and each column, of original data values.
!
!    Output, real ( kind = 8 ) XFAT(M2,N2), the fattened data, where
!    M2 = (M-1)*(MFAT+1)+1,
!    N2 = (N-1)*(NFAT+1)+1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mfat
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jp1
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) xfat((m-1)*(mfat+1)+1,(n-1)*(nfat+1)+1)

  do i = 1, m

    if ( i < m ) then
      ihi = mfat
    else
      ihi = 0
    end if

    do j = 1, n

      if ( j < n ) then
        jhi = nfat
      else
        jhi = 0
      end if

      if ( i < m ) then
        ip1 = i + 1
      else
        ip1 = i
      end if

      if ( j < n ) then
        jp1 = j + 1
      else
        jp1 = j
      end if

      x00 = x(i,j)
      x10 = x(ip1,j)
      x01 = x(i,jp1)
      x11 = x(ip1,jp1)

      do ii = 0, ihi

        s = real ( ii, kind = 8 ) &
          / real ( ihi + 1, kind = 8 )

        do jj = 0, jhi

          t = real ( jj, kind = 8 ) &
            / real ( jhi + 1, kind = 8 )

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj

          xfat(iii,jjj) = &
                                            x00   &
              + s     * (       x10       - x00 ) &
              + t     * (             x01 - x00 ) &
              + s * t * ( x11 - x10 - x01 + x00 )

        end do

      end do

    end do

  end do

  return
end
subroutine r8mat_expand_linear2 ( m, n, a, m2, n2, a2 )

!*****************************************************************************80
!
!! R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In this version of the routine, the expansion is indicated
!    by specifying the dimensions of the expanded array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), a "small" M by N array.
!
!    Input, integer ( kind = 4 ) M2, N2, the number of rows and columns in A2.
!
!    Output, real ( kind = 8 ) A2(M2,N2), the expanded array, which
!    contains an interpolated version of the data in A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a2(m2,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) s
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2

  do i = 1, m2

    if ( m2 == 1 ) then
      r = 0.5D+00
    else
      r = real ( i - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
    end if

    i1 = 1 + int ( r * real ( m - 1, kind = 8 ) )
    i2 = i1 + 1

    if ( m < i2 ) then
      i1 = m - 1
      i2 = m
    end if

    r1 = real ( i1 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    r2 = real ( i2 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    do j = 1, n2

      if ( n2 == 1 ) then
        s = 0.5D+00
      else
        s = real ( j - 1, kind = 8 ) &
          / real ( n2 - 1, kind = 8 )
      end if

      j1 = 1 + int ( s * real ( n - 1, kind = 8 ) )
      j2 = j1 + 1

      if ( n < j2 ) then
        j1 = n - 1
        j2 = n
      end if

      s1 = real ( j1 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      s2 = real ( j2 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      a2(i,j) = &
        ( ( r2 - r ) * ( s2 - s ) * a(i1,j1) &
        + ( r - r1 ) * ( s2 - s ) * a(i2,j1) &
        + ( r2 - r ) * ( s - s1 ) * a(i1,j2) &
        + ( r - r1 ) * ( s - s1 ) * a(i2,j2) ) &
        / ( ( r2 - r1 ) * ( s2 - s1 ) )

    end do

  end do

  return
end
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8MAT_FS factors and solves a system with one right hand side.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine differs from R8MAT_FSS in two ways:
!    * only one right hand side is allowed;
!    * the input matrix A is not modified.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) temp

  a2(1:n,1:n) = a(1:n,1:n)

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a2(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a2(i,jcol) ) ) then
        piv = abs ( a2(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a2(jcol,1:n)
      a2(jcol,1:n) = a2(ipiv,1:n)
      a2(ipiv,1:n) = row(1:n)

      t       = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = t

    end if
!
!  Scale the pivot row.
!
    a2(jcol,jcol+1:n) = a2(jcol,jcol+1:n) / a2(jcol,jcol)
    b(jcol) = b(jcol) / a2(jcol,jcol)
    a2(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a2(i,jcol) /= 0.0D+00 ) then
        temp = - a2(i,jcol)
        a2(i,jcol) = 0.0D+00
        a2(i,jcol+1:n) = a2(i,jcol+1:n) + temp * a2(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a2(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine r8mat_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! R8MAT_FSS factors and solves a system with multiple right hand sides.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, real ( kind = 8 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t(nb)
  real ( kind = 8 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
subroutine r8mat_givens_post ( n, a, row, col, g )

!*****************************************************************************80
!
!! R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Givens post-multiplier matrix G(ROW,COL) has the property that
!    the (ROW,COL)-th entry of A*G is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of A*G which is to be zeroed out.
!
!    Output, real ( kind = 8 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) row
  real ( kind = 8 ) theta

  call r8mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(row,row) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r8mat_givens_pre ( n, a, row, col, g )

!*****************************************************************************80
!
!! R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Givens premultiplier rotation matrix G(ROW,COL) has the
!    property that the (ROW,COL)-th entry of G*A is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of the G*A which is to be zeroed out.
!
!    Output, real ( kind = 8 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) row
  real ( kind = 8 ) theta

  call r8mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(col,col) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r8mat_hess ( fx, n, x, h )

!*****************************************************************************80
!
!! R8MAT_HESS approximates a Hessian matrix via finite differences.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    H(I,J) = d2 F / d X(I) d X(J)
!
!    The values returned by this routine will be only approximate.
!    In some cases, they will be so poor that they are useless.
!    However, one of the best applications of this routine is for
!    checking your own Hessian calculations, since as Heraclitus
!    said, you'll never get the same result twice when you differentiate
!    a complicated expression by hand.
!
!    The user function routine, here called "FX", should have the form:
!
!      subroutine fx ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FX, the name of the user function routine.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the approximated N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) f00
  real ( kind = 8 ) fmm
  real ( kind = 8 ) fmp
  real ( kind = 8 ) fpm
  real ( kind = 8 ) fpp
  external fx
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi
  real ( kind = 8 ) xj
!
!  Choose the stepsizes.
!
  eps = ( epsilon ( eps ) )**0.33D+00

  do i = 1, n
    s(i) = eps * max ( abs ( x(i) ), 1.0D+00 )
  end do
!
!  Calculate the diagonal elements.
!
  do i = 1, n

    xi = x(i)

    call fx ( n, x, f00 )

    x(i) = xi + s(i)
    call fx ( n, x, fpp )

    x(i) = xi - s(i)
    call fx ( n, x, fmm )

    h(i,i) = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s(i)**2

    x(i) = xi

  end do
!
!  Calculate the off diagonal elements.
!
  do i = 1, n

    xi = x(i)

    do j = i + 1, n

      xj = x(j)

      x(i) = xi + s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fpp )

      x(i) = xi + s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fpm )

      x(i) = xi - s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fmp )

      x(i) = xi - s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fmm )

      h(j,i) = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0D+00 * s(i) * s(j) )

      h(i,j) = h(j,i)

      x(j) = xj

    end do

    x(i) = xi

  end do

  return
end
subroutine r8mat_house_axh ( n, a, v, ah )

!*****************************************************************************80
!
!! R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be postmultiplied.
!
!    Input, real ( kind = 8 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 8 ) AH(N,N), the product A*H.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ah(n,n)
  real ( kind = 8 ) ah_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0D+00 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine r8mat_house_form ( n, v, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    H(v) = I - 2 * v * v' / ( v' * v )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector defining the Householder matrix.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) beta
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(n)
!
!  Compute the L2 norm of V.
!
  beta = sum ( v(1:n)**2 )
!
!  Form the matrix H.
!
  call r8mat_identity ( n, h )

  do i = 1, n
    do j = 1, n
      h(i,j) = h(i,j) - 2.0D+00 * v(i) * v(j) / beta
    end do
  end do

  return
end
subroutine r8mat_house_hxa ( n, a, v, ha )

!*****************************************************************************80
!
!! R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be premultiplied.
!
!    Input, real ( kind = 8 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 8 ) HA(N,N), the product H*A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ha(n,n)
  real ( kind = 8 ) ha_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ha_temp(i,j) = a(i,j)
      do k = 1, n
        ha_temp(i,j) = ha_temp(i,j) - 2.0D+00 * v(i) * v(k) * a(k,j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into HA.
!  Doing it this way means the user can identify the input arguments A and HA.
!
  ha(1:n,1:n) = ha_temp(1:n,1:n)

  return
end
subroutine r8mat_house_post ( n, a, row, col, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    H(ROW,COL) has the property that the ROW-th column of
!    A*H(ROW,COL) is zero from entry COL+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose Householder matrix
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same row, but higher column, will be zeroed out if
!    A is postmultiplied by H.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) row
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(n)
!
!  Set up the vector V.
!
  w(1:col-1) = 0.0D+00
  w(col:n) = a(row,col:n)

  call r8vec_house_column ( n, w, col, v )
!
!  Form the matrix H(V).
!
  call r8mat_house_form ( n, v, h )

  return
end
subroutine r8mat_house_pre ( n, a, row, col, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    H(ROW,COL) has the property that the COL-th column of
!    H(ROW,COL)*A is zero from entry ROW+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose Householder matrix
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same column, but higher rows, will be zeroed out if A is
!    premultiplied by H.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) row
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(n)
!
!  Set up the vector V.
!
  w(1:row-1) = 0.0D+00
  w(row:n) = a(row:n,col)

  call r8vec_house_column ( n, w, row, v )
!
!  Form the matrix H(V).
!
  call r8mat_house_form ( n, v, h )

  return
end
subroutine r8mat_identity ( n, a )

!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  return
end
function r8mat_in_01 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, logical R8MAT_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  logical r8mat_in_01

  if ( any ( a(1:m,1:n) < 0.0D+00 .or. 1.0D+00 < a(1:m,1:n) ) ) then
    r8mat_in_01 = .false.
  else
    r8mat_in_01 = .true.
  end if

  return
end
subroutine r8mat_indicator ( m, n, table )

!*****************************************************************************80
!
!! R8MAT_INDICATOR sets up an "indicator" R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The value of each entry suggests its location, as in:
!
!      11  12  13  14
!      21  22  23  24
!      31  32  33  34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  real ( kind = 8 ) table(m,n)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      table(i,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
function r8mat_insignificant ( m, n, r, s )

!*****************************************************************************80
!
!! R8MAT_INSIGNIFICANT determines if an R8MAT is insignificant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the dimension of the matrices.
!
!    Input, real ( kind = 8 ) R(M,N), the vector to be compared against.
!
!    Input, real ( kind = 8 ) S(M,N), the vector to be compared.
!
!    Output, logical R8MAT_INSIGNIFICANT, is TRUE if S is insignificant
!    compared to R.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(m,n)
  logical r8mat_insignificant
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  logical value

  value = .true.

  do j = 1, n
    do i = 1, m

      t = r(i,j) + s(i,j)
      tol = epsilon ( r(i,j) ) * abs ( r(i,j) )

      if ( tol < abs ( r(i,j) - t ) ) then 
        value = .false.
        exit
      end if

    end do
  end do
  
  r8mat_insignificant = value

  return
end
subroutine r8mat_inverse_2d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) det
  real ( kind = 8 ) r8mat_det_2d
!
!  Compute the determinant of A.
!
  det = r8mat_det_2d ( a )

  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

  else

    b(1,1) =  a(2,2) / det
    b(1,2) = -a(1,2) / det
    b(2,1) = -a(2,1) / det
    b(2,2) =  a(1,1) / det

  end if

  return
end
subroutine r8mat_inverse_3d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_3D inverts a 3 by 3 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(3,3), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) b(3,3)
  real ( kind = 8 ) det
  real ( kind = 8 ) r8mat_det_3d
!
!  Compute the determinant of A.
!
  det = r8mat_det_3d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    b(1:3,1:3) = 0.0D+00
    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine r8mat_inverse_4d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_4D inverts a 4 by 4 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(4,4), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) b(4,4)
  real ( kind = 8 ) det
  real ( kind = 8 ) r8mat_det_4d
!
!  Compute the determinant of A.
!
  det = r8mat_det_4d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:4,1:4) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = +( &
        + a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,1) = -( &
        + a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,1) = +( &
        + a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,1) = -( &
        + a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(2,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,2) = -( &
        + a(1,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(1,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,2) = +( &
        + a(1,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,2) = -( &
        + a(1,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(1,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,2) = +( &
        + a(1,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(1,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(1,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,3) = +( &
        + a(1,2) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,2) - a(2,2) * a(4,4) ) &
        + a(1,4) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        ) / det

  b(2,3) = -( &
        + a(1,1) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,3) - a(2,3) * a(4,1) ) &
        ) / det

  b(3,3) = +( &
        + a(1,1) * ( a(2,2) * a(4,4) - a(2,4) * a(4,2) ) &
        + a(1,2) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(4,3) = -( &
        + a(1,1) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        + a(1,2) * ( a(2,3) * a(4,1) - a(2,1) * a(4,3) ) &
        + a(1,3) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(1,4) = -( &
        + a(1,2) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,2) - a(2,2) * a(3,4) ) &
        + a(1,4) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        ) / det

  b(2,4) = +( &
        + a(1,1) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
        ) / det

  b(3,4) = -( &
        + a(1,1) * ( a(2,2) * a(3,4) - a(2,4) * a(3,2) ) &
        + a(1,2) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  b(4,4) = +( &
        + a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  return
end
subroutine r8mat_is_identity ( n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The routine returns the Frobenius norm of A - I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - I, which would be exactly zero
!    if A were the identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  error_frobenius = 0.0D+00

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        error_frobenius = error_frobenius + ( a(i,j) - 1.0D+00 )**2
      else
        error_frobenius = error_frobenius + a(i,j)**2
      end if
    end do 
  end do

  error_frobenius = sqrt ( error_frobenius )

  return
end
subroutine r8mat_is_nonnegative ( m, n, a, ival )

!*****************************************************************************80
!
!! R8MAT_IS_NONNEGATIVE checks whether an R8MAT is nonnegative.
!
!  Discussion:
!
!    An R8MAT is a matrix of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the row and column dimensions of 
!    the matrix.  M and N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, logical IVAL:
!    TRUE, the matrix is nonnegative.
!    FALSE, at least one element of A is less than 0.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) ival

  ival = all ( 0.0D+00 <= a(1:m,1:n) )

  return
end
subroutine r8mat_is_symmetric ( m, n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
!
!  Discussion:
!
!    An R8MAT is a matrix of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, measures the 
!    Frobenius norm of ( A - A' ), which would be zero if the matrix
!    were exactly symmetric.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) r8_huge

  if ( m /= n ) then
    error_frobenius = r8_huge ( )
    return
  end if

  error_frobenius = sqrt ( &
                      sum ( &
                        ( &
                          abs ( a(1:m,1:n) - transpose ( a(1:m,1:n) ) ) &
                         )**2 &
                       ) &
                     )

  return
end
subroutine r8mat_jac ( m, n, eps, fx, x, fprime )

!*****************************************************************************80
!
!! R8MAT_JAC estimates a dense jacobian matrix of the function FX.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    FPRIME(I,J) = d F(I) / d X(J).
!
!    The jacobian is assumed to be dense, and the LINPACK/LAPACK
!    double precision general matrix storage mode ("DGE") is used.
!
!    Forward differences are used, requiring N+1 function evaluations.
!
!    Values of EPS have typically been chosen between
!    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
!    machine tolerance.
!
!    If EPS is too small, then F(X+EPS) will be the same as
!    F(X), and the jacobian will be full of zero entries.
!
!    If EPS is too large, the finite difference estimate will
!    be inaccurate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) EPS, a tolerance to be used for shifting the
!    X values during the finite differencing.  No single value
!    of EPS will be reliable for all vectors X and functions FX.
!
!    Input, external FX, the name of the user written
!    routine which evaluates the function at a given point X, of the form:
!      subroutine fx ( m, n, x, f )
!      integer m
!      integer n
!      real ( kind = 8 ) f(m)
!      real ( kind = 8 ) x(n)
!      f(1:m) = ...
!      return
!      end
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian
!    is to be estimated.
!
!    Output, real ( kind = 8 ) FPRIME(M,N), the M by N estimated jacobian
!    matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  real ( kind = 8 ) fprime(m,n)
  external fx
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xsave
  real ( kind = 8 ) work1(m)
  real ( kind = 8 ) work2(m)
!
!  Evaluate the function at the base point, X.
!
  call fx ( m, n, x, work2 )
!
!  Now, one by one, vary each component J of the base point X, and
!  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
!
  do j = 1, n

    xsave = x(j)
    del = eps * ( 1.0D+00 + abs ( x(j) ) )
    x(j) = x(j) + del
    call fx ( m, n, x, work1 )
    x(j) = xsave
    fprime(1:m,j) = ( work1(1:m) - work2(1:m) ) / del

  end do

  return
end
subroutine r8mat_l_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_L_INVERSE inverts a lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A lower triangular matrix is a matrix whose only nonzero entries
!    occur on or below the diagonal.
!
!    The inverse of a lower triangular matrix is a lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the lower triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    do i = 1, n

      if ( i < j ) then
        b(i,j) = 0.0D+00
      else if ( j == i ) then
        b(i,j) = 1.0D+00 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r8mat_l_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_L_PRINT prints a lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Example:
!
!    M = 5, N = 5
!    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)
!
!    11
!    21 22
!    31 32 33
!    41 42 43 44
!    51 52 53 54 55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(*), the M by N matrix.  Only the lower
!    triangular elements are stored, in column major order.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(10)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) size
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  jmax = min ( n, m )

  if ( m <= n ) then
    size = ( m * ( m + 1 ) ) / 2
  else if ( n < m ) then
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n
  end if

  if ( all ( a(1:size) == aint ( a(1:size) ) ) ) then

    nn = 10

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a8,10i8)' ) '  Col   ', ( j, j = jlo, jhi )
      write ( *, '(a6)' ) '  Row '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,10i8)' ) i, int ( a(indx(1:jhi+1-jlo)) )
      end do
    end do

  else if ( maxval ( abs ( a(1:size) ) ) < 1000000.0D+00 ) then

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5f14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  else

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5g14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  end if

  return
end
subroutine r8mat_l_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_L_SOLVE solves a lower triangular linear system.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine r8mat_l1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_L1_INVERSE inverts a unit lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of a unit lower triangular matrix is also
!    a unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r8mat_l1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the unit lower triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n

    do j = 1, n

      if ( i < j ) then
        b(i,j) = 0.0D+00
      else if ( j == i ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = -dot_product ( a(i,1:i-1), b(1:i-1,j) )
      end if

    end do
  end do

  return
end
subroutine r8mat_lt_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_LT_SOLVE solves a transposed lower triangular linear system.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    Given the lower triangular matrix A, the linear system to be solved is:
!
!      A' * x = b
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end
subroutine r8mat_lu ( m, n, a, l, p, u )

!*****************************************************************************80
!
!! R8MAT_LU computes the LU factorization of a rectangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The routine is given an M by N matrix A, and produces
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix, and
!      P, an M by M permutation matrix P,
!
!    so that
!
!      A = P' * L * U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix to be factored.
!
!    Output, real ( kind = 8 ) L(M,M), the M by M unit lower triangular factor.
!
!    Output, real ( kind = 8 ) P(M,M), the M by M permutation matrix.
!
!    Output, real ( kind = 8 ) U(M,N), the M by N upper triangular factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  real ( kind = 8 ) l(m,m)
  real ( kind = 8 ) p(m,m)
  real ( kind = 8 ) pivot
  real ( kind = 8 ) u(m,n)

!  Initialize:
!
!    U:=A
!    L:=Identity
!    P:=Identity
!
  u(1:m,1:n) = a(1:m,1:n)

  call r8mat_identity ( m, l )

  p(1:m,1:m) = l(1:m,1:m)
!
!  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
!
  do j = 1, min ( m - 1, n )

    pivot = 0.0D+00
    ipiv = 0

    do i = j, m

      if ( pivot < abs ( u(i,j) ) ) then
        pivot = abs ( u(i,j) )
        ipiv = i
      end if

    end do
!
!  Unless IPIV is zero, swap rows J and IPIV.
!
    if ( ipiv /= 0 ) then

      call r8row_swap ( m, n, u, j, ipiv )

      call r8row_swap ( m, m, l, j, ipiv )

      call r8row_swap ( m, m, p, j, ipiv )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j + 1, m

        if ( u(i,j) /= 0.0D+00 ) then

          l(i,j) = u(i,j) / u(j,j)

          u(i,j) = 0.0D+00

          u(i,j+1:n) = u(i,j+1:n) - l(i,j) * u(j,j+1:n)

        end if

      end do

    end if

  end do

  return
end
function r8mat_max ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAX returns the maximum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAX, the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_max

  r8mat_max = maxval ( a(1:m,1:n) )

  return
end
subroutine r8mat_max_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = -1
  j = -1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(i,j) < a(ii,jj) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r8mat_maxcol_minrow ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    R8MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAXCOL_MINROW, the maximum column
!    minimum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_maxcol_minrow
  real ( kind = 8 ) r8mat_minrow

  r8mat_maxcol_minrow = 0.0D+00

  do i = 1, m

    r8mat_minrow = minval ( a(i,1:n) )

    if ( i == 1 ) then
      r8mat_maxcol_minrow = r8mat_minrow
    else
      r8mat_maxcol_minrow = max ( r8mat_maxcol_minrow, r8mat_minrow )
    end if

  end do

  return
end
function r8mat_maxrow_mincol ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    R8MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAXROW_MINCOL, the maximum row
!    minimum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_maxrow_mincol
  real ( kind = 8 ) r8mat_mincol

  r8mat_maxrow_mincol = 0.0D+00

  do j = 1, n

    r8mat_mincol = minval ( a(1:m,j) )

    if ( j == 1 ) then
      r8mat_maxrow_mincol = r8mat_mincol
    else
      r8mat_maxrow_mincol = max ( r8mat_maxrow_mincol, r8mat_mincol )
    end if

  end do

  return
end
function r8mat_min ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MIN returns the minimum entry of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MIN, the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_min

  r8mat_min = minval ( a(1:m,1:n) )

  return
end
subroutine r8mat_min_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = -1
  j = -1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) < a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r8mat_mincol_maxrow ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    R8MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MINCOL_MAXROW, the minimum column
!    maximum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_mincol_maxrow
  real ( kind = 8 ) r8mat_maxrow

  r8mat_mincol_maxrow = 0.0D+00

  do i = 1, m

    r8mat_maxrow = maxval ( a(i,1:n) )

    if ( i == 1 ) then
      r8mat_mincol_maxrow = r8mat_maxrow
    else
      r8mat_mincol_maxrow = min ( r8mat_mincol_maxrow, r8mat_maxrow )
    end if

  end do

  return
end
function r8mat_minrow_maxcol ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    R8MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MINROW_MAXCOL, the minimum row
!    maximum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_minrow_maxcol
  real ( kind = 8 ) r8mat_maxcol

  r8mat_minrow_maxcol = 0.0D+00

  do j = 1, n

    r8mat_maxcol = maxval ( a(1:m,j) )

    if ( j == 1 ) then
      r8mat_minrow_maxcol = r8mat_maxcol
    else
      r8mat_minrow_maxcol = min ( r8mat_minrow_maxcol, r8mat_maxcol )
    end if

  end do

  return
end
subroutine r8mat_minvm ( n1, n2, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MINVM computes inverse(A) * B for R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N1,N1), B(N1,N2), the matrices.
!
!    Output, real ( kind = 8 ) C(N1,N2), the result, C = inverse(A) * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a(n1,n1)
  real ( kind = 8 ) alu(n1,n1)
  real ( kind = 8 ) b(n1,n2)
  real ( kind = 8 ) c(n1,n2)
  integer ( kind = 4 ) info

  alu(1:n1,1:n1) = a(1:n1,1:n1)
  c(1:n1,1:n2) = b(1:n1,1:n2)

  call r8mat_fss ( n1, alu, n2, c, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_MINVM - Fatal error!'
    write ( *, '(a)' ) '  The matrix A was numerically singular.'
    stop
  end if

  return
end
subroutine r8mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MM multiplies two R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N1,N2), B(N2,N3), the matrices to multiply.
!
!    Output, real ( kind = 8 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real ( kind = 8 ) a(n1,n2)
  real ( kind = 8 ) b(n2,n3)
  real ( kind = 8 ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( a(1:n1,1:n2), b(1:n2,1:n3) )

  return
end
subroutine r8mat_mmt ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MMT computes C = A * B' for two R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = matmul ( A(1:N1,1;N2) ), transpose ( B(1:N3,1:N2) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N1,N2), B(N3,N2), the matrices to multiply.
!
!    Output, real ( kind = 8 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real ( kind = 8 ) a(n1,n2)
  real ( kind = 8 ) b(n3,n2)
  real ( kind = 8 ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( &
                                      a(1:n1,1:n2), &
                          transpose ( b(1:n3,1:n2) ) &
                        )

  return
end
subroutine r8mat_mtm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MTM computes C = A' * B for two R8MAT's.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = matmul ( transpose ( A(1:N2,1;N1) ), B(1:N2,1:N3) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N2,N1), B(N2,N3), the matrices to multiply.
!
!    Output, real ( kind = 8 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real ( kind = 8 ) a(n2,n1)
  real ( kind = 8 ) b(n2,n3)
  real ( kind = 8 ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( transpose ( a(1:n2,1:n1) ), b(1:n2,1:n3) )

  return
end
subroutine r8mat_mtv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R8MAT_MTV multiplies a transposed matrix times a vector
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the product A'*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) y(n)

  y(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8mat_mv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R8MAT_MV multiplies a matrix times a vector.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In FORTRAN90, this operation can be more efficiently carried
!    out by the command
!
!      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  y(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine r8mat_nint ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NINT rounds the entries of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, real ( kind = 8 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)

  a(1:m,1:n) = real ( nint ( a(1:m,1:n) ), kind = 8 )

  return
end
function r8mat_norm_eis ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The EISPACK norm is defined as:
!
!      R8MAT_NORM_EIS =
!        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose EISPACK norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_EIS, the EISPACK norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_norm_eis

  r8mat_norm_eis = sum ( abs ( a(1:m,1:n) ) )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
function r8mat_norm_fro_affine ( m, n, a1, a2 )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A1(M,N), A2(M,N), the matrices for whose 
!    difference the Frobenius norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO_AFFINE, the Frobenius 
!    norm of A1 - A2.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(m,n)
  real ( kind = 8 ) a2(m,n)
  real ( kind = 8 ) r8mat_norm_fro_affine

  r8mat_norm_fro_affine = sqrt ( sum ( ( a1(1:m,1:n) - a2(1:m,1:n) )**2 ) )

  return
end
function r8mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The matrix L1 norm is defined as:
!
!      R8MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L1 norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) col_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_norm_l1

  r8mat_norm_l1 = 0.0D+00

  do j = 1, n
    col_sum = sum ( abs ( a(1:m,j) ) )
    r8mat_norm_l1 = max ( r8mat_norm_l1, col_sum )
  end do

  return
end
function r8mat_norm_l2 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The matrix L2 norm is defined as:
!
!      R8MAT_NORM_L2 = sqrt ( max ( 1 <= I <= M ) LAMBDA(I) )
!
!    where LAMBDA contains the eigenvalues of A * A'.
!
!    The matrix L2 norm is derived from the vector L2 norm, and
!    satisifies:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_l2 ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,m)
  real ( kind = 8 ) diag(m)
  real ( kind = 8 ) r8mat_norm_l2
!
!  Compute B = A * A'.
!
  b(1:m,1:m) = matmul ( a(1:m,1:n), transpose ( a(1:m,1:n) ) )
!
!  Diagonalize B.
!
  call r8mat_symm_jacobi ( m, b )
!
!  Find the maximum eigenvalue, and take its square root.
!
  call r8mat_diag_get_vector ( m, b, diag )

  r8mat_norm_l2 = sqrt ( maxval ( diag(1:m) ) )

  return
end
function r8mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The matrix L-oo norm is defined as:
!
!      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-oo norm is derived from the vector L-oo norm,
!    and satisifies:
!
!      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L-oo
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_norm_li
  real ( kind = 8 ) row_sum

  r8mat_norm_li = 0.0D+00

  do i = 1, m
    row_sum = sum ( abs ( a(i,1:n) ) )
    r8mat_norm_li = max ( r8mat_norm_li, row_sum )
  end do

  return
end
subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )

  return
end
subroutine r8mat_nullspace ( m, n, a, nullspace_size, nullspace )

!*****************************************************************************80
!
!! R8MAT_NULLSPACE computes the nullspace of a matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    Let A be an MxN matrix.
!
!    If X is an N-vector, and A*X = 0, then X is a null vector of A.
!
!    The set of all null vectors of A is called the nullspace of A.
!
!    The 0 vector is always in the null space.
!
!    If the 0 vector is the only vector in the nullspace of A, then A
!    is said to have maximum column rank.  (Because A*X=0 can be regarded
!    as a linear combination of the columns of A).  In particular, if A
!    is square, and has maximum column rank, it is nonsingular.
!
!    The dimension of the nullspace is the number of linearly independent
!    vectors that span the nullspace.  If A has maximum column rank,
!    its nullspace has dimension 0.
!
!    This routine uses the reduced row echelon form of A to determine
!    a set of NULLSPACE_SIZE independent null vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be analyzed.
!
!    Input, integer ( kind = 4 ) NULLSPACE_SIZE, the size of the nullspace.
!
!    Output, real ( kind = 8 ) NULLSPACE(N,NULLSPACE_SIZE), vectors that
!    span the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullspace_size

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  real ( kind = 8 ) nullspace(n,nullspace_size)
  integer ( kind = 4 ) row(m)
  real ( kind = 8 ) rref(m,n)
!
!  Make a copy of A.
!
  rref(1:m,1:n) = a(1:m,1:n)
!
!  Get the reduced row echelon form of A.
!
  call r8mat_rref ( m, n, rref )
!
!  Note in ROW the columns of the leading nonzeros.
!  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
!
  row(1:m) = 0

  do j = 1, n
    col(j) = - j
  end do

  do i = 1, m
    do j = 1, n
      if ( rref(i,j) == 1.0D+00 ) then
        row(i) = j
        col(j) = j
        exit
      end if
    end do
  end do

  nullspace(1:n,1:nullspace_size) = 0.0D+00

  j2 = 0
!
!  If column J does not contain a leading 1, then it contains
!  information about a null vector.
!
  do j = 1, n

    if ( col(j) < 0 ) then

      j2 = j2 + 1

      do i = 1, m
        if ( rref(i,j) /= 0.0D+00 ) then
          i2 = row(i)
          nullspace(i2,j2) = - rref(i,j)
        end if
      end do

      nullspace(j,j2) = 1.0D+00

    end if

  end do

  return
end
subroutine r8mat_nullspace_size ( m, n, a, nullspace_size )

!*****************************************************************************80
!
!! R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    Let A be an MxN matrix.
!
!    If X is an N-vector, and A*X = 0, then X is a null vector of A.
!
!    The set of all null vectors of A is called the nullspace of A.
!
!    The 0 vector is always in the null space.
!
!    If the 0 vector is the only vector in the nullspace of A, then A
!    is said to have maximum column rank.  (Because A*X=0 can be regarded
!    as a linear combination of the columns of A).  In particular, if A
!    is square, and has maximum column rank, it is nonsingular.
!
!    The dimension of the nullspace is the number of linearly independent
!    vectors that span the nullspace.  If A has maximum column rank,
!    its nullspace has dimension 0.
!
!    This routine ESTIMATES the dimension of the nullspace.  Cases of
!    singularity that depend on exact arithmetic will probably be missed.
!
!    The nullspace will be estimated by counting the leading 1's in the
!    reduced row echelon form of A, and subtracting this from N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be analyzed.
!
!    Output, integer ( kind = 4 ) NULLSPACE_SIZE, the estimated size
!    of the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leading
  integer ( kind = 4 ) nullspace_size
  real ( kind = 8 ) rref(m,n)
!
!  Get the reduced row echelon form of A.
!
  rref(1:m,1:n) = a(1:m,1:n)

  call r8mat_rref ( m, n, rref )
!
!  Count the leading 1's in A.
!
  leading = 0
  do i = 1, m
    do j = 1, n
      if ( rref(i,j) == 1.0D+00 ) then
        leading = leading + 1
        exit
      end if
    end do
  end do

  nullspace_size = n - leading

  return
end
subroutine r8mat_orth_uniform ( n, seed, a )

!*****************************************************************************80
!
!! R8MAT_ORTH_UNIFORM returns a random orthogonal R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Properties:
!
!    The inverse of A is equal to A'.
!
!    A * A'  = A' * A = I.
!
!    Columns and rows of A have unit Euclidean norm.
!
!    Distinct pairs of columns of A are orthogonal.
!
!    Distinct pairs of rows of A are orthogonal.
!
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!
!    The determinant of A is +1 or -1.
!
!    All the eigenvalues of A have modulus 1.
!
!    All singular values of A are 1.
!
!    All entries of A are between -1 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pete Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) A(N,N), the orthogonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n - 1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    x(1:j-1) = 0.0D+00

    do i = j, n
      x(i) = r8_normal_01 ( seed )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call r8vec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call r8mat_house_axh ( n, a, v, a )

  end do

  return
end
subroutine r8mat_plot ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PLOT "plots" an R8MAT, with an optional title.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character r8mat_plot_symbol
  character ( len = 70 ) string
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do jlo = 1, n, 70
    jhi = min ( jlo + 70-1, n )
    write ( *, '(a)' ) ' '
    write ( *, '(8x,2x,70i1)' ) ( mod ( j, 10 ), j = jlo, jhi )
    write ( *, '(a)' ) ' '

    do i = 1, m
      do j = jlo, jhi
        string(j+1-jlo:j+1-jlo) = r8mat_plot_symbol ( a(i,j) )
      end do
      write ( *, '(i8,2x,a)' ) i, string(1:jhi+1-jlo)
    end do
  end do

  return
end
function r8mat_plot_symbol ( r )

!*****************************************************************************80
!
!! R8MAT_PLOT_SYMBOL returns a symbol for an element of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, a value whose symbol is desired.
!
!    Output, character R8MAT_PLOT_SYMBOL, is
!    '-' if R is negative,
!    '0' if R is zero,
!    '+' if R is positive.
!
  implicit none

  character r8mat_plot_symbol
  real ( kind = 8 ) r

  if ( r < 0.0D+00 ) then
    r8mat_plot_symbol = '-'
  else if ( r == 0.0D+00 ) then
    r8mat_plot_symbol = '0'
  else if ( 0.0D+00 < r ) then
    r8mat_plot_symbol = '+'
  end if

  return
end
subroutine r8mat_poly_char ( n, a, p )

!*****************************************************************************80
!
!! R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(N) contains the coefficient of X^N
!    (which will be 1), P(I) contains the coefficient of X^I,
!    and P(0) contains the constant term.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) r8mat_trace
  real ( kind = 8 ) trace
  real ( kind = 8 ) work1(n,n)
  real ( kind = 8 ) work2(n,n)
!
!  Initialize WORK1 to the identity matrix.
!
  call r8mat_identity ( n, work1 )

  p(n) = 1.0D+00

  do order = n - 1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    trace = r8mat_trace ( n, work2 )
!
!  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
!
    p(order) = -trace / real ( n - order, kind = 8 )
!
!  WORK1 := WORK2 + P(ORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    do i = 1, n
      work1(i,i) = work1(i,i) + p(order)
    end do

  end do

  return
end
subroutine r8mat_power ( n, a, npow, b )

!*****************************************************************************80
!
!! R8MAT_POWER computes a nonnegative power of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The algorithm is:
!
!      B = I
!      do NPOW times:
!        B = A * B
!      end
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be raised to a power.
!
!    Input, integer ( kind = 4 ) NPOW, the power to which A is to be raised.
!    NPOW must be nonnegative.
!
!    Output, real ( kind = 8 ) B(N,N), the value of A^NPOW.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) npow

  if ( npow < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_POWER - Fatal error!'
    write ( *, '(a)' ) '  Input value of NPOW < 0.'
    write ( *, '(a,i8)' ) '  NPOW = ', npow
    stop
  end if

  call r8mat_identity ( n, b )

  do ipow = 1, npow
    b(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )
  end do

  return
end
subroutine r8mat_power_method ( n, a, r, v )

!*****************************************************************************80
!
!! R8MAT_POWER_METHOD applies the power method to an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the power method has not converged, then calling the routine
!    again immediately with the output from the previous call will
!    continue the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) R, the estimated eigenvalue.
!
!    Input/output, real ( kind = 8 ) V(N), on input, an estimate
!    for the eigenvector.  On output, an improved estimate for the
!    eigenvector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) av(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) it
  real ( kind = 8 ), parameter :: it_eps = 0.0001D+00
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ), parameter :: it_min = 10
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) r_old
  real ( kind = 8 ) v(n)

  eps = sqrt ( epsilon ( 1.0D+00 ) )

  r = sqrt ( sum ( v(1:n)**2 ) )

  if ( r == 0.0D+00 ) then
    v(1:n) = 1.0D+00
    r = sqrt ( real ( n, kind = 8 ) )
  end if

  v(1:n) = v(1:n) / r

  do it = 1, it_max

    av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

    r_old = r
    r = sqrt ( sum ( av(1:n)**2 ) )

    if ( it_min < it ) then
      if ( abs ( r - r_old ) <= it_eps * ( 1.0D+00 + abs ( r ) ) ) then
        exit
      end if
    end if

    v(1:n) = av(1:n)

    if ( r /= 0.0D+00 ) then
      v(1:n) = v(1:n) / r
    end if
!
!  Perturb V a bit, to avoid cases where the initial guess is exactly
!  the eigenvector of a smaller eigenvalue.
!
    if ( it < it_max / 2 ) then
      j = 1 + mod ( it - 1, n )
      v(j) = v(j) + eps * ( 1.0D+00 + abs ( v(j) ) )
      r2 = sqrt ( sum ( v(1:n)**2 ) )
      v(1:n) = v(1:n) / r2
    end if

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_print2 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_PRINT2 prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amin
  integer ( kind = 4 ) i
  character ( len = 10 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical integ
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) npline
  real ( kind = 8 ) r8_log_10
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, m
    do j = 1, n

      if ( integ ) then
        if ( a(i,j) /= real ( int ( a(i,j) ), kind = 8 ) ) then
          integ = .false.
        end if
      end if

    end do
  end do
!
!  Find the maximum and minimum entries.
!
  amax = maxval ( a(1:m,1:n) )
  amin = minval ( a(1:m,1:n) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real matrices.
!
  lmax = int ( r8_log_10 ( amax ) )

  if ( integ ) then
    npline = 79 / ( lmax + 3 )
    write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
  else
    npline = 5
    iform = ' '
  end if
!
!  Print a scalar quantity.
!
  if ( m == 1 .and. n == 1 ) then

    if ( integ ) then
      write ( *, iform ) int ( a(1,1) )
    else
      write ( *, '(2x,g14.6)' ) a(1,1)
    end if
!
!  Column vector of length M,
!
  else if ( n == 1 ) then

    do ilo = 1, m, npline

      ihi = min ( ilo+npline-1, m )

      if ( integ ) then
        write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
      else
        write ( *, '(2x,5g14.6)' ) a(ilo:ihi,1)
      end if

    end do
!
!  Row vector of length N,
!
  else if ( m == 1 ) then

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( integ ) then
        write ( *, iform ) int ( a(1,jlo:jhi) )
      else
        write ( *, '(2x,5g14.6)' ) a(1,jlo:jhi)
      end if

    end do
!
!  M by N Array
!
  else

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( npline < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a,i8)' ) 'Matrix columns ', jlo, ' to ', jhi
        write ( *, '(a)' ) ' '
      end if

      do i = 1, m

        if ( integ ) then
          write ( *, iform ) int ( a(i,jlo:jhi) )
        else
          write ( *, '(2x,5g14.6)' ) a(i,jlo:jhi)
        end if

      end do
    end do

  end if

  return
end
subroutine r8mat_ref ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_REF computes the row echelon form of a matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A matrix is in row echelon form if:
!
!    * The first nonzero entry in each row is 1.
!
!    * The leading 1 in a given row occurs in a column to
!      the right of the leading 1 in the previous row.
!
!    * Rows which are entirely zero must occur last.
!
!  Example:
!
!    Input matrix:
!
!     1.0  3.0  0.0  2.0  6.0  3.0  1.0
!    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
!     3.0  9.0  0.0  0.0  6.0  6.0  2.0
!    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
!
!    Output matrix:
!
!     1.0  3.0  0.0  2.0  6.0  3.0  1.0
!     0.0  0.0  0.0  1.0  2.0  4.5  1.5
!     0.0  0.0  0.0  0.0  0.0  1.0  0.3
!     0.0  0.0  0.0  0.0  0.0  0.0  0.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the matrix to be
!    analyzed.  On output, the REF form of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) r
  real ( kind = 8 ) temp

  lead = 1

  do r = 1, m

    if ( n < lead ) then
      exit
    end if

    i = r

    do while ( a(i,lead) == 0.0D+00 )

      i = i + 1

      if ( m < i ) then
        i = r
        lead = lead + 1
        if ( n < lead ) then
          lead = -1
          exit
        end if
      end if

    end do

    if ( lead < 0 ) then
      exit
    end if

    do j = 1, n
      temp   = a(i,j)
      a(i,j) = a(r,j)
      a(r,j) = temp
    end do

    a(r,1:n) = a(r,1:n) / a(r,lead)

    do i = r + 1, m
      a(i,1:n) = a(i,1:n) - a(i,lead) * a(r,1:n)
    end do

    lead = lead + 1

  end do

  return
end
function r8mat_rms ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RMS returns the RMS norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The matrix RMS norm is defined as:
!
!      R8MAT_RMS = sqrt ( 
!        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) A(I,J)^2 / M / N ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the dimensions of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_RMS, the RMS norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_rms

  r8mat_rms = sqrt ( sum ( a(1:m,1:n)**2 ) / m / n )

  return
end
subroutine r8mat_rref ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RREF computes the reduced row echelon form of a matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A matrix is in row echelon form if:
!
!    * The first nonzero entry in each row is 1.
!
!    * The leading 1 in a given row occurs in a column to
!      the right of the leading 1 in the previous row.
!
!    * Rows which are entirely zero must occur last.
!
!    The matrix is in reduced row echelon form if, in addition to
!    the first three conditions, it also satisfies:
!
!    * Each column containing a leading 1 has no other nonzero entries.
!
!  Example:
!
!    Input matrix:
!
!     1.0  3.0  0.0  2.0  6.0  3.0  1.0
!    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
!     3.0  9.0  0.0  0.0  6.0  6.0  2.0
!    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
!
!    Output matrix:
!
!     1.0  3.0  0.0  0.0  2.0  0.0  0.0
!     0.0  0.0  0.0  1.0  2.0  0.0  0.0
!     0.0  0.0  0.0  0.0  0.0  1.0  0.3
!     0.0  0.0  0.0  0.0  0.0  0.0  0.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the matrix to be
!    analyzed.  On output, the RREF form of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) r
  real ( kind = 8 ) temp

  lead = 1

  do r = 1, m

    if ( n < lead ) then
      exit
    end if

    i = r

    do while ( a(i,lead) == 0.0D+00 )

      i = i + 1

      if ( m < i ) then
        i = r
        lead = lead + 1
        if ( n < lead ) then
          lead = -1
          exit
        end if
      end if

    end do

    if ( lead < 0 ) then
      exit
    end if

    do j = 1, n
      temp   = a(i,j)
      a(i,j) = a(r,j)
      a(r,j) = temp
    end do

    a(r,1:n) = a(r,1:n) / a(r,lead)

    do i = 1, m
      if ( i /= r ) then
        a(i,1:n) = a(i,1:n) - a(i,lead) * a(r,1:n)
      end if
    end do

    lead = lead + 1

  end do

  return
end
subroutine r8mat_scale ( m, n, s, a )

!*****************************************************************************80
!
!! R8MAT_SCALE multiplies an R8MAT by a scalar.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) S, the scale factor.
!
!    Input/output, real ( kind = 8 ) A(M,N), the matrix to be scaled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) s

  a(1:m,1:n) = a(1:m,1:n) * s

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  real ( kind = 8 ) t(n+rhs_num)

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  The pivot row moves into the J-th row.
!
    if ( ipivot /= j ) then
      t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
      a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
      a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
    end if
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then
        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
      end if

    end do

  end do

  return
end
subroutine r8mat_solve_2d ( a, b, det, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant DET is returned as zero, then the matrix A is
!    singular, and does not have an inverse.  In that case, X is
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix.
!
!    Input, real ( kind = 8 ) B(2), the right hand side.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 8 ) X(2), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ) det
  real ( kind = 8 ) x(2)
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    x(1:2) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  x(1) = (  a(2,2) * b(1) - a(1,2) * b(2) ) / det
  x(2) = ( -a(2,1) * b(1) + a(1,1) * b(2) ) / det

  return
end
subroutine r8mat_solve_3d ( a, b, det, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant DET is returned as zero, then the matrix A is
!    singular, and does not have an inverse.  In that case, X is
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix.
!
!    Input, real ( kind = 8 ) B(3), the right hand side.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 8 ) X(3), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) det
  real ( kind = 8 ) x(3)
!
!  Compute the determinant.
!
  det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    x(1:3) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  x(1) = (   ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) * b(1) &
           - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) * b(2) &
           + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) * b(3) ) / det

  x(2) = ( - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) * b(1) &
           + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) * b(2) &
           - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) * b(3) ) / det

  x(3) = (   ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) * b(1) &
           - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) * b(2) &
           + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) * b(3) ) / det

  return
end
subroutine r8mat_solve2 ( n, a, b, x, ierror )

!*****************************************************************************80
!
!! R8MAT_SOLVE2 computes the solution of an N by N linear system.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) ipiv(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  ierror = 0

  ipiv(1:n) = 0
  x(1:n) = 0.0D+00
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!      IMAX has not already been used as a pivot;
!      A(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0D+00
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( amax < abs ( a(i,k) ) ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0D+00

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0D+00
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0D+00

      if ( b(j) == 0.0D+00 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Warning:'
        write ( *, '(a,i8)' ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Error:'
        write ( *, '(a,i8)' ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

  end do

  return
end
function r8mat_sum ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_SUM returns the sum of the entries of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      R8MAT_SUM ( M, N, A ) = SUM ( A(1:M,1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8MAT_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_sum

  r8mat_sum = sum ( a(1:m,1:n) )

  return
end
subroutine r8mat_symm_eigen ( n, x, q, a )

!*****************************************************************************80
!
!! R8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The user must supply the desired eigenvalue vector, and the desired
!    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
!    suitable random orthogonal matrix can be generated by R8MAT_ORTH_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) X(N), the desired eigenvalues for the matrix.
!
!    Input, real ( kind = 8 ) Q(N,N), the eigenvector matrix of A.
!
!    Output, real ( kind = 8 ) A(N,N), a symmetric matrix with
!    eigenvalues X and eigenvectors the columns of Q.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) q(n,n)
  real ( kind = 8 ) x(n)
!
!  Set A = Q * Lambda * Q'.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n
      do k = 1, n
        a(i,j) = a(i,j) + q(i,k) * x(k) * q(j,k)
      end do
    end do
  end do

  return
end
subroutine r8mat_symm_jacobi ( n, a )

!*****************************************************************************80
!
!! R8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This code was modified so that it treats as zero the off-diagonal
!    elements that are sufficiently close to, but not exactly, zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, real ( kind = 8 ) A(N,N), a symmetric N by N matrix.
!    On output, the matrix has been overwritten by an approximately
!    diagonal matrix, with the eigenvalues on the diagonal.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ), parameter :: eps = 0.00001D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) norm_fro
  real ( kind = 8 ) s
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u

  norm_fro = r8mat_norm_fro ( n, n, a )

  it = 0

  do

    it = it + 1

    do i = 1, n
      do j = 1, i - 1

        if ( eps * norm_fro < abs ( a(i,j) ) + abs ( a(j,i) ) ) then

          u = ( a(j,j) - a(i,i) ) / ( a(i,j) + a(j,i) )

          t = sign ( 1.0D+00, u ) / ( abs ( u ) + sqrt ( u * u + 1.0D+00 ) )
          c = 1.0D+00 / sqrt ( t * t + 1.0D+00 )
          s = t * c
!
!  A -> A * Q.
!
          do k = 1, n
            t1 = a(i,k)
            t2 = a(j,k)
            a(i,k) = t1 * c - t2 * s
            a(j,k) = t1 * s + t2 * c
          end do
!
!  A -> QT * A
!
          do k = 1, n
            t1 = a(k,i)
            t2 = a(k,j)
            a(k,i) = c * t1 - s * t2
            a(k,j) = s * t1 + c * t2
          end do

        end if
      end do
    end do
!
!  Test the size of the off-diagonal elements.
!
    sum2 = 0.0D+00
    do i = 1, n
      do j = 1, i - 1
        sum2 = sum2 + abs ( a(i,j) )
      end do
    end do

    if ( sum2 <= eps * ( norm_fro + 1.0D+00 ) ) then
      exit
    end if

    if ( it_max <= it ) then
      exit
    end if

  end do

  return
end
subroutine r8mat_to_r8plu ( n, a, pivot, lu, info )

!*****************************************************************************80
!
!! R8MAT_TO_R8PLU factors a general R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be factored.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, real ( kind = 8 ) LU(N,N), an upper triangular matrix U and
!    the multipliers L which were used to obtain it.  The factorization
!    can be written A = L * U, where L is a product of permutation and
!    unit lower triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) lu(n,n)
  real ( kind = 8 ) temp

  lu(1:n,1:n) = a(1:n,1:n)

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( lu(l,k) ) < abs ( lu(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( lu(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_TO_R8PLU - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      temp    = lu(l,k)
      lu(l,k) = lu(k,k)
      lu(k,k) = temp
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    lu(k+1:n,k) = -lu(k+1:n,k) / lu(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        temp    = lu(l,j)
        lu(l,j) = lu(k,j)
        lu(k,j) = temp
      end if

      lu(k+1:n,j) = lu(k+1:n,j) + lu(k+1:n,k) * lu(k,j)

    end do

  end do

  pivot(n) = n

  if ( lu(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_TO_R8PLU - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
function r8mat_trace ( n, a )

!*****************************************************************************80
!
!! R8MAT_TRACE computes the trace of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The trace of a square matrix is the sum of the diagonal elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose trace is desired.
!
!    Output, real ( kind = 8 ) R8MAT_TRACE, the trace of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8mat_trace

  r8mat_trace = 0.0D+00
  do i = 1, n
    r8mat_trace = r8mat_trace + a(i,i)
  end do

  return
end
subroutine r8mat_transpose ( m, n, a, at )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE makes a transposed copy of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    FORTRAN90 provides the transpose ( ) function which should be preferred
!    over this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be transposed.
!
!    Output, real ( kind = 8 ) AT(N,M), the matrix to be transposed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) at(n,m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  at = transpose ( a )

  return
end
subroutine r8mat_transpose_in_place ( n, a )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_IN_PLACE transposes an R8MAT in place.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns
!    of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be transposed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t

  do j = 1, n
    do i = 1, j - 1
      t      = a(i,j)
      a(i,j) = a(j,i)
      a(j,i) = t
    end do
  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_u_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_U_INVERSE inverts an upper triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    An upper triangular matrix is a matrix whose only nonzero entries
!    occur on or above the diagonal.
!
!    The inverse of an upper triangular matrix is an upper triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the upper triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0D+00
      else if ( i == j ) then
        b(i,j) = 1.0D+00 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,i+1:j), b(i+1:j,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r8mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of a unit upper triangular matrix is also
!    a unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r8mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt,
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the unit upper triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0D+00
      else if ( i == j ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = -dot_product ( a(i,i+1:j), b(i+1:j,j) )
      end if

    end do
  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8mat_uniform_ab ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
!
!  Discussion:
!
!    A <= R(I,J) <= B.
!
!    An R8MAT is an array of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8mat_uniform_abvec ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.
!
!  Discussion:
!
!    A(I) <= R(I,J) <= B(I)
!
!    An R8MAT is an array of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_ABVEC - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a(i) + ( b(i) - a(i) ) * real ( seed, kind = 8 ) &
        * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_uniform_unit ( m, seed, w )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_UNIT generates a uniformly random unit vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) W(M), a random direction vector,
!    with unit norm.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) norm
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(m)
!
!  Get M values from a standard normal distribution.
!
  call r8vec_normal_01 ( m, seed, w )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( w(1:m)**2 ) )
!
!  Normalize the vector.
!
  w(1:m) = w(1:m) / norm

  return
end
subroutine r8mat_vand2 ( n, x, a )

!*****************************************************************************80
!
!! R8MAT_VAND2 returns the N by N row Vandermonde matrix A.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The row Vandermonde matrix returned by this routine reads "across"
!    rather than down.  In particular, each row begins with a 1, followed by
!    some value X, followed by successive powers of X.
!
!  Formula:
!
!    A(I,J) = X(I)^(J-1)
!
!  Properties:
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    The determinant of A is
!
!      det(A) = product ( 2 <= I <= N ) (
!        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
!
!    The matrix A is generally ill-conditioned.
!
!  Example:
!
!    N = 5, X = (2, 3, 4, 5, 6)
!
!    1 2  4   8   16
!    1 3  9  27   81
!    1 4 16  64  256
!    1 5 25 125  625
!    1 6 36 216 1296
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix desired.
!
!    Input, real ( kind = 8 ) X(N), the values that define A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N row Vandermonde matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n
    do j = 1, n

      if ( j == 1 .and. x(i) == 0.0D+00 ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = x(i)**(j-1)
      end if

    end do
  end do

  return
end
function r8mat_vtmv ( m, n, x, a, y )

!*****************************************************************************80
!
!! R8MAT_VTMV multiplies computes the scalar x' * A * y.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) X(N), the first vector factor.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 8 ) Y(M), the second vector factor.
!
!    Output, real ( kind = 8 ) R8MAT_VTMV, the value of X' * A * Y.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_vtmv
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) y(n)

  r8mat_vtmv = dot_product ( x(1:m), matmul ( a(1:m,1:n), y(1:n) ) )

  return
end
subroutine r8mat_zero ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_ZERO zeroes an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix of zeroes.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)

  a(1:m,1:n) = 0.0D+00

  return
end
subroutine r8plu_det ( n, pivot, lu, det )

!*****************************************************************************80
!
!! R8PLU_DET computes the determinant of an R8PLU matrix.
!
!  Discussion:
!
!    The matrix should have been factored by R8MAT_TO_R8PLU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors computed
!    by R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * lu(i,i)
    if ( pivot(i) /= i ) then
      det = -det
    end if
  end do

  return
end
subroutine r8plu_inverse ( n, pivot, lu, a_inverse )

!*****************************************************************************80
!
!! R8PLU_INVERSE computes the inverse of an R8PLU matrix.
!
!  Discussion:
!
!    The matrix should have been factored by R8MAT_TO_R8PLU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from
!    R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors computed by
!    R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) A_INVERSE(N,N), the inverse of the original
!    matrix A that was factored by R8MAT_TO_R8PLU.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_inverse(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(n)

  a_inverse(1:n,1:n) = lu(1:n,1:n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a_inverse(k,k)     = 1.0D+00 / a_inverse(k,k)
    a_inverse(1:k-1,k) = -a_inverse(1:k-1,k) * a_inverse(k,k)

    do j = k + 1, n

      temp             = a_inverse(k,j)
      a_inverse(k,j)   = 0.0D+00
      a_inverse(1:k,j) = a_inverse(1:k,j) + temp * a_inverse(1:k,k)

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a_inverse(k+1:n,k)
    a_inverse(k+1:n,k) = 0.0D+00

    do j = k + 1, n
      a_inverse(1:n,k) = a_inverse(1:n,k) + a_inverse(1:n,j) * work(j)
    end do

    if ( pivot(k) /= k ) then

      do i = 1, n
        temp                  = a_inverse(i,k)
        a_inverse(i,k)        = a_inverse(i,pivot(k))
        a_inverse(i,pivot(k)) = temp
      end do

    end if

  end do

  return
end
subroutine r8plu_mul ( n, pivot, lu, x, b )

!*****************************************************************************80
!
!! R8PLU_MUL computes A * x using the PLU factors of A.
!
!  Discussion:
!
!    It is assumed that R8MAT_TO_R8PLU has computed the PLU factors of
!    the matrix A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the matrix factors computed by
!    R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)
!
!  Y = U * X.
!
  do j = 1, n
    b(1:j-1) = b(1:j-1) + lu(1:j-1,j) * b(j)
    b(j) = lu(j,j) * b(j)
  end do
!
!  B = PL * Y = PL * U * X = A * x.
!
  do j = n - 1, 1, -1

    b(j+1:n) = b(j+1:n) - lu(j+1:n,j) * b(j)

    k = pivot(j)

    if ( k /= j ) then
      temp = b(k)
      b(k) = b(j)
      b(j) = temp
    end if

  end do

  return
end
subroutine r8plu_sol ( n, pivot, lu, b, x )

!*****************************************************************************80
!
!! R8PLU_SOL solves a linear system A*x=b from the PLU factors.
!
!  Discussion:
!
!    The PLU factors should have been computed by R8MAT_TO_R8PLU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors from R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lu(n,n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)
!
!  Solve PL * Y = B.
!
  x(1:n) = b(1:n)

  do k = 1, n - 1

    j = pivot(k)

    if ( j /= k ) then
      temp = x(j)
      x(j) = x(k)
      x(k) = temp
    end if

    x(k+1:n) = x(k+1:n) + lu(k+1:n,k) * x(k)

  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1
    x(k) = x(k) / lu(k,k)
    x(1:k-1) = x(1:k-1) - lu(1:k-1,k) * x(k)
  end do

  return
end
subroutine r8plu_to_r8mat ( n, pivot, lu, a )

!*****************************************************************************80
!
!! R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the matrix factors computed by
!    R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix whose factors are
!    represented by LU and PIVOT.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp

  a(1:n,1:n) = 0.0D+00
  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  do j = 1, n

    do i = 1, n
      a(1:i-1,j) = a(1:i-1,j) + lu(1:i-1,i) * a(i,j)
      a(i,j) = lu(i,i) * a(i,j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do i = n - 1, 1, -1

      a(i+1:n,j) = a(i+1:n,j) - lu(i+1:n,i) * a(i,j)

      k = pivot(i)

      if ( k /= i ) then
        temp   = a(k,j)
        a(k,j) = a(i,j)
        a(i,j) = temp
      end if

    end do

  end do

  return
end
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_deriv ( n, c, p, cp )

!*****************************************************************************80
!
!! R8POLY_DERIV returns the derivative of a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    Input, integer ( kind = 4 ) P, the order of the derivative.
!    0 means no derivative is taken.
!    1 means first derivative,
!    2 means second derivative and so on.
!    Values of P less than 0 are meaningless.  Values of P greater
!    than N are meaningful, but the code will behave as though the
!    value of P was N+1.
!
!    Output, real ( kind = 8 ) CP(0:N-P), the polynomial coefficients of
!    the derivative.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cp(0:*)
  real ( kind = 8 ) cp_temp(0:n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p

  if ( n < p ) then
    return
  end if

  cp_temp(0:n) = c(0:n)

  do d = 1, p
    do i = 0, n - d
      cp_temp(i) = real ( i + 1, kind = 8 ) * cp_temp(i+1)
    end do
    cp_temp(n-d+1) = 0.0D+00
  end do

  cp(0:n-p) = cp_temp(0:n-p)

  return
end
subroutine r8poly_lagrange_0 ( npol, xpol, xval, wval )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which
!    should be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) WVAL, the value of the Lagrange factor at XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real ( kind = 8 ) wval
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  wval = product ( xval - xpol(1:npol) )

  return
end
subroutine r8poly_lagrange_1 ( npol, xpol, xval, dwdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
!
!  Formula:
!
!    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(XPOL(1:NPOL))(X)
!      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )
!
!    We also have the recursion:
!
!      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
!                    = W(XPOL(1:NPOL-1))(X)
!                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) DWDX, the derivative of W with respect to X.
!
  implicit none

  integer ( kind = 4 ) npol

  real ( kind = 8 ) dwdx
  integer ( kind = 4 ) i
  real ( kind = 8 ) w
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  dwdx = 0.0D+00
  w = 1.0D+00

  do i = 1, npol

    dwdx = w + ( xval - xpol(i) ) * dwdx
    w = w * ( xval - xpol(i) )

  end do

  return
end
subroutine r8poly_lagrange_2 ( npol, xpol, xval, dw2dx2 )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
!
!  Formula:
!
!    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(X) = Sum ( 1 <= J <= NPOL )
!            Product ( I /= J ) ( X - XPOL(I) )
!
!    W"(X) = Sum ( 1 <= K <= NPOL )
!            Sum ( J =/ K )
!            Product ( I /= K, J ) ( X - XPOL(I) )
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) DW2DX2, the second derivative of W
!    with respect to XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real ( kind = 8 ) dw2dx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) term
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  dw2dx2 = 0.0D+00

  do k = 1, npol

    do j = 1, npol

      if ( j /= k ) then
        term = 1.0D+00

        do i = 1, npol
          if ( i /= j .and. i /= k ) then
            term = term * ( xval - xpol(i) )
          end if
        end do

        dw2dx2 = dw2dx2 + term

      end if

    end do

  end do

  return
end
subroutine r8poly_lagrange_coef ( npol, ipol, xpol, pcof )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
!
!  Discussion:
!
!    Given distinct abscissas XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!    However sometimes it is desirable to be able to write down
!    the standard polynomial coefficients of L(IPOL)(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer ( kind = 4 ) IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the
!    Lagrange polynomials.  The entries in XPOL must be distinct.
!
!    Output, real ( kind = 8 ) PCOF(0:NPOL-1), the standard polynomial
!    coefficients of the IPOL-th Lagrange polynomial:
!      L(IPOL)(X) = SUM ( 0 <= I <= NPOL-1 ) PCOF(I) * X^I
!
  implicit none

  integer ( kind = 4 ) npol

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) ipol
  integer ( kind = 4 ) j
  real ( kind = 8 ) pcof(0:npol-1)
  logical r8vec_distinct
  real ( kind = 8 ) xpol(npol)
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    write ( *, '(a,i8)' ) '  NPOL = ', npol
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop
  end if

  pcof(0) = 1.0D+00
  pcof(1:npol-1) = 0.0D+00

  indx = 0

  do i = 1, npol

    if ( i /= ipol ) then

      indx = indx + 1

      do j = indx, 0, -1

        pcof(j) = -xpol(i) * pcof(j) / ( xpol(ipol) - xpol(i) )

        if ( 0 < j ) then
          pcof(j) = pcof(j) + pcof(j-1) / ( xpol(ipol) - xpol(i) )
        end if

      end do

    end if

  end do

  return
end
subroutine r8poly_lagrange_factor ( npol, xpol, xval, wval, dwdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    Suppose F(X) is at least N times continuously differentiable in the
!    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
!    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
!    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
!    Then in the interval [A,B], the maximum error
!
!      abs ( F(X) - P(X) )
!
!    is bounded by:
!
!      C * FNMAX * W(X)
!
!    where
!
!      C is a constant,
!      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
!      W(X) is the Lagrange factor.
!
!    Thus, the value of W(X) is useful as part of an estimated bound
!    for the interpolation error.
!
!    Note that the Chebyshev abscissas have the property that they minimize
!    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
!    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
!    other choices.
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) WVAL, the value of the Lagrange factor at XVAL.
!
!    Output, real ( kind = 8 ) DWDX, the derivative of W with respect to XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real ( kind = 8 ) dwdx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) term
  real ( kind = 8 ) wval
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  wval = product ( xval - xpol(1:npol) )

  dwdx = 0.0D+00

  do i = 1, npol

    term = 1.0D+00

    do j = 1, npol
      if ( i /= j ) then
        term = term * ( xval - xpol(j) )
      end if
    end do

    dwdx = dwdx + term

  end do

  return
end
subroutine r8poly_lagrange_val ( npol, ipol, xpol, xval, pval, dpdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
!
!  Discussion:
!
!    Given NPOL distinct abscissas, XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer ( kind = 4 ) IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the Lagrange
!    polynomials.  The entries in XPOL must be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the IPOL-th
!    Lagrange polynomial is to be evaluated.
!
!    Output, real ( kind = 8 ) PVAL, the value of the IPOL-th Lagrange
!    polynomial at XVAL.
!
!    Output, real ( kind = 8 ) DPDX, the derivative of the IPOL-th
!    Lagrange polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real ( kind = 8 ) dpdx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipol
  integer ( kind = 4 ) j
  real ( kind = 8 ) p2
  real ( kind = 8 ) pval
  logical r8vec_distinct
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_VAL - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_VAL - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop
  end if
!
!  Evaluate the polynomial.
!
  pval = 1.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      pval = pval * ( xval - xpol(i) ) / ( xpol(ipol) - xpol(i) )

    end if

  end do
!
!  Evaluate the derivative, which can be found by summing up the result
!  of differentiating one factor at a time, successively.
!
  dpdx = 0.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      p2 = 1.0D+00
      do j = 1, npol

        if ( j == i ) then
          p2 = p2                      / ( xpol(ipol) - xpol(j) )
        else if ( j /= ipol ) then
          p2 = p2 * ( xval - xpol(j) ) / ( xpol(ipol) - xpol(j) )
        end if

      end do

      dpdx = dpdx + p2

    end if

  end do

  return
end
subroutine r8poly_order ( na, a, order )

!*****************************************************************************80
!
!! R8POLY_ORDER returns the order of a polynomial.
!
!  Discussion:
!
!    The order of a polynomial is one more than the degree.
!
!    The order of a constant polynomial is 1.  The order of the
!    zero polynomial is debatable, but this routine returns the
!    order as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) ORDER, the order of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) order

  order = na + 1

  do while ( 1 < order )

    if ( a(order-1) /= 0.0D+00 ) then
      return
    end if

    order = order - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( n2 <= 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_shift ( scale, shift, n, poly_cof )

!*****************************************************************************80
!
!! R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
!
!  Discussion:
!
!    Assuming P(X) is a polynomial in the argument X, of the form:
!
!      P(X) =
!          C(N) * X^N
!        + ...
!        + C(1) * X
!        + C(0),
!
!    and that Z is related to X by the formula:
!
!      Z = SCALE * X + SHIFT
!
!    then this routine computes coefficients C for the polynomial Q(Z):
!
!      Q(Z) =
!          C(N) * Z^N
!        + ...
!        + C(1) * Z
!        + C(0)
!
!    so that:
!
!      Q(Z(X)) = P(X)
!
!  Example:
!
!    P(X) = 2 * X^2 - X + 6
!
!    Z = 2.0 * X + 3.0
!
!    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12
!
!    Q(Z(X)) = 0.5 * ( 4.0 * X^2 + 12.0 * X +  9 )
!            - 3.5 * (              2.0 * X +  3 )
!                                           + 12
!
!            = 2.0         * X^2 -  1.0 * X +  6
!
!            = P(X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 1999
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes: The Art of Scientific Computing,
!    Cambridge University Press.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SHIFT, SCALE, the shift and scale applied to X,
!    so that Z = SCALE * X + SHIFT.
!
!    Input, integer ( kind = 4 ) N, the number of coefficients.
!
!    Input/output, real ( kind = 8 ) POLY_COF(0:N).
!    On input, the coefficient array in terms of the X variable.
!    On output, the coefficient array in terms of the Z variable.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) poly_cof(0:n)
  real ( kind = 8 ) scale
  real ( kind = 8 ) shift

  do i = 1, n
    poly_cof(i:n) = poly_cof(i:n) / scale
  end do

  do i = 0, n - 1
    do j = n - 1, i, -1
      poly_cof(j) = poly_cof(j) - shift * poly_cof(j+1)
    end do
  end do

  return
end
subroutine r8poly_value ( m, c, n, x, p )

!*****************************************************************************80
!
!! R8POLY_VALUE evaluates a polynomial.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c1 + c2 * x + c3 * x^2 + ... + cm * x^(m-1)
!
!    is to be evaluated at the vector of values X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = 8 ) C(0:M), the polynomial coefficients.  
!    C(1) is the constant term.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) P(N), the value of the polynomial at the 
!    evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) x(n)

  p(1:n) = c(m)
  do i = m - 1, 0, -1
    p(1:n) = p(1:n) * x(1:n) + c(i)
  end do

  return
end
subroutine r8poly_value_horner ( n, c, x, cx )

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of C.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cx
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  cx = c(n)
  do i = n - 1, 0, -1
    cx = cx * x + c(i)
  end do

  return
end
function r8poly_value_old ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_VALUE_OLD evaluates an R8POLY.
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE_OLD, the value of the polynomial 
!    at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value_old
  real ( kind = 8 ) x

  r8poly_value_old = a(n)
  do i = n - 1, 1, -1
    r8poly_value_old = r8poly_value_old * x + a(i)
  end do

  return
end
subroutine r8poly_value_2d ( m, c, n, x, y, p )

!*****************************************************************************80
!
!! R8POLY_VALUE_2D evaluates a polynomial in 2 variables, X and Y.
!
!  Discussion:
!
!    We assume the polynomial is of total degree M, and has the form:
!
!      p(x,y) = c00 
!             + c10 * x                + c01 * y
!             + c20 * x^2   + c11 * xy + c02 * y^2
!             + ...
!             + cm0 * x^(m) + ...      + c0m * y^m.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) C(T(M+1)), the polynomial coefficients.  
!    C(1) is the constant term.  T(M+1) is the M+1-th triangular number.
!    The coefficients are stored consistent with the following ordering
!    of monomials: 1, X, Y, X^2, XY, Y^2, X^3, X^2Y, XY^2, Y^3, X^4, ...
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the evaluation points.
!
!    Output, real ( kind = 8 ) P(N), the value of the polynomial at the 
!    evaluation points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(*)
  integer ( kind = 4 ) ex
  integer ( kind = 4 ) ey
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) p(n)
  integer ( kind = 4 ) s
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  p(1:n) = 0.0D+00

  j = 0
  do s = 0, m
    do ex = s, 0, -1
      ey = s - ex
      j = j + 1
      p(1:n) = p(1:n) + c(j) * x(1:n) ** ex * y(1:n) ** ey
    end do
  end do

  return
end
subroutine r8poly2_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )

!*****************************************************************************80
!
!! R8POLY2_EX finds the extremal point of a parabola determined by three points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = 8 ) X, Y, the X coordinate of the extremal point
!    of the parabola, and the value of the parabola at that point.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none

  real ( kind = 8 ) bot
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if

  bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3

  if ( bot == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = 0.5D+00 * ( &
          x1**2 * ( y3 - y2 ) &
        + x2**2 * ( y1 - y3 ) &
        + x3**2 * ( y2 - y1 ) ) / bot

  y = ( &
         ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
       - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
       + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
       ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )

!*****************************************************************************80
!
!! R8POLY2_EX2 finds extremal point of a parabola determined by three points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = 8 ) X, Y, the X coordinate of the extremal
!    point of the parabola, and the value of the parabola at that point.
!
!    Output, real ( kind = 8 ) A, B, C, the coefficients that define the
!    parabola: P(X) = A * X * X + B * X + C.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) det
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) v(3,3)
  real ( kind = 8 ) w(3,3)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if
!
!  Set up the Vandermonde matrix.
!
  v(1,1) = 1.0D+00
  v(1,2) = x1
  v(1,3) = x1 * x1

  v(2,1) = 1.0D+00
  v(2,2) = x2
  v(2,3) = x2 * x2

  v(3,1) = 1.0D+00
  v(3,2) = x3
  v(3,3) = x3 * x3
!
!  Get the inverse.
!
  call r8mat_inverse_3d ( v, w, det )
!
!  Compute the parabolic coefficients.
!
  c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
  b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
  a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
!
!  Determine the extremal point.
!
  if ( a == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = -b / ( 2.0D+00 * a )
  y = a * x * x + b * x + c

  return
end
subroutine r8poly2_root ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R8POLY2_ROOT returns the two roots of a quadratic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X * X + B * X + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, the roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  complex ( kind = 8 ) disc
  complex ( kind = 8 ) q
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_ROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0D+00 * a * c
  q = -0.5D+00 * ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine r8poly2_rroot ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
!
!  Example:
!
!     A    B    C       roots              R1   R2
!    --   --   --     ------------------   --   --
!     1   -4    3     1          3          1    3
!     1    0    4     2*i      - 2*i        0    0
!     2   -6    5     3 +   i    3 -   i    3    3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the coefficients of the quadratic
!    polynomial A * X * X + B * X + C = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, real ( kind = 8 ) R1, R2, the real parts of the roots
!    of the polynomial.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) disc
  real ( kind = 8 ) q
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_RROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0D+00 * a * c
  disc = max ( disc, 0.0D+00 )

  q = ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = -0.5D+00 * q / a
  r2 = -2.0D+00 * c / q

  return
end
subroutine r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )

!*****************************************************************************80
!
!! R8POLY2_VAL evaluates a parabola defined by three data values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, three pairs of data.
!    If the X values are distinct, then all the Y values represent
!    actual values of the parabola.
!
!    Three special cases are allowed:
!
!      X1 == X2 /= X3: Y2 is the derivative at X1;
!      X1 /= X2 == X3: Y3 is the derivative at X3;
!      X1 == X2 == X3: Y2 is the derivative at X1, and
!                      Y3 is the second derivative at X1.
!
!    Input, real ( kind = 8 ) X, an abscissa at which the parabola is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) Y, YP, YPP, the values of the parabola and
!    its first and second derivatives at X.
!
  implicit none

  integer ( kind = 4 ) distinct
  real ( kind = 8 ) dif1
  real ( kind = 8 ) dif2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yp
  real ( kind = 8 ) ypp
!
!  If any X's are equal, put them and the Y data first.
!
  if ( x1 == x2 .and. x2 == x3 ) then
    distinct = 1
  else if ( x1 == x2 ) then
    distinct = 2
  else if ( x1 == x3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL - Fatal error!'
    write ( *, '(a)' ) '  X1 = X3 =/= X2.'
    write ( *, '(a,g14.6)' ) '  X1 = ', x1
    write ( *, '(a,g14.6)' ) '  X2 = ', x2
    write ( *, '(a,g14.6)' ) '  X3 = ', x3
    stop
  else if ( x2 == x3 ) then
    distinct = 2
    call r8_swap ( x1, x2 )
    call r8_swap ( x2, x3 )
    call r8_swap ( y1, y2 )
    call r8_swap ( y2, y3 )
  else
    distinct = 3
  end if
!
!  Set up the coefficients.
!
  if ( distinct == 1 ) then

    dif1 = y2
    dif2 = 0.5D+00 * y3

  else if ( distinct == 2 ) then

    dif1 = y2
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 ) - y2 ) / ( x3 - x2 )

  else if ( distinct == 3 ) then

    dif1 = ( y2 - y1 ) / ( x2 - x1 )
    dif2 =  ( ( y3 - y1 ) / ( x3 - x1 ) &
            - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 )

  end if
!
!  Evaluate.
!
  y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2
  yp = dif1 + ( 2.0D+00 * x - x1 - x2 ) * dif2
  ypp = 2.0D+00 * dif2

  return
end
subroutine r8poly2_val2 ( dim_num, ndata, tdata, ydata, left, tval, yval )

!*****************************************************************************80
!
!! R8POLY2_VAL2 evaluates a parabolic interpolant through tabular data.
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of a single data point.
!    DIM_NUM must be at least 1.
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The values in TDATA must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) YDATA(DIM_NUM,NDATA), the data points
!    corresponding to the abscissas.
!
!    Input, integer ( kind = 4 ) LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real ( kind = 8 ) TVAL, the value of T at which the parabolic
!    interpolant is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA),
!    and the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real ( kind = 8 ) YVAL(DIM_NUM), the value of the parabolic
!    interpolant at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) dif1
  real ( kind = 8 ) dif2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) tval
  real ( kind = 8 ) tdata(ndata)
  real ( kind = 8 ) ydata(dim_num,ndata)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yval(dim_num)
!
!  Check.
!
  if ( left < 1 .or. ndata-2 < left ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  LEFT < 1 or NDATA-2 < LEFT.'
    write ( *, '(a,i8)' ) '  LEFT = ', left
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t2 <= t1 .or. t3 <= t2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  T2 <= T1 or T3 <= T2.'
    write ( *, '(a,g14.6)' ) '  T1 = ', t1
    write ( *, '(a,g14.6)' ) '  T2 = ', t2
    write ( *, '(a,g14.6)' ) '  T3 = ', t3
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, dim_num

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
           - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine r8poly3_root ( a, b, c, d, r1, r2, r3 )

!*****************************************************************************80
!
!! R8POLY3_ROOT returns the three roots of a cubic polynomial.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A * X^3 + B * X^2 + C * X + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  complex ( kind = 8 ) i
  complex ( kind = 8 ) one
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) temp
  real ( kind = 8 ) theta

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY3_ROOT - Fatal error!'
    write ( *, '(a)' ) '  A must not be zero!'
    stop
  end if

  one = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  i = sqrt ( -one )

  q = ( ( b / a )**2 - 3.0D+00 * ( c / a ) ) / 9.0D+00

  r = ( 2.0D+00 * ( b / a )**3 - 9.0D+00 * ( b / a ) * ( c / a ) &
      + 27.0D+00 * ( d / a ) ) / 54.0D+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q**3 ) )
    r1 = -2.0D+00 * sqrt ( q ) * cos (   theta                  / 3.0D+00 )
    r2 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 2.0D+00 * pi ) / 3.0D+00 )
    r3 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 4.0D+00 * pi ) / 3.0D+00 )

  else if ( q * q * q <= r * r ) then

    temp = -r + sqrt ( r**2 - q**3 )
    s1 = sign ( 1.0D+00, temp ) * ( abs ( temp ) )**(1.0D+00/3.0D+00)

    temp = -r - sqrt ( r**2 - q**3 )
    s2 = sign ( 1.0D+00, temp ) * ( abs ( temp ) )**(1.0D+00/3.0D+00)

    r1 = s1 + s2
    r2 = -0.5D+00 * ( s1 + s2 ) + i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )
    r3 = -0.5D+00 * ( s1 + s2 ) - i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0D+00 * a )
  r2 = r2 - b / ( 3.0D+00 * a )
  r3 = r3 - b / ( 3.0D+00 * a )

  return
end
subroutine r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

!*****************************************************************************80
!
!! R8POLY4_ROOT returns the four roots of a quartic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, R3, R4, the roots of the polynomial.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a3
  real ( kind = 8 ) a4
  real ( kind = 8 ) b
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) c
  real ( kind = 8 ) c3
  real ( kind = 8 ) c4
  real ( kind = 8 ) d
  real ( kind = 8 ) d3
  real ( kind = 8 ) d4
  real ( kind = 8 ) e
  complex ( kind = 8 ) p
  complex ( kind = 8 ) q
  complex ( kind = 8 ) r
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r4
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY4_ROOT - Fatal error!'
    write ( *, '(a)') '  A must not be zero!'
    stop
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0D+00
  b3 = -b4
  c3 = a4 * c4 - 4.0D+00 * d4
  d3 = -a4 * a4 * d4 + 4.0D+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25D+00 * A4**2 - B4 + R1 )
!
  r = sqrt ( 0.25D+00 * a4**2 - b4 + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75D+00 * a4**2 - r**2 - 2.0D+00 * b4 &
        + 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4**3 ) / r )

    q = sqrt ( 0.75D+00 * a4**2 - r**2 - 2.0D+00 * b4 &
        - 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4**3 ) / r )

  else

    p = sqrt ( 0.75D+00 * a4**2 - 2.0D+00 * b4 &
      + 2.0D+00 * sqrt ( r1**2 - 4.0D+00 * d4 ) )

    q = sqrt ( 0.75D+00 * a4**2 - 2.0D+00 * b4 &
      - 2.0D+00 * sqrt ( r1**2 - 4.0D+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = -0.25D+00 * a4 + 0.5D+00 * r + 0.5D+00 * p
  r2 = -0.25D+00 * a4 + 0.5D+00 * r - 0.5D+00 * p
  r3 = -0.25D+00 * a4 - 0.5D+00 * r + 0.5D+00 * q
  r4 = -0.25D+00 * a4 - 0.5D+00 * r - 0.5D+00 * q

  return
end
function r8r8_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! R8R8_COMPARE compares two R8R8's.
!
!  Discussion:
!
!    An R8R8 is simply a pair of R8 values, stored separately.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, the first vector.
!
!    Input, real ( kind = 8 ) X2, Y2, the second vector.
!
!    Output, integer ( kind = 4 ) R8R8_COMPARE:
!    -1, (X1,Y1) < (X2,Y2);
!     0, (X1,Y1) = (X2,Y2);
!    +1, (X1,Y1) > (X2,Y2).
!
  implicit none

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8_compare
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( x1 < x2 ) then
    compare = -1
  else if ( x2 < x1 ) then
    compare = +1
  else if ( y1 < y2 ) then
    compare = -1
  else if ( y2 < y1 ) then
    compare = +1
  else
    compare = 0
  end if

  r8r8_compare = compare

  return
end
subroutine r8r8_print ( a1, a2, title )

!*****************************************************************************80
!
!! R8R8_PRINT prints an R8R8.
!
!  Discussion:
!
!    An R8R8 is simply a pair of R8R8's, stored separately.
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, A2, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', a1, ',', a2, ')'
  else
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1 )' ) '(', a1, ',', a2, ')'
  end if

  return
end
function r8r8r8_compare ( x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! R8R8R8_COMPARE compares two R8R8R8's.
!
!  Discussion:
!
!    An R8R8R8 is simply 3 R8 values, stored as scalars.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, Z1, the first vector.
!
!    Input, real ( kind = 8 ) X2, Y2, Z2, the second vector.
!
!    Output, integer ( kind = 4 ) R8R8R8_COMPARE:
!    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
!     0, (X1,Y1,Z1) = (X2,Y2,Z2);
!    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
!
  implicit none

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8r8_compare
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  if ( x1 < x2 ) then
    compare = -1
  else if ( x2 < x1 ) then
    compare = +1
  else if ( y1 < y2 ) then
    compare = -1
  else if ( y2 < y1 ) then
    compare = +1
  else if ( z1 < z2 ) then
    compare = -1
  else if ( z2 < z1 ) then
    compare = +1
  else
    compare = 0
  end if

  r8r8r8_compare = compare

  return
end
subroutine r8r8r8vec_index_insert_unique ( n_max, n, x, y, z, indx, &
  xval, yval, zval, ival, ierror )

!*****************************************************************************80
!
!! R8R8R8VEC_INDEX_INSERT_UNIQUE inserts unique R8R8R in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8R8VEC is set of N R8R8R8 items.
!
!    An R8R8R8 is simply 3 R8 values, stored as scalars.
!
!    If the input value does not occur in the current list, it is added,
!    and N, X, Y, Z and INDX are updated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), Z(N), the R8R8R8 vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be inserted
!    if it is not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in X, Y, Z corresponding
!    to the value XVAL, YVAL, ZVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error
!    occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n_max)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(n_max)
  real ( kind = 8 ) zval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    z(1) = zval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
!
  call r8r8r8vec_index_search ( n, x, y, z, indx, xval, yval, zval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    z(n+1) = zval
    ival = n + 1
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = indx(equal)

  end if

  return
end
subroutine r8r8r8vec_index_search ( n, x, y, z, indx, xval, yval, &
  zval, less, equal, more )

!*****************************************************************************80
!
!! R8R8R8VEC_INDEX_SEARCH searches for R8R8R8 value in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8R8VEC is set of N R8R8R8 items.
!
!    An R8R8R8 is simply 3 R8 values, stored as scalars.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the list.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8r8_compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ymid
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(n)
  real ( kind = 8 ) zhi
  real ( kind = 8 ) zlo
  real ( kind = 8 ) zmid
  real ( kind = 8 ) zval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))
  zlo = z(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))
  zhi = z(indx(hi))

  compare = r8r8r8_compare ( xval, yval, zval, xlo, ylo, zlo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = r8r8r8_compare ( xval, yval, zval, xhi, yhi, zhi )

  if ( compare == 1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))
    zmid = z(indx(mid))

    compare = r8r8r8_compare ( xval, yval, zval, xmid, ymid, zmid )

    if ( compare == 0 ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8r8vec_index_insert_unique ( n_max, n, x, y, indx, xval, yval, &
  ival, ierror )

!*****************************************************************************80
!
!! R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8VEC is set of N R8R8 items.
!
!    An R8R8 is simply 2 R8 values, stored as scalars.
!
!    If the input value does not occur in the current list, it is added,
!    and N, X, Y and INDX are updated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), the list of R8R8 vectors.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in X, Y corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an
!    error occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n_max)
  real ( kind = 8 ) yval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in ( X, Y )?
!
  call r8r8vec_index_search ( n, x, y, indx, xval, yval, less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    ival = n + 1
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = indx(equal)

  end if

  return
end
subroutine r8r8vec_index_search ( n, x, y, indx, xval, yval, less, equal, &
  more )

!*****************************************************************************80
!
!! R8R8VEC_INDEX_SEARCH searches for an R8R8 in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8VEC is set of N R8R8 items.
!
!    An R8R8 is simply 2 R8 values, stored as scalars.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8_compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ymid
  real ( kind = 8 ) yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  compare = r8r8_compare ( xval, yval, xlo, ylo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = r8r8_compare ( xval, yval, xhi, yhi )

  if ( compare == 1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    compare = r8r8_compare ( xval, yval, xmid, ymid )

    if ( compare == 0 ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8row_compare ( m, n, a, i, j, value )

!*****************************************************************************80
!
!! R8ROW_COMPARE compares rows in an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
!    each of length N.
!
!  Example:
!
!    Input:
!
!      M = 4, N = 3, I = 2, J = 4
!
!      A = (
!        1  5  9
!        2  6 10
!        3  7 11
!        4  8 12 )
!
!    Output:
!
!      VALUE = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ( kind = 4 ) VALUE, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) value
!
!  Check.
!
  if ( i < 1 .or. m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  value = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      value = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      value = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8row_max ( m, n, a, amax )

!*****************************************************************************80
!
!! R8ROW_MAX returns the maximums of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MAX =
!      3
!      7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMAX(M), the maximums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amax(i) = a(i,1)
    do j = 2, n
      if ( amax(i) < a(i,j) ) then
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! R8ROW_MEAN returns the means of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      2
!      5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) MEAN(M), the means, or averages, of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(m)

  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )
  end do

  return
end
subroutine r8row_min ( m, n, a, amin )

!*****************************************************************************80
!
!! R8ROW_MIN returns the minimums of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MIN =
!      1
!      2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMIN(M), the minimums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_part_quick_a ( m, n, a, l, r )

!*****************************************************************************80
!
!! R8ROW_PART_QUICK_A reorders the rows of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
!    each of length N.
!
!    The routine reorders the rows of A.  Using A(1,1:N) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      M = 8, N = 2
!      A = ( 2 4
!            8 8
!            6 2
!            0 2
!           10 6
!           10 0
!            0 6
!            5 8 )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( 0 2    LEFT
!            0 6
!            ----
!            2 4    KEY
!            ----
!            8 8    RIGHT
!            6 2
!           10 6
!           10 0
!            5 8 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row dimension of A.
!
!    Input, integer ( kind = 4 ) N, the column dimension of A, and the
!    length of a row.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three
!    segments.  Let KEY = the input value of A(1,1:N).  Then
!    I <= L                 A(I,1:N) < KEY;
!         L < I < R         A(I,1:N) = KEY;
!                 R <= I    KEY < A(I,1:N).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) key(n)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    return
  end if

  if ( m == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:n) = a(1,1:n)
  k = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = m + 1

  do j = 2, m

    if ( r8vec_gt ( n, a(l+1,1:n), key(1:n) ) ) then
      r = r - 1
      call r8vec_swap ( n, a(r,1:n), a(l+1,1:n) )
    else if ( r8vec_eq ( n, a(l+1,1:n), key(1:n) ) ) then
      k = k + 1
      call r8vec_swap ( n, a(k,1:n), a(l+1,1:n) )
      l = l + 1
    else if ( r8vec_lt ( n, a(l+1,1:n), key(1:n) ) ) then
      l = l + 1
    end if

  end do
!
!  Shift small elements to the left.
!
  do j = 1, l - k
    a(j,1:n) = a(j+k,1:n)
  end do
!
!  Shift KEY elements to center.
!
  do j = l - k + 1, l
    a(j,1:n) = key(1:n)
  end do
!
!  Update L.
!
  l = l - k

  return
end
subroutine r8row_reverse ( m, n, a )

!****************************************************************************80
!
!! R8ROW_REVERSE reverses the order of the rows of an R8ROW.
!
!  Discussion:
!
!    To reverse the rows is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      51 52 53 54 55
!      41 42 43 44 45
!      31 32 33 34 35
!      21 22 23 24 25
!      11 12 13 14 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) j
  real ( kind = 8 ) t(n)

  ihi = m / 2

  do i = 1, ihi
    t(1:n) = a(i,1:n)
    a(i,1:n) = a(m+1-i,1:n)
    a(m+1-i,1:n) = t(1:n)
  end do

  return
end
subroutine r8row_sort_heap_a ( m, n, a )

!*****************************************************************************80
!
!! R8ROW_SORT_HEAP_A ascending heapsorts an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
!    each of length N.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( m <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8row_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8row_sort_heap_index_a ( m, n, a, indx )

!*****************************************************************************80
!
!! R8ROW_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
!    each of length N.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    A(I1,*) < A(I1,*) if the first nonzero entry of A(I1,*)-A(I2,*)
!    is negative.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:M),1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) INDX(M), the sort index.  The I-th element
!    of the sorted array is row INDX(I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(m)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) row(n)

  if ( n < 1 ) then
    return
  end if

  do i = 1, m
    indx(i) = i
  end do

  if ( m == 1 ) then
    return
  end if

  l = ( m / 2 ) + 1
  ir = m

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      row(1:n) = a(indxt,1:n)

    else

      indxt = indx(ir)
      row(1:n) = a(indxt,1:n)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        call r8row_compare ( m, n, a, indx(j), indx(j+1), isgn )

        if ( isgn < 0 ) then
          j = j + 1
        end if

      end if

      call r8vec_compare ( n, row, a(indx(j),1:n), isgn )

      if ( isgn < 0 ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8row_sort_quick_a ( m, n, a )

!*****************************************************************************80
!
!! R8ROW_SORT_QUICK_A ascending quick sorts an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
!    each of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A,
!    and the length of a row.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 30
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) m_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n <= 0 ) then
    return
  end if

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    write ( *, '(a,i8)' ) '  M = ', m
    stop
  end if

  if ( m == 1 ) then
    return
  end if

  level = 1
  rsave(level) = m + 1
  base = 1
  m_segment = m

  do
!
!  Partition the segment.
!
    call r8row_part_quick_a ( m_segment, n, a(base:base+m_segment-1,1:n), &
      l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8ROW_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      m_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < m_segment ) then

      m_segment = m_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        m_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < m_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8row_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8ROW_SORTED_UNIQUE_COUNT counts unique elements in an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    M rows of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      unique_num = unique_num + 1
      i1 = i2
    end if

  end do

  return
end
subroutine r8row_sum ( m, n, a, rowsum )

!*****************************************************************************80
!
!! R8ROW_SUM returns the sums of the rows of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Output, real ( kind = 8 ) ROWSUM(M), the sum of the entries of
!    each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rowsum(m)

  do i = 1, m
    rowsum(i) = sum ( a(i,1:n) )
  end do

  return
end
subroutine r8row_swap ( m, n, a, i1, i2 )

!*****************************************************************************80
!
!! R8ROW_SWAP swaps two rows of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I1, I2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) row(n)

  if ( i1 < 1 .or. m < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I1 is out of range.'
    write ( *, '(a,i8)' ) '  I1 = ', i1
    stop
  end if

  if ( i2 < 1 .or. m < i2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I2 is out of range.'
    write ( *, '(a,i8)' ) '  I2 = ', i2
    stop
  end if

  if ( i1 == i2 ) then
    return
  end if

  row(1:n) = a(i1,1:n)
  a(i1,1:n) = a(i2,1:n)
  a(i2,1:n) = row(1:n)

  return
end
subroutine r8row_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Output, real ( kind = 8 ) X(M*N), a vector containing the M rows of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m*n)

  j = 1
  do i = 1, m
    x(j:j+n-1) = a(i,1:n)
    j = j + n
  end do

  return
end
subroutine r8row_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! R8ROW_VARIANCE returns the variances of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of R8 values, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array whose variances are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variances of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(m)

  do i = 1, m

    mean = sum ( a(i,1:n) ) / real ( n, kind = 8 )

    variance(i) = 0.0D+00
    do j = 1, n
      variance(i) = variance(i) + ( a(i,j) - mean )**2
    end do

    if ( 1 < n ) then
      variance(i) = variance(i) / real ( n - 1, kind = 8 )
    else
      variance(i) = 0.0D+00
    end if

  end do

  return
end
subroutine r8slmat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8SLMAT_PRINT prints a strict lower triangular R8MAT.
!
!  Example:
!
!    M = 5, N = 5
!    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)
!
!    21
!    31 32
!    41 42 43
!    51 52 53 54
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(*), the M by N matrix.  Only the strict
!    lower triangular elements are stored, in column major order.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(10)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) size
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  jmax = min ( n, m - 1 )

  if ( m-1 <= n ) then
    size = ( m * ( m - 1 ) ) / 2
  else if ( n < m-1 ) then
    size = ( n * ( n - 1 ) ) / 2 + ( m - n - 1 ) * n
  end if

  if ( all ( a(1:size) == aint ( a(1:size) ) ) ) then

    nn = 10

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a8,10i8)' ) '     Col', ( j, j = jlo, jhi )
      write ( *, '(a8)' )      '     Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,10i8)' ) i, int ( a(indx(1:jhi+1-jlo)) )
      end do
    end do

  else if ( maxval ( abs ( a(1:size) ) ) < 1000000.0D+00 ) then

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a10,5(i8,6x))' ) '       Col', ( j, j = jlo, jhi )
      write ( *, '(a10)' )          '       Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,5f14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  else

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a10,5(i8,6x))' ) '       Col', ( j, j = jlo, jhi )
      write ( *, '(a10)' ) '       Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,5g14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  end if

  return
end
subroutine r8vec_01_to_ab ( n, a, amax, amin )

!*****************************************************************************80
!
!! R8VEC_01_TO_AB shifts and rescales an R8VEC to lie within given bounds.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    On input, A contains the original data, which is presumed to lie
!    between 0 and 1.  However, it is not necessary that this be so.
!
!    On output, A has been shifted and rescaled so that all entries which
!    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
!    be mapped in a corresponding way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be rescaled.
!
!    Input, real ( kind = 8 ) AMAX, AMIN, the maximum and minimum values
!    allowed for A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amax2
  real ( kind = 8 ) amax3
  real ( kind = 8 ) amin
  real ( kind = 8 ) amin2
  real ( kind = 8 ) amin3

  if ( amax == amin ) then
    a(1:n) = amin
    return
  end if

  amax2 = max ( amax, amin )
  amin2 = min ( amax, amin )

  amin3 = minval ( a(1:n) )
  amax3 = maxval ( a(1:n) )

  if ( amax3 /= amin3 ) then

    a(1:n) = ( ( amax3 - a(1:n)         ) * amin2   &
             + (         a(1:n) - amin3 ) * amax2 ) &
             / ( amax3          - amin3 )

  else

    a(1:n) = 0.5D+00 * ( amax2 + amin2 )

  end if

  return
end
subroutine r8vec_ab_to_01 ( n, a )

!*****************************************************************************80
!
!! R8VEC_AB_TO_01 shifts and rescales an R8VEC to lie within [0,1].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    On input, A contains the original data.  On output, A has been shifted
!    and scaled so that all entries lie between 0 and 1.
!
!  Formula:
!
!    A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input/output, real ( kind = 8 ) A(N), the data to be rescaled.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amin

  amax = maxval ( a(1:n) )
  amin = minval ( a(1:n) )

  if ( amin == amax ) then
    a(1:n) = 0.5D+00
  else
    a(1:n) = ( a(1:n) - amin ) / ( amax - amin )
  end if

  return
end
subroutine r8vec_ab_to_cd ( n, a, bmin, bmax, b )

!*****************************************************************************80
!
!! R8VEC_AB_TO_CD shifts and rescales an R8VEC from one interval to another.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The mininum entry of A is mapped to BMIN, the maximum entry
!    to BMAX, and values in between are mapped linearly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(N), the data to be remapped.
!
!    Input, real ( kind = 8 ) BMIN, BMAX, the values to which min(A) and max(A)
!    are to be assigned.
!
!    Output, real ( kind = 8 ) B(N), the remapped data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amin
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bmax
  real ( kind = 8 ) bmin

  if ( bmax == bmin ) then
    b(1:n) = bmin
    return
  end if

  amin = minval ( a(1:n) )
  amax = maxval ( a(1:n) )

  if ( amax == amin ) then
    b(1:n) = 0.5D+00 * ( bmax + bmin )
    return
  end if

  b(1:n) = ( ( amax - a(1:n)        ) * bmin   &
         + (          a(1:n) - amin ) * bmax ) &
           / ( amax          - amin )

  return
end
function r8vec_all_nonpositive ( n, a )

!*****************************************************************************80
!
!! R8VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, double ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_ALL_NONPOSITIVE is TRUE if all entries
!    of A are less than or equal to zero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_all_nonpositive

  r8vec_all_nonpositive = all ( a(1:n) <= 0.0D+00 )

  return
end
subroutine r8vec_amax ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_AMAX returns the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMAX, the value of the entry
!    of largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax

  amax = maxval ( abs ( a(1:n) ) )

  return
end
subroutine r8vec_amax_index ( n, a, amax_index )

!*****************************************************************************80
!
!! R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMAX_INDEX, the index of the entry of
!    largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  integer ( kind = 4 ) amax_index
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    amax_index = -1

  else

    amax_index = 1
    amax = abs ( a(1) )

    do i = 2, n
      if ( amax < abs ( a(i) ) ) then
        amax_index = i
        amax = abs ( a(i) )
      end if
    end do

  end if

  return
end
subroutine r8vec_amin ( n, a, amin )

!*****************************************************************************80
!
!! R8VEC_AMIN returns the minimum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 )A(N), the array.
!
!    Output, real ( kind = 8 ) AMIN, the value of the entry
!    of smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amin

  amin = minval ( abs ( a(1:n) ) )

  return
end
subroutine r8vec_amin_index ( n, a, amin_index )

!*****************************************************************************80
!
!! R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN_INDEX, the index of the entry of
!    smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amin
  integer ( kind = 4 ) amin_index
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    amin_index = 0

  else

    amin_index = 1
    amin = abs ( a(1) )

    do i = 2, n
      if ( abs ( a(i) ) < amin ) then
        amin_index = i
        amin = abs ( a(i) )
      end if
    end do

  end if

  return
end
function r8vec_any_negative ( n, a )

!*****************************************************************************80
!
!! R8VEC_ANY_NEGATIVE: ( any A < 0 ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_ANY_NEGATIVE is TRUE if any entry is negative.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_any_negative

  r8vec_any_negative = any ( a(1:n) < 0.0D+00 )

  return
end
function r8vec_any_nonzero ( n, a )

!*****************************************************************************80
!
!! R8VEC_ANY_NONZERO: ( any A nonzero ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_ANY_NONZERO is TRUE if any entry is nonzero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_any_nonzero

  r8vec_any_nonzero = any ( a(1:n) /= 0.0D+00 )

  return
end
subroutine r8vec_any_normal ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_ANY_NORMAL returns some normal vector to V1.
!
!  Discussion:
!
!    If DIM_NUM < 2, then no normal vector can be returned.
!
!    If V1 is the zero vector, then any unit vector will do.
!
!    No doubt, there are better, more robust algorithms.  But I will take
!    just about ANY reasonable unit vector that is normal to V1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) V1(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) V2(DIM_NUM), a vector that is
!    normal to V2, and has unit Euclidean length.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) vj
  real ( kind = 8 ) vk

  if ( dim_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
    write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
    stop
  end if

  if ( r8vec_norm ( dim_num, v1 ) == 0.0D+00 ) then
    v2(1) = 1.0D+00
    v2(2:dim_num) = 0.0D+00
    return
  end if
!
!  Seek the largest entry in V1, VJ = V1(J), and the
!  second largest, VK = V1(K).
!
!  Since V1 does not have zero norm, we are guaranteed that
!  VJ, at least, is not zero.
!
  j = - 1
  vj = 0.0D+00

  k = - 1
  vk = 0.0D+00

  do i = 1, dim_num

    if ( abs ( vk ) < abs ( v1(i) ) .or. k < 1 ) then

      if ( abs ( vj ) < abs ( v1(i) ) .or. j < 1 ) then
        k = j
        vk = vj
        j = i
        vj = v1(i)
      else
        k = i
        vk = v1(i)
      end if

    end if

  end do
!
!  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
!  will just about do the trick.
!
  v2(1:dim_num) = 0.0D+00

  v2(j) = - vk / sqrt ( vk * vk + vj * vj )
  v2(k) =   vj / sqrt ( vk * vk + vj * vj )

  return
end
function r8vec_ascends ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For example, if:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!
!    then
!
!      R8VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS, is TRUE if the
!    entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends
  real ( kind = 8 ) x(n)

  do i = 1, n - 1
    if ( x(i+1) < x(i) ) then
      r8vec_ascends = .false.
      return
    end if
  end do

  r8vec_ascends = .true.

  return
end
function r8vec_ascends_strictly ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Notice the effect of entry number 6 in the following results:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
!
!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
!    entries of X strictly ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends_strictly
  real ( kind = 8 ) x(n)

  do i = 1, n - 1
    if ( x(i+1) <= x(i) ) then
      r8vec_ascends_strictly = .false.
      return
    end if
  end do

  r8vec_ascends_strictly = .true.

  return
end
subroutine r8vec_bin ( n, x, bin_num, bin_min, bin_max, bin, bin_limit )

!*****************************************************************************80
!
!! R8VEC_BIN computes bins based on a given R8VEC.
!
!  Discussion:
!
!    The user specifies minimum and maximum bin values, BIN_MIN and
!    BIN_MAX, and the number of bins, BIN_NUM.  This determines a
!    "bin width":
!
!      H = ( BIN_MAX - BIN_MIN ) / BIN_NUM
!
!    so that bin I will count all entries X(J) such that
!
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!
!    The array X does NOT have to be sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of X.
!
!    Input, real ( kind = 8 ) X(N), an (unsorted) array to be binned.
!
!    Input, integer ( kind = 4 ) BIN_NUM, the number of bins.  Two extra bins,
!    #0 and #BIN_NUM+1, count extreme values.
!
!    Input, real ( kind = 8 ) BIN_MIN, BIN_MAX, define the range and size
!    of the bins.  BIN_MIN and BIN_MAX must be distinct.
!    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
!    this, but proper results will be computed if BIN_MIN > BIN_MAX.
!
!    Output, integer ( kind = 4 ) BIN(0:BIN_NUM+1).
!    BIN(0) counts entries of X less than BIN_MIN.
!    BIN(BIN_NUM+1) counts entries greater than or equal to BIN_MAX.
!    For 1 <= I <= BIN_NUM, BIN(I) counts the entries X(J) such that
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!    where H is the bin spacing.
!
!    Output, real ( kind = 8 ) BIN_LIMIT(0:BIN_NUM), the "limits" of the bins.
!    BIN(I) counts the number of entries X(J) such that
!      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) bin_num

  integer ( kind = 4 ) bin(0:bin_num+1)
  real ( kind = 8 ) bin_limit(0:bin_num)
  real ( kind = 8 ) bin_max
  real ( kind = 8 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  if ( bin_max == bin_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BIN - Fatal error!'
    write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
    stop
  end if

  bin(0:bin_num+1) = 0

  do i = 1, n

    t = ( x(i) - bin_min ) / ( bin_max - bin_min )

    if ( t < 0.0D+00 ) then
      j = 0
    else if ( 1.0D+00 <= t ) then
      j = bin_num + 1
    else
      j = 1 + int ( real ( bin_num, kind = 8 ) * t )
    end if

    bin(j) = bin(j) + 1

  end do
!
!  Compute the bin limits.
!
  do i = 0, bin_num
    bin_limit(i) = (   real ( bin_num - i, kind = 8 ) * bin_min   &
                     + real (           i, kind = 8 ) * bin_max ) &
                     / real ( bin_num,     kind = 8 )
  end do

  return
end
subroutine r8vec_blend ( n, t1, x1, t2, x2, x )

!*****************************************************************************80
!
!! R8VEC_BLEND performs weighted interpolation of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The formula used is:
!
!      x(i) = t * x1(i) + (1-t) * x2(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in each  vector.
!
!    Input, real ( kind = 8 ) T1, the weight factor for vector 1.
!
!    Input, real ( kind = 8 ) X1(N), the first vector.
!
!    Input, real ( kind = 8 ) T2, the weight factor for vector 2.
!
!    Input, real ( kind = 8 ) X2(N), the second vector.
!
!    Output, real ( kind = 8 ) X(N), the interpolated or extrapolated value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1(n)
  real ( kind = 8 ) x2(n)

  x(1:n) = t1 * x1(1:n) + t2 * x2(1:n)

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    This is an inefficient implementation that uses linear search.
!
!    An R8VEC is a vector of R8's.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!    It is always true that RIGHT = LEFT+1.
!
!    If XVAL < X(1), then LEFT = 1, RIGHT = 2, and
!      XVAL   < X(1) < X(2);
!    If X(1) <= XVAL < X(N), then
!      X(LEFT) <= XVAL < X(RIGHT);
!    If X(N) <= XVAL, then LEFT = N-1, RIGHT = N, and
!      X(LEFT) <= X(RIGHT) <= XVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine r8vec_bracket2 ( n, x, xval, start, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET2 searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    containing the given value.
!
!    R8VEC_BRACKET2 is a variation on R8VEC_BRACKET.  It seeks to reduce
!    the search time by allowing the user to suggest an interval that
!    probably contains the value.  The routine will look in that interval
!    and the intervals to the immediate left and right.  If this does
!    not locate the point, a binary search will be carried out on
!    appropriate subportion of the sorted array.
!
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and X(LEFT) <= XVAL <= X(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
!      Value is greater than all data values:
!    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
!      Value is equal to a data value:
!    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed by entries of X.
!
!    Input, integer ( kind = 4 ) START, between 1 and N, specifies that XVAL
!    is likely to be in the interval:
!      [ X(START), X(START+1) ]
!    or, if not in that interval, then either
!      [ X(START+1), X(START+2) ]
!    or
!      [ X(START-1), X(START) ].
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) right
  integer ( kind = 4 ) start
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET2 - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( start < 1 .or. n < start ) then
    start = ( n + 1 ) / 2
  end if
!
!  XVAL = X(START)?
!
  if ( x(start) == xval ) then

    left = start
    right = start
    return
!
!  X(START) < XVAL?
!
  else if ( x(start) < xval ) then
!
!  X(START) = X(N) < XVAL < oo?
!
    if ( n < start + 1 ) then

      left = start
      right = -1
      return
!
!  XVAL = X(START+1)?
!
    else if ( xval == x(start+1) ) then

      left = start + 1
      right = start + 1
      return
!
!  X(START) < XVAL < X(START+1)?
!
    else if ( xval < x(start+1) ) then

      left = start
      right = start + 1
      return
!
!  X(START+1) = X(N) < XVAL < oo?
!
    else if ( n < start + 2 ) then

      left = start + 1
      right = -1
      return
!
!  XVAL = X(START+2)?
!
    else if ( xval == x(start+2) ) then

      left = start + 2
      right = start + 2
      return
!
!  X(START+1) < XVAL < X(START+2)?
!
    else if ( xval < x(start+2) ) then

      left = start + 1
      right = start + 2
      return
!
!  Binary search for XVAL in [ X(START+2), X(N) ],
!  where XVAL is guaranteed to be greater than X(START+2).
!
    else

      low = start + 2
      high = n
      call r8vec_bracket ( high + 1 - low, x(low), xval, left, right )
      left = left + low - 1
      right = right + low - 1

    end if
!
!  -oo < XVAL < X(START) = X(1).
!
  else if ( start == 1 ) then

    left = -1
    right = start
    return
!
!  XVAL = X(START-1)?
!
  else if ( xval == x(start-1) ) then

    left = start - 1
    right = start - 1
    return
!
!  X(START-1) < XVAL < X(START)?
!
  else if ( x(start-1) <= xval ) then

    left = start - 1
    right = start
    return
!
!  Binary search for XVAL in [ X(1), X(START-1) ],
!  where XVAL is guaranteed to be less than X(START-1).
!
  else

    low = 1
    high = start - 1
    call r8vec_bracket ( high + 1 - low, x(1), xval, left, right )

  end if

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) T(N), an array that has been sorted
!    into ascending order.
!
!    Input, real ( kind = 8 ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer ( kind = 4 ) LEFT.
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
subroutine r8vec_bracket4 ( nt, t, ns, s, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET4 finds the nearest interval to each of a vector of values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(NT).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, length of the input array.
!
!    Input, real ( kind = 8 ) T(NT), an array that has been sorted
!    into ascending order.
!
!    Input, integer ( kind = 4 ) NS, the number of points to be bracketed.
!
!    Input, real ( kind = 8 ) S(NS), values to be bracketed by entries of T.
!
!    Output, integer ( kind = 4 ) LEFT(NS).
!    LEFT(I) is set so that the interval [ T(LEFT(I)), T(LEFT(I)+1) ]
!    is the closest to S(I); it either contains S(I), or else S(I)
!    lies outside the interval [ T(1), T(NT) ].
!
  implicit none

  integer ( kind = 4 ) ns
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) high
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left(ns)
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real ( kind = 8 ) s(ns)
  real ( kind = 8 ) t(nt)
!
!  Check the input data.
!
  if ( nt < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET4 - Fatal error!'
    write ( *, '(a)' ) '  NT must be at least 2.'
    stop
  end if

  do i = 1, ns

    left(i) = ( nt + 1 ) / 2
!
!  CASE 1: S < T(LEFT):
!  Search for S in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
    if ( s(i) < t(left(i)) ) then

      if ( left(i) == 1 ) then
        cycle
      else if ( left(i) == 2 ) then
        left(i) = 1
        cycle
      else if ( t(left(i)-1) <= s(i) ) then
        left(i) = left(i) - 1
        cycle
      else if ( s(i) <= t(2) ) then
        left(i) = 1
        cycle
      end if
!
!  ...Binary search for S in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
      low = 2
      high = left(i) - 2

      do

        if ( low == high ) then
          left(i) = low
          exit
        end if

        mid = ( low + high + 1 ) / 2

        if ( t(mid) <= s(i) ) then
          low = mid
        else
          high = mid - 1
        end if

      end do
!
!  CASE2: T(LEFT+1) < S:
!  Search for S in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
    else if ( t(left(i)+1) < s(i) ) then

      if ( left(i) == nt - 1 ) then
        cycle
      else if ( left(i) == nt - 2 ) then
        left(i) = left(i) + 1
        cycle
      else if ( s(i) <= t(left(i)+2) ) then
        left(i) = left(i) + 1
        cycle
      else if ( t(nt-1) <= s(i) ) then
        left(i) = nt - 1
        cycle
      end if
!
!  ...Binary search for S in [T(I), T(I+1)] for intervals I = LEFT+2 to NT-2.
!
      low = left(i) + 2
      high = nt - 2

      do

        if ( low == high ) then
          left(i) = low
          exit
        end if

        mid = ( low + high + 1 ) / 2

        if ( t(mid) <= s(i) ) then
          low = mid
        else
          high = mid - 1
        end if

      end do
!
!  CASE3: T(LEFT) <= S <= T(LEFT+1):
!  S is in [T(LEFT), T(LEFT+1)].
!
    else

    end if

  end do

  return
end
function r8vec_bracket5 ( nd, xd, xi )

!*****************************************************************************80
!
!! R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
!
!  Discussion:
!
!    We assume XD is sorted.
!
!    If XI is contained in the interval [XD(1),XD(N)], then the returned 
!    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
!
!    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
!
!    This code implements a version of binary search which is perhaps more
!    understandable than the usual ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) XD(N), the sorted data.
!
!    Input, real ( kind = 8 ) XD, the query value.
!
!    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) b
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi

  if ( xi < xd(1) .or. xd(nd) < xi ) then

    b = -1

  else

    l = 1
    r = nd

    do while ( l + 1 < r )
      m = ( l + r ) / 2
      if ( xi < xd(m) ) then
        r = m
      else
        l = m
      end if
    end do

    b = l

  end if

  r8vec_bracket5 = b

  return
end
subroutine r8vec_bracket6 ( nd, xd, ni, xi, b )

!*****************************************************************************80
!
!! R8VEC_BRACKET6 brackets data between successive entries of a sorted R8VEC.
!
!  Discussion:
!
!    We assume XD is sorted.
!
!    If XI(I) is contained in the interval [XD(1),XD(N)], then the value of
!    B(I) indicates that XI(I) is contained in [ XD(B(I)), XD(B(I)+1) ].
!
!    If XI(I) is not contained in the interval [XD(1),XD(N)], then B(I) = -1.
!
!    This code implements a version of binary search which is perhaps more
!    understandable than the usual ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) XD(N), the sorted data.
!
!    Input, integer ( kind = 4 ) NI, the number of inquiry values.
!
!    Input, real ( kind = 8 ) XD(NI), the query values.
!
!    Output, integer ( kind = 4 ) B(NI), the bracket information.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) b(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)

  do i = 1, ni

    if ( xi(i) < xd(1) .or. xd(nd) < xi(i) ) then

      b(i) = -1

    else

      l = 1
      r = nd

      do while ( l + 1 < r )
        m = ( l + r ) / 2
        if ( xi(i) < xd(m) ) then
          r = m
        else
          l = m
        end if
      end do

      b(i) = l

    end if

  end do

  return
end
subroutine r8vec_ceiling ( n, r8vec, ceilingvec )

!*****************************************************************************80
!
!! R8VEC_CEILING rounds "up" (towards +oo) entries of an R8VEC.
!
!  Example:
!
!    R8    Value
!
!   -1.1  -1.0
!   -1.0  -1.0
!   -0.9   0.0
!    0.0   0.0
!    5.0   5.0
!    5.1   6.0
!    5.9   6.0
!    6.0   6.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the vector.
!
!    Output, real ( kind = 8 ) CEILINGVEC(N), the rounded values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ceilingvec(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec(n)
  real ( kind = 8 ) value

  do i = 1, n

    value = real ( int ( r8vec(i) ), kind = 8 )

    if ( value < r8vec(i) ) then
      value = value + 1.0D+00
    end if

    ceilingvec(i) = value

  end do

  return
end
subroutine r8vec_chebyspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_CHEBYSPACE creates a vector of Chebyshev spaced values in [A,B].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the interval.
!
!    Output, real ( kind = 8 ) X(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( n - i, kind = 8 ) * pi &
            / real ( n - 1, kind = 8 )

      c = cos ( theta )

      if ( mod ( n, 2 ) == 1 ) then
        if ( 2 * i - 1 == n ) then
          c = 0.0D+00
        end if
      end if

      x(i) = ( ( 1.0D+00 - c ) * a  &
             + ( 1.0D+00 + c ) * b ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_cheby1space ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_CHEBY1SPACE creates Type 1 Chebyshev spaced values in [A,B].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( 2 * ( n - i ) + 1, kind = 8 ) * pi &
        / real ( 2 * n, kind = 8 )

      c = cos ( theta )

      if ( mod ( n, 2 ) == 1 ) then
        if ( 2 * i - 1 == n ) then
          c = 0.0D+00
        end if
      end if

      x(i) = ( ( 1.0D+00 - c ) * a   &
             + ( 1.0D+00 + c ) * b ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_cheby2space ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_CHEBY2SPACE creates Type 2 Chebyshev spaced values in [A,B].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( n - i, kind = 8 ) * pi / real ( n - 1, kind = 8 )

      c = cos ( theta )

      if ( mod ( n, 2 ) == 1 ) then
        if ( 2 * i - 1 == n ) then
          c = 0.0D+00
        end if
      end if

      x(i) = ( ( 1.0D+00 - c ) * a   &
             + ( 1.0D+00 + c ) * b ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_circular_variance ( n, x, circular_variance )

!*****************************************************************************80
!
!! R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) CIRCULAR VARIANCE, the circular variance
!    of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) circular_variance
  real ( kind = 8 ) mean
  real ( kind = 8 ) x(n)

  call r8vec_mean ( n, x, mean )

  circular_variance = &
      ( sum ( cos ( x(1:n) - mean ) ) )**2 &
    + ( sum ( sin ( x(1:n) - mean ) ) )**2

  circular_variance = sqrt ( circular_variance ) / real ( n, kind = 8 )

  circular_variance = 1.0D+00 - circular_variance

  return
end
subroutine r8vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! R8VEC_COMPARE compares two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8vec_convolution ( m, x, n, y, z )

!*****************************************************************************80
!
!! R8VEC_CONVOLUTION returns the convolution of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The I-th entry of the convolution can be formed by summing the products 
!    that lie along the I-th diagonal of the following table:
!
!    Y3 | 3   4   5   6   7
!    Y2 | 2   3   4   5   6
!    Y1 | 1   2   3   4   5
!       +------------------
!        X1  X2  X3  X4  X5
!
!    which will result in:
!
!    Z = ( X1 * Y1,
!          X1 * Y2 + X2 * Y1,
!          X1 * Y3 + X2 * Y2 + X3 * Y1,
!                    X2 * Y3 + X3 * Y2 + X4 * Y1,
!                              X3 * Y3 + X4 * Y2 + X5 * Y1,
!                                        X4 * Y3 + X5 * Y2,
!                                                  X5 * Y3 )
!            
!  Example:
!
!    Input:
!
!      X = (/ 1, 2, 3, 4 /)
!      Y = (/ -1, 5, 3 /)
!
!    Output:
!
!      Z = (/ -1, 3, 10, 17, 29, 12 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of X.
!
!    Input, real ( kind = 8 ) X(M), the first vector to be convolved.
!
!    Input, integer ( kind = 4 ) N, the dimension of Y.
!
!    Input, real ( kind = 8 ) Y(N), the second vector to be convolved.
!
!    Output, real ( kind = 8 ) Z(M+N-1), the convolution of X and Y.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(m+n-1)

  z(1:m+n-1) = 0.0D+00

  do j = 1, n
    z(j:j+m-1) = z(j:j+m-1) + x(1:m) * y(j)
  end do

  return
end
subroutine r8vec_convolution_circ ( n, x, y, z )

!*****************************************************************************80
!
!! R8VEC_CONVOLUTION_CIRC: discrete circular convolution of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The formula used is:
!
!      z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
!
!    Here, if the index of Y becomes nonpositive, it is "wrapped around"
!    by having N added to it.
!
!    The circular convolution is equivalent to multiplication of Y by a
!    circulant matrix formed from the vector X.
!
!  Example:
!
!    Input:
!
!      X = (/ 1, 2, 3, 4 /)
!      Y = (/ 1, 2, 4, 8 /)
!
!    Output:
!
!      Circulant form:
!
!      Z = ( 1 4 3 2 )   ( 1 )
!          ( 2 1 4 3 )   ( 2 )
!          ( 3 2 1 4 ) * ( 4 )
!          ( 4 3 2 1 )   ( 8 )
!
!      The formula:
!
!      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
!             1*2 + 2*1 + 3*8 + 4*4,
!             1*4 + 2*2 + 3*1 + 4*8,
!             1*8 + 2*4 + 3*2 + 4*1 /)
!
!      Result:
!
!      Z = (/ 37, 44, 43, 26 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vectors to be convolved.
!
!    Output, real ( kind = 8 ) Z(N), the circular convolution of X and Y.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  do m = 1, n
    z(m) = dot_product ( x(1:m), y(m:1:-1) ) &
         + dot_product ( x(m+1:n), y(n:m+1:-1) )
  end do

  return
end
subroutine r8vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_COPY copies an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), the vector to be copied.
!
!    Output, real ( kind = 8 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine r8vec_correlation ( n, x, y, correlation )

!*****************************************************************************80
!
!! R8VEC_CORRELATION returns the correlation of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If X and Y are two nonzero vectors of length N, then
!
!      correlation = (x/||x||)' (y/||y||)
!
!    It is the cosine of the angle between the two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vectors to be convolved.
!
!    Output, real ( kind = 8 ) CORRELATION, the correlation of X and Y.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) correlation
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_norm
  real ( kind = 8 ) xy_dot
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y_norm

  x_norm = r8vec_norm ( n, x )
  y_norm = r8vec_norm ( n, y )
  xy_dot = dot_product ( x(1:n), y(1:n) )

  if ( x_norm == 0.0D+00 .or. y_norm == 0.0D+00 ) then
    correlation = 0.0D+00
  else
    correlation = xy_dot / x_norm / y_norm
  end if

  return
end
function r8vec_covar ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_COVAR computes the covariance of two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2013
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), Y(N), the two vectors.
!
!    Input, integer ( kind = 4 ) N, the dimension of the two vectors.
!
!    Output, real ( kind = 8 ) R4VEC_COVAR, the covariance of the vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec_covar
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_average
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y_average

  x_average = sum ( x(1:n) ) / real ( n, kind = 8 )
  y_average = sum ( y(1:n) ) / real ( n, kind = 8 )
 
  value = 0.0D+00
  do i = 1, n
    value = value + ( x(i) - x_average ) * ( y(i) - y_average )
  end do

  r8vec_covar = value / real ( n - 1, kind = 8 )

  return
end
function r8vec_cross_product_2d ( v1, v2 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of vectors in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors V1 and V2 should be considered
!    to lie in a 3D space, both having Z coordinate zero.  The cross
!    product value V3 then represents the standard cross product vector
!    (0,0,V3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(2), V2(2), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_CROSS_PRODUCT_2D, the cross product.
!
  implicit none

  real ( kind = 8 ) r8vec_cross_product_2d
  real ( kind = 8 ) v1(2)
  real ( kind = 8 ) v2(2)

  r8vec_cross_product_2d = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
function r8vec_cross_product_affine_2d ( v0, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors V1 and V2 should be considered
!    to lie in a 3D space, both having Z coordinate zero.  The cross
!    product value V3 then represents the standard cross product vector
!    (0,0,V3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V0(2), the base vector.
!
!    Input, real ( kind = 8 ) V1(2), V2(2), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_CROSS_PRODUCT_AFFINE_2D,
!    the cross product (V1-V0) x (V2-V0).
!
  implicit none

  real ( kind = 8 ) r8vec_cross_product_affine_2d
  real ( kind = 8 ) v0(2)
  real ( kind = 8 ) v1(2)
  real ( kind = 8 ) v2(2)

  r8vec_cross_product_affine_2d = &
      ( v1(1) - v0(1) ) * ( v2(2) - v0(2) ) &
    - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

  return
end
subroutine r8vec_cross_product_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
subroutine r8vec_cross_product_affine_3d ( v0, v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!    Here, we use V0 as the base of an affine system so we compute
!    the cross product of (V1-V0) and (V2-V0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V0(3), the base vector.
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector
!    ( V1-V0) x (V2-V0).
!
  implicit none

  real ( kind = 8 ) v0(3)
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  v3(1) = ( v1(2) - v0(2) ) * ( v2(3) - v0(3) ) &
        - ( v2(2) - v0(2) ) * ( v1(3) - v0(3) )

  v3(2) = ( v1(3) - v0(3) ) * ( v2(1) - v0(1) ) &
        - ( v2(3) - v0(3) ) * ( v1(1) - v0(1) )

  v3(3) = ( v1(1) - v0(1) ) * ( v2(2) - v0(2) ) &
        - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

  return
end
subroutine r8vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! R8VEC_CUM computes the cumulutive sums of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Input:
!
!      A = (/ 1.0, 2.0, 3.0, 4.0 /)
!
!    Output:
!
!      A_CUM = (/ 1.0, 3.0, 6.0, 10.0 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be summed.
!
!    Output, real ( kind = 8 ) A_CUM(1:N), the cumulative sums.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_cum(n)
  integer ( kind = 4 ) i

  a_cum(1) = a(1)

  do i = 2, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine r8vec_cum0 ( n, a, a_cum )

!*****************************************************************************80
!
!! R8VEC_CUM0 computes the cumulutive sums of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Input:
!
!      A = (/ 1.0, 2.0, 3.0, 4.0 /)
!
!    Output:
!
!      A_CUM = (/ 0.0, 1.0, 3.0, 6.0, 10.0 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be summed.
!
!    Output, real ( kind = 8 ) A_CUM(0:N), the cumulative sums.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0.0D+00

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine r8vec_dif ( n, h, cof )

!*****************************************************************************80
!
!! R8VEC_DIF computes coefficients for estimating the N-th derivative.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine computes the N+1 coefficients for a centered finite difference
!    estimate of the N-th derivative of a function.
!
!    The estimate has the form
!
!      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
!
!    To understand the computation of the coefficients, it is enough
!    to realize that the first difference approximation is
!
!      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
!
!    and that the second difference approximation can be regarded as
!    the first difference approximation repeated:
!
!      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
!         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
!
!    and so on for higher order differences.
!
!    Thus, the next thing to consider is the integer coefficients of
!    the sampled values of F, which are clearly the Pascal coefficients,
!    but with an alternating negative sign.  In particular, if we
!    consider row I of Pascal's triangle to have entries j = 0 through I,
!    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
!    and P(0,0) = 1.
!
!       1
!      -1  1
!       1 -2   1
!      -1  3  -3   1
!       1 -4   6  -4   1
!      -1  5 -10  10  -5  1
!       1 -6  15 -20  15 -6 1
!
!    Next, note that the denominator of the approximation for the
!    N-th derivative will be (2*DX)^N.
!
!    And finally, consider the location of the N+1 sampling
!    points for F:
!
!      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
!
!    Thus, a formula for evaluating FDIF(N,X) is
!
!      fdif = 0.0
!      do i = 0, n
!        xi = x + (2*i-n) * h
!        fdif = fdif + cof(i) * f(xi)
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the derivative to be
!    approximated.  N must be 0 or greater.
!
!    Input, real ( kind = 8 ) H, the half spacing between points.
!    H must be positive.
!
!    Output, real ( kind = 8 ) COF(0:N), the coefficients needed to approximate
!    the N-th derivative of a function F.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cof(0:n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_DIF - Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative order N = ', n
    write ( *, '(a)' ) '  but N must be at least 0.'
    stop
  end if

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_DIF - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The half sampling spacing is H = ', h
    write ( *, '(a)' ) '  but H must be positive.'
    stop
  end if

  do i = 0, n

    cof(i) = 1.0D+00

    do j = i - 1, 1, -1
      cof(j) = -cof(j) + cof(j-1)
    end do

    if ( 0 < i ) then
      cof(0) = -cof(0)
    end if

  end do

  cof(0:n) = cof(0:n) / ( 2.0D+00 * h )**n

  return
end
function r8vec_diff_dot_product ( n, u1, v1, u2, v2 )

!*****************************************************************************80
!
!! R8VEC_DIFF_DOT_PRODUCT: dot product of a pair of R8VEC differences.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) U1(N), V1(N), defines the vector U1-V1.
!
!    Input, real ( kind = 8 ) U2(N), V2(N), defines the vector U2-V2.
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_DOT_PRODUCT, the dot product 
!    of (U1-V1)*(U2-V2).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec_diff_dot_product
  real ( kind = 8 ) u1(n)
  real ( kind = 8 ) u2(n)
  real ( kind = 8 ) v1(n)
  real ( kind = 8 ) v2(n)
  real ( kind = 8 ) value

  value = 0.0D+00
  do i = 1, n
    value = value + ( u1(i) - v1(i) ) * ( u2(i) - v2(i) )
  end do

  r8vec_diff_dot_product = value

  return
end
function r8vec_diff_norm ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM, the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm

  r8vec_diff_norm = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

  return
end
function r8vec_diff_norm_l1 ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM_L1 returns the L1 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L1 norm is defined as:
!
!      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM_L1, the L1 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm_l1

  r8vec_diff_norm_l1 = sum ( abs ( a(1:n) - b(1:n) ) )

  return
end
function r8vec_diff_norm_l2 ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM_L2 returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM_L2, the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm_l2

  r8vec_diff_norm_l2 = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

  return
end
function r8vec_diff_norm_li ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L-oo norm is defined as:
!
!      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM_LI, the L-oo norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm_li

  r8vec_diff_norm_li = maxval ( abs ( a(1:n) - b(1:n) ) )

  return
end
function r8vec_diff_norm_squared ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A(I) - B(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM_SQUARED, the square of 
!    the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm_squared

  r8vec_diff_norm_squared = sum ( ( a(1:n) - b(1:n) )**2 )

  return
end
subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of 
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
function r8vec_distance ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) V1(DIM_NUM), V2(DIM_NUM), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DISTANCE, the Euclidean distance
!    between the vectors.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) r8vec_distance
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)

  r8vec_distance = sqrt ( sum ( ( v1(1:dim_num) - v2(1:dim_num) )**2 ) )

  return
end
function r8vec_distinct ( n, a )

!*****************************************************************************80
!
!! R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be checked.
!
!    Output, logical R8VEC_DISTINCT is TRUE if the elements of A are distinct.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct

  r8vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        return
      end if
    end do
  end do

  r8vec_distinct = .true.

  return
end
function r8vec_dot_product ( n, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In FORTRAN90, the system routine DOT_PRODUCT should be called
!    directly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) V1(N), V2(N), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DOT_PRODUCT, the dot product.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_dot_product
  real ( kind = 8 ) v1(n)
  real ( kind = 8 ) v2(n)

  r8vec_dot_product = dot_product ( v1(1:n), v2(1:n) )

  return
end
function r8vec_dot_product_affine ( n, v0, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product V1-V0 * V2-V0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) V0(N), the base vector.
!
!    Input, real ( kind = 8 ) V1(N), V2(N), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DOT_PRODUCT_AFFINE, the dot product.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_dot_product_affine
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)
  real ( kind = 8 ) v2(n)

  r8vec_dot_product_affine = dot_product ( &
    v1(1:n) - v0(1:n),  &
    v2(1:n) - v0(1:n) )

  return
end
function r8vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_EQ is true if two R8VECs are equal.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ, is TRUE if every pair of elements A1(I)
!    and A2(I) are equal, and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = 8 ) * xlo   &
           + real (     ival - 1, kind = 8 ) * xhi ) &
           / real ( n        - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_even2 ( maxval, nfill, nold, xold, nval, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN2 linearly interpolates new numbers into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The number of values created between two old values can vary from
!    one pair of values to the next.
!
!    The interpolated values are evenly spaced.
!
!    This routine is a generalization of R8VEC_EVEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXVAL, the size of the XVAL array, as declared
!    by the user.  MAXVAL must be large enough to hold the NVAL values computed
!    by this routine.  In other words, MAXVAL must be at least equal to
!    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).
!
!    Input, integer ( kind = 4 ) NFILL(NOLD-1), the number of values
!    to be interpolated between XOLD(I) and XOLD(I+1).
!    NFILL(I) does not count the endpoints.  Thus, if
!    NFILL(I) is 1, there will be one new point generated
!    between XOLD(I) and XOLD(I+1).
!    NFILL(I) must be nonnegative.
!
!    Input, integer ( kind = 4 ) NOLD, the number of values XOLD,
!    between which extra values are to be interpolated.
!
!    Input, real ( kind = 8 ) XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, integer ( kind = 4 ) NVAL, the number of values computed
!    in the XVAL array.
!    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)
!
!    Output, real ( kind = 8 ) XVAL(MAXVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as the interpolated
!    values, making a total of NVAL values.
!
  implicit none

  integer ( kind = 4 ) maxval
  integer ( kind = 4 ) nold

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nadd
  integer ( kind = 4 ) nfill(nold-1)
  integer ( kind = 4 ) nval
  real ( kind = 8 ) xold(nold)
  real ( kind = 8 ) xval(maxval)

  nval = 1

  do i = 1, nold - 1

    if ( nfill(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_EVEN2 - Fatal error!'
      write ( *, '(a,i8)' ) '  NFILL(I) is negative for I = ', i
      write ( *, '(a,i8)' ) '  NFILL(I) = ', nfill(i)
      stop
    end if

    if ( maxval < nval + nfill(i) + 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_EVEN2 - Fatal error!'
      write ( *, '(a)' ) '  MAXVAL is not large enough.  '
      write ( *, '(a,i8)' ) '  MAXVAL = ', maxval
      write ( *, '(a)' ) '  which is exceeded by storage requirements'
      write ( *, '(a,i8)' ) '  for interpolating in interval ', i
      stop
    end if

    nadd = nfill(i) + 2

    do j = 1, nadd
      xval(nval+j-1) = ( real ( nadd - j,     kind = 8 ) * xold(i)   &
                       + real (        j - 1, kind = 8 ) * xold(i+1) ) &
                       / real ( nadd     - 1, kind = 8 )
    end do

    nval = nval + nfill(i) + 1

  end do

  return
end
subroutine r8vec_even2_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN2_SELECT returns the I-th of N evenly spaced midpoint values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This function returns the I-th of N evenly spaced midpoints of N
!    equal subintervals of [XLO,XHI].
!
!    XVAL = ( ( 2 * N - 2 * IVAL + 1 ) * XLO 
!           + (         2 * IVAL - 1 ) * XHI ) 
!           / ( 2 * N                )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced midpoints
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  xval = ( real ( 2 * n - 2 * ival + 1, kind = 8 ) * xlo   &
         + real (         2 * ival - 1, kind = 8 ) * xhi ) &
         / real ( 2 * n, kind = 8 )

  return
end
subroutine r8vec_even3 ( nold, nval, xold, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN3 evenly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This routine accepts a short vector of numbers, and returns a longer
!    vector of numbers, created by interpolating new values between
!    the given values.
!
!    Between any two original values, new values are evenly interpolated.
!
!    Over the whole vector, the new numbers are interpolated in
!    such a way as to try to minimize the largest distance interval size.
!
!    The algorithm employed is not "perfect".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NOLD, the number of values XOLD, between
!    which extra values are to be interpolated.
!
!    Input, integer ( kind = 4 ) NVAL, the number of values to be computed
!    in the XVAL array.  NVAL should be at least NOLD.
!
!    Input, real ( kind = 8 ) XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, real ( kind = 8 ) XVAL(NVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as interpolated
!    values, making a total of NVAL values.
!
  implicit none

  integer ( kind = 4 ) nval
  integer ( kind = 4 ) nold

  real ( kind = 8 ) density
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nmaybe
  integer ( kind = 4 ) npts
  integer ( kind = 4 ) ntemp
  integer ( kind = 4 ) ntot
  real ( kind = 8 ) xlen
  real ( kind = 8 ) xleni
  real ( kind = 8 ) xlentot
  real ( kind = 8 ) xold(nold)
  real ( kind = 8 ) xval(nval)

  xlen = 0.0D+00
  do i = 1, nold - 1
    xlen = xlen + abs ( xold(i+1) - xold(i) )
  end do

  ntemp = nval - nold

  density = real ( ntemp, kind = 8 ) / xlen

  ival = 1
  ntot = 0
  xlentot = 0.0D+00

  do i = 1, nold - 1

    xleni = abs ( xold(i+1) - xold(i) )
    npts = int ( density * xleni )
    ntot = ntot + npts
!
!  Determine if we have enough left-over density that it should
!  be changed into a point.  A better algorithm would agonize
!  more over where that point should go.
!
    xlentot = xlentot + xleni
    nmaybe = nint ( xlentot * density )

    if ( ntot < nmaybe ) then
      npts = npts + nmaybe - ntot
      ntot = nmaybe
    end if

    do j = 1, npts + 2
      xval(ival+j-1) = ( real ( npts+2 - j,     kind = 8 ) * xold(i)   &
                       + real (          j - 1, kind = 8 ) * xold(i+1) ) &
                       / real ( npts+2     - 1, kind = 8 )
    end do

    ival = ival + npts + 1

  end do

  return
end
subroutine r8vec_expand_linear ( n, x, fat, xfat )

!*****************************************************************************80
!
!! R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This routine copies the old data, and inserts NFAT new values
!    between each pair of old data values.  This would be one way to
!    determine places to evenly sample a curve, given the (unevenly
!    spaced) points at which it was interpolated.
!
!  Example:
!
!    N = 3
!    NFAT = 2
!
!    X(1:N)        = (/ 0.0,           6.0,             7.0 /)
!    XFAT(1:2*3+1) = (/ 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of input data values.
!
!    Input, real ( kind = 8 ) X(N), the original data.
!
!    Input, integer ( kind = 4 ) FAT, the number of data values to interpolate
!    between each pair of original data values.
!
!    Output, real ( kind = 8 ) XFAT((N-1)*(FAT+1)+1), the "fattened" data.
!
  implicit none

  integer ( kind = 4 ) fat
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xfat((n-1)*(fat+1)+1)

  k = 0

  do i = 1, n - 1

    k = k + 1
    xfat(k) = x(i)

    do j = 1, fat
      k = k + 1
      xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(i)     &
                + real (       j,     kind = 8 ) * x(i+1) ) &
                / real ( fat     + 1, kind = 8 )
    end do

  end do

  k = k + 1
  xfat(k) = x(n)

  return
end
subroutine r8vec_expand_linear2 ( n, x, before, fat, after, xfat )

!*****************************************************************************80
!
!! R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This routine starts with a vector of data.
!
!    The intent is to "fatten" the data, that is, to insert more points
!    between successive values of the original data.
!
!    There will also be extra points placed BEFORE the first original
!    value and AFTER that last original value.
!
!    The "fattened" data is equally spaced between the original points.
!
!    The BEFORE data uses the spacing of the first original interval,
!    and the AFTER data uses the spacing of the last original interval.
!
!  Example:
!
!    N = 3
!    BEFORE = 3
!    FAT = 2
!    AFTER = 1
!
!    X    = (/                   0.0,           6.0,             7.0       /)
!    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
!            3 "BEFORE's"        Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of input data values.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the original data.
!
!    Input, integer ( kind = 4 ) BEFORE, the number of "before" values.
!
!    Input, integer ( kind = 4 ) FAT, the number of data values to interpolate
!    between each pair of original data values.
!
!    Input, integer ( kind = 4 ) AFTER, the number of "after" values.
!
!    Output, real ( kind = 8 ) XFAT(BEFORE+(N-1)*(FAT+1)+1+AFTER), the
!    "fattened" data.
!
  implicit none

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xfat(before+(n-1)*(fat+1)+1+after)

  k = 0
!
!  Points BEFORE.
!
  do j = 1 - before + fat, fat
    k = k + 1
    xfat(k) = ( real ( fat - j + 1, kind = 8 ) * ( x(1) - ( x(2) - x(1) ) ) &
              + real (       j,     kind = 8 ) *   x(1)          ) &
              / real ( fat     + 1, kind = 8 )
  end do
!
!  Original points and FAT points.
!
  do i = 1, n - 1

    k = k + 1
    xfat(k) = x(i)

    do j = 1, fat
      k = k + 1
      xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(i)     &
                + real (       j,     kind = 8 ) * x(i+1) ) &
                / real ( fat     + 1, kind = 8 )
    end do

  end do

  k = k + 1
  xfat(k) = x(n)
!
!  Points AFTER.
!
  do j = 1, after
    k = k + 1
    xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(n)     &
              + real (       j,     kind = 8 ) &
              * ( x(n) + ( x(n) - x(n-1) ) ) ) &
              / real ( fat     + 1, kind = 8 )
  end do

  return
end
subroutine r8vec_first_index ( n, a, tol, first_index )

!*****************************************************************************80
!
!! R8VEC_FIRST_INDEX indexes the first occurrence of values in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
!    the first occurrence of the value A(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) FIRST_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) first_index(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol

  first_index(1:n) = -1

  do i = 1, n

    if ( first_index(i) == -1 ) then

      first_index(i) = i

      do j = i + 1, n
        if ( abs ( a(i) - a(j) ) <= tol ) then
          first_index(j) = i
        end if
      end do

    end if

  end do

  return
end
subroutine r8vec_floor ( n, r8vec, floorvec )

!*****************************************************************************80
!
!! R8VEC_FLOOR rounds "down" (towards -oo) entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    R8    Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the values to be rounded down.
!
!    Output, integer ( kind = 4 ) FLOORVEC(N), the rounded value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) floorvec(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) value

  do i = 1, n

    value = int ( r8vec(i) )

    if ( r8vec(i) < real ( value, kind = 8 ) ) then
      value = value - 1
    end if

    floorvec(i) = value

  end do

  return
end
subroutine r8vec_frac ( n, a, k, frac )

!*****************************************************************************80
!
!! R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Hoare's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, A is the array to search.
!    On output, the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.  Other values
!    of K search for the entry which is K-th in size.  K must be at
!    least 1, and no greater than N.
!
!    Output, real ( kind = 8 ) FRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) frac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  real ( kind = 8 ) temp
  real ( kind = 8 ) x

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      frac = a(k)
      exit
    end if

    x = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then
        if ( j < k ) then
          left = i
        end if
        if ( k < i ) then
          iryt = j
        end if
        exit
      end if
!
!  Find I so that X <= A(I).
!
      do while ( a(i) < x )
        i = i + 1
      end do
!
!  Find J so that A(J) <= X.
!
      do while ( x < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then

        temp = a(i)
        a(i) = a(j)
        a(j) = temp

        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine r8vec_fraction ( n, x, fraction )

!*****************************************************************************80
!
!! R8VEC_FRACTION returns the fraction parts of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If we regard a real number as
!
!      R8 = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R8    R8_FRACTION
!
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(N), the arguments.
!
!    Output, real ( kind = 8 ) FRACTION(N), the fraction parts.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fraction(n)
  real ( kind = 8 ) x(n)

  fraction(1:n) = abs ( x(1:n) ) - real ( int ( abs ( x(1:n) ) ), kind = 8 )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>                              A1(1) > A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  logical r8vec_gt

  r8vec_gt = .false.

  do i = 1, n

    if ( a2(i) < a1(i) ) then
      r8vec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r8vec_gt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_heap_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_HEAP_A reorders an R8VEC into an ascending heap.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real ( kind = 8 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine r8vec_heap_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_HEAP_D reorders an R8VEC into an descending heap.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real ( kind = 8 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine r8vec_heap_d_extract ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_HEAP_D_EXTRACT: extract maximum from a heap descending sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In other words, the routine finds the maximum value in the
!    heap, returns that value to the user, deletes that value from
!    the heap, and restores the heap to its proper form.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, real ( kind = 8 ) A(N), the heap.
!
!    Output, real ( kind = 8 ) VALUE, the item of maximum value, which has
!    been removed from the heap.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) n
  real ( kind = 8 ) value

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the maximum value.
!
  value = a(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last value down.
!
  a(1) = a(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call r8vec_sort_heap_d ( n, a )

  return
end
subroutine r8vec_heap_d_insert ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_HEAP_D_INSERT inserts a value into a heap descending sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, real ( kind = 8 ) A(N), the heap.
!
!    Input, real ( kind = 8 ) VALUE, the value to be inserted.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent
  real ( kind = 8 ) value

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( value <= a(parent) ) then
      exit
    end if

    a(i) = a(parent)
    i = parent

  end do

  a(i) = value

  return
end
subroutine r8vec_heap_d_max ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_HEAP_D_MAX returns the maximum value in a heap descending sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input, real ( kind = 8 ) A(N), the heap.
!
!    Output, real ( kind = 8 ) VALUE, the maximum value in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) value

  value = a(1)

  return
end
subroutine r8vec_histogram ( n, a, a_lo, a_hi, histo_num, histo_gram )

!*****************************************************************************80
!
!! R8VEC_HISTOGRAM histograms an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Values between A_LO and A_HI will be histogrammed into the bins
!    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
!    and values greater than A_HI are counted in bin HISTO_NUM+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array to examine.
!
!    Input, real ( kind = 8 ) A_LO, A_HI, the lowest and highest
!    values to be histogrammed.  These values will also define the bins.
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the number of bins to use.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM+1), contains the
!    number of entries of A in each bin.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_hi
  real ( kind = 8 ) a_lo
  real ( kind = 8 ) delta
  integer ( kind = 4 ) histo_gram(0:histo_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  histo_gram(0:histo_num+1) = 0

  delta = ( a_hi - a_lo ) / real ( 2 * histo_num, kind = 8 )

  do i = 1, n

    if ( a(i) < a_lo ) then

      histo_gram(0) = histo_gram(0) + 1

    else if ( a(i) <= a_hi ) then

      j = nint ( &
        ( ( a_hi -           delta - a(i)        ) &
        * real ( 1,         kind = 8 )   &
        + (      -           delta + a(i) - a_lo ) &
        * real ( histo_num, kind = 8 ) ) &
        / ( a_hi - 2.0D+00 * delta        - a_lo ) )

      histo_gram(j) = histo_gram(j) + 1

    else if ( a_hi < a(i) ) then

      histo_gram(histo_num+1) = histo_gram(histo_num+1) + 1

    end if

  end do

  return
end
subroutine r8vec_house_column ( n, a, k, v )

!*****************************************************************************80
!
!! R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine returns a vector V that defines a Householder
!    premultiplier matrix H(V) that zeros out the subdiagonal entries of
!    column K of the matrix A.
!
!       H(V) = I - 2 * v * v'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N), column K of the matrix A.
!
!    Input, integer ( kind = 4 ) K, the column of the matrix to be modified.
!
!    Output, real ( kind = 8 ) V(N), a vector of unit L2 norm which defines an
!    orthogonal Householder premultiplier matrix H with the property
!    that the K-th column of H*A is zero below the diagonal.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) s
  real ( kind = 8 ) v(n)

  v(1:n) = 0.0D+00

  if ( k < 1 .or. n <= k ) then
    return
  end if

  s = sqrt ( dot_product ( a(k:n), a(k:n) ) )

  if ( s == 0.0D+00 ) then
    return
  end if

  v(k) = a(k) + sign ( s, a(k) )
  v(k+1:n) = a(k+1:n)

  v(k:n) = v(k:n) / sqrt ( dot_product ( v(k:n), v(k:n) ) )

  return
end
function r8vec_i4vec_dot_product ( n, r8vec, i4vec )

!*****************************************************************************80
!
!! R8VEC_I4VEC_DOT_PRODUCT finds the dot product of an R8VEC and an I4VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) R8VEC(N), the first vector.
!
!    Input, integer ( kind = 4 ) I4VEC(N), the second vector.
!
!    Output, real ( kind = 8 ) R8VEC_I4VEC_DOT_PRODUCT, the dot product.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i4vec(n)
  real ( kind = 8 ) r8vec(n)
  real ( kind = 8 ) r8vec_i4vec_dot_product

  r8vec_i4vec_dot_product = dot_product ( r8vec(1:n), &
                                   real ( i4vec(1:n), kind = 8 ) )

  return
end
function r8vec_in_01 ( n, a )

!*****************************************************************************80
!
!! R8VEC_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_in_01

  if ( any ( a(1:n) < 0.0D+00 .or. 1.0D+00 < a(1:n) ) ) then
    r8vec_in_01 = .false.
  else
    r8vec_in_01 = .true.
  end if

  return
end
function r8vec_in_ab ( n, x, a, b )

!*****************************************************************************80
!
!! R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in X.
!
!    Input, real ( kind = 8 ) X(N), the vector.
!
!    Input, real ( kind = 8 ) A, B, the limits of the range.
!
!    Output, logical R8VEC_IN_AB, is TRUE if every entry of A is
!    between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical r8vec_in_ab
  real ( kind = 8 ) x(n)

  if ( any ( x(1:n) < a .or. b < x(1:n) ) ) then
    r8vec_in_ab = .false.
  else
    r8vec_in_ab = .true.
  end if

  return
end
subroutine r8vec_index_delete_all ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_ALL deletes a value from an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Note that the value of N is adjusted because of the deletions!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) equal1
  integer ( kind = 4 ) equal2
  integer ( kind = 4 ) get
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) put
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xval

  if ( n < 1 ) then
    n = 0
    return
  end if

  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    return
  end if

  equal1 = equal

  do

    if ( equal1 <= 1 ) then
      exit
    end if

    if ( x(indx(equal1-1)) /= xval ) then
      exit
    end if

    equal1 = equal1 - 1

  end do

  equal2 = equal

  do

    if ( n <= equal2 ) then
      exit
    end if

    if ( x(indx(equal2+1)) /= xval ) then
      exit
    end if

    equal2 = equal2 + 1

  end do
!
!  Discard certain X values.
!
  put = 0

  do get = 1, n

    if ( x(get) /= xval ) then
      put = put + 1
      x(put) = x(get)
    end if

  end do

  x(put+1:n) = 0.0D+00
!
!  Adjust the INDX values.
!
  do equal = equal1, equal2
    do i = 1, n
      if ( indx(equal) < indx(i) ) then
        indx(i) = indx(i) - 1
      end if
    end do
  end do
!
!  Discard certain INDX values.
!
  indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
  indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
  n = put

  return
end
subroutine r8vec_index_delete_dupes ( n, x, indx, n2, x2, indx2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The output quantities N2, X2, and INDX2 are computed from the
!    input quantities by sorting, and eliminating duplicates.
!
!    The output arrays should be dimensioned of size N, unless the user
!    knows in advance what the value of N2 will be.
!
!    The output arrays may be identified with the input arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input list.
!
!    Input, real ( kind = 8 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique entries in X.
!
!    Output, real ( kind = 8 ) X2(N2), a copy of the list which has
!    been sorted, and made unique.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the new list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n)
  real ( kind = 8 ) x3(n)

  i = 0
  n3 = 0

  do

    i = i + 1

    if ( n < i ) then
      exit
    end if

    if ( 1 < i ) then
      if ( x(indx(i)) == x3(n3) ) then
        cycle
      end if
    end if

    n3 = n3 + 1
    x3(n3) = x(indx(i))

  end do
!
!  Copy data into output arrays.
!
  n2 = n3
  x2(1:n2) = x3(1:n3)
  call i4vec_indicator ( n2, indx2 )

  return
end
subroutine r8vec_index_delete_one ( n, x, indx, xval, n2, x2, indx2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_ONE deletes one copy of a value from indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!    Note that the value of N is adjusted because of the deletions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) N2, the size of the current list.
!
!    Output, real ( kind = 8 ) X2(N2), the list.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n)
  real ( kind = 8 ) xval

  if ( n < 1 ) then
    n2 = 0
    return
  end if

  n2 = n
  indx2(1:n2) = indx(1:n2)
  x2(1:n2) = x(1:n2)

  call r8vec_index_search ( n2, x2, indx2, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx2(equal)
    x2(j:n2-1) = x2(j+1:n2)
    indx2(equal:n2-1) = indx2(equal+1:n2)
    do i = 1, n2-1
      if ( j < indx2(i) ) then
        indx2(i) = indx2(i) - 1
      end if
    end do
    n2 = n2 - 1
  end if

  return
end
subroutine r8vec_index_insert ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_INSERT inserts a value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine r8vec_index_insert_unique ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If the value does not occur in the list, it is included in the list,
!    and N, X and INDX are updated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine r8vec_index_order ( n, x, indx )

!*****************************************************************************80
!
!! R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine r8vec_index_search ( n, x, indx, xval, less, equal, more )

!*****************************************************************************80
!
!! R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xhi < xval ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xmid < xval ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8vec_index_sort_unique ( n, x, indx, n2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_SORT_UNIQUE creates a sorted unique index for an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.  On output, X contains only
!    unique elements.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique elements in X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) n2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  n2 = 0

  do i = 1, n
    call r8vec_index_insert_unique ( n2, y, indx, x(i) )
  end do

  x(1:n2) = y(1:n2)

  x(n2+1:n) = 0.0D+00
  indx(n2+1:n) = 0

  return
end
subroutine r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, i_lo, i_hi )

!*****************************************************************************80
!
!! R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!
!    Input, real ( kind = 8 ) R(N), the index sorted vector.
!
!    Input, integer ( kind = 4 ) INDX(N), the vector used to sort R.
!    The vector R(INDX(*)) is sorted.
!
!    Input, real ( kind = 8 ) R_LO, R_HI, the limits of the range.
!
!    Output, integer ( kind = 4 ) I_LO, I_HI, the range of indices
!    so that I_LO <= I <= I_HI => R_LO <= R(INDX(I)) <= R_HI.  If no
!    values in R lie in the range, then I_HI < I_LO will be returned.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
!
!  Cases we can handle immediately.
!
  if ( r(indx(n)) < r_lo ) then
    i_lo = n + 1
    i_hi = n
    return
  end if

  if ( r_hi < r(indx(1)) ) then
    i_lo = 1
    i_hi = 0
    return
  end if
!
!  Are there are least two intervals?
!
  if ( n == 1 ) then
    if ( r_lo <= r(indx(1)) .and. r(indx(1)) <= r_hi ) then
      i_lo = 1
      i_hi = 1
    else
      i_lo = 0
      i_hi = -1
    end if
    return
  end if
!
!  Bracket R_LO.
!
  if ( r_lo <= r(indx(1)) ) then

    i_lo = 1

  else
!
!  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
!  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
!  Does R_LO lie here, or below or above?
!
    j1 = 1
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_lo < r(indx(i1)) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(indx(i2)) < r_lo ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_lo = i1
        exit
      end if

    end do

  end if
!
!  Bracket R_HI.
!
  if ( r(indx(n)) <= r_hi ) then

    i_hi = n

  else

    j1 = i_lo
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_hi < r(indx(i1)) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(indx(i2)) < r_hi ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_hi = i2
        exit
      end if

    end do

  end if
!
!  We expect to have computed the largest I_LO and smallest I_HI such that
!    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
!  but what we want is actually
!    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
!  which we can usually get simply by incrementing I_LO and decrementing I_HI.
!
  if ( r(indx(i_lo)) < r_lo ) then
    i_lo = i_lo + 1
    if ( n < i_lo ) then
      i_hi = i_lo - 1
    end if
  end if

  if ( r_hi < r(indx(i_hi)) ) then
    i_hi = i_hi - 1
    if ( i_hi < 1 ) then
      i_lo = i_hi + 1
    end if
  end if

  return
end
subroutine r8vec_indexed_heap_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D creates a descending heap from an indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
!    each referencing an entry of the data vector.
!
!    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
!    we have:
!      A(INDX(2*J))   <= A(INDX(J))
!    and
!      A(INDX(2*J+1)) <= A(INDX(J))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the index array.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index array.
!    Each entry of INDX must be a valid index for the array A.
!    On output, the indices have been reordered into a descending heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = indx(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(indx(m)) < a(indx(m+1)) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(indx(m)) <= a(key) ) then
        exit
      end if

      indx(ifree) = indx(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    indx(ifree) = key

  end do

  return
end
subroutine r8vec_indexed_heap_d_extract ( n, a, indx, indx_extract )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
!    each referencing an entry of the data vector.
!
!    The routine finds the maximum value in the heap, returns that value to the
!    user, deletes that value from the heap, and restores the heap to its
!    proper form.
!
!    Note that the argument N must be a variable, which will be decremented
!    before return, and that INDX will hold one less value on output than it
!    held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the
!    index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_EXTRACT, the index in A of the item of
!    maximum value, which has now been removed from the heap.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_extract
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the index of the maximum value.
!
  indx_extract = indx(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last index down.
!
  indx(1) = indx(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call r8vec_indexed_heap_d ( n, a, indx )

  return
end
subroutine r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
!    each referencing an entry of the data vector.
!
!    Note that the argument N must be a variable, and will be incremented before
!    return, and that INDX must be able to hold one more entry on output than
!    it held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the
!    index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Input, integer ( kind = 4 ) INDX_INSERT, the index in A of the value
!    to be inserted into the heap.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_insert
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( a(indx_insert) <= a(indx(parent)) ) then
      exit
    end if

    indx(i) = indx(parent)
    i = parent

  end do

  indx(i) = indx_insert

  return
end
subroutine r8vec_indexed_heap_d_max ( n, a, indx, indx_max )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
!    each referencing an entry of the data vector.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_MAX, the index in A of the maximum value
!    in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx_max

  indx_max = indx(1)

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_insert ( n, a, pos, value )

!*****************************************************************************80
!
!! R8VEC_INSERT inserts a value into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the array on input.
!
!    Input/output, real ( kind = 8 ) A(N+1), the array.  On input, A is
!    assumed to contain only N entries, while on output, A actually
!    contains N+1 entries.
!
!    Input, integer ( kind = 4 ) POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, real ( kind = 8 ) VALUE, the value to be inserted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pos
  real ( kind = 8 ) value

  if ( pos < 1 .or. n + 1 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_INSERT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n + 1, pos + 1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
function r8vec_insignificant ( n, r, s )

!*****************************************************************************80
!
!! R8VEC_INSIGNIFICANT determines if an R8VEC is insignificant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) R(N), the vector to be compared against.
!
!    Input, real ( kind = 8 ) S(N), the vector to be compared.
!
!    Output, logical R8VEC_INSIGNIFICANT, is TRUE if S is insignificant
!    compared to R.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n)
  logical r8vec_insignificant
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  logical value

  value = .true.

  do i = 1, n

    t = r(i) + s(i)
    tol = epsilon ( r(i) ) * abs ( r(i) )

    if ( tol < abs ( r(i) - t ) ) then 
      value = .false.
      exit
    end if

  end do
  
  r8vec_insignificant = value

  return
end
function r8vec_is_int ( n, a )

!*****************************************************************************80
!
!! R8VEC_IS_INT is TRUE if the entries of an R8VEC are integers.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IS_INT, is TRUE if every entry of A is
!    integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_is_int

  r8vec_is_int = all ( a(1:n) == aint ( a(1:n) ) )

  return
end
function r8vec_is_nonnegative ( n, a )

!*****************************************************************************80
!
!! R8VEC_IS_NONNEGATIVE is TRUE if all the entries of an R8VEC are nonnegative.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IS_NONNEGATIVE, the value of the condition.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_is_nonnegative

  r8vec_is_nonnegative = all ( 0.0D+00 <= a(1:n) )

  return
end
function r8vec_is_zero ( n, a )

!*****************************************************************************80
!
!! R8VEC_IS_ZERO is TRUE if all the entries of an R8VEC are zero.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IS_ZERO, the value of the condition.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_is_zero

  r8vec_is_zero = all ( a(1:n) == 0.0D+00 )

  return
end
subroutine r8vec_legendre ( n, x_first, x_last, x )

!*****************************************************************************80
!
!! R8VEC_LEGENDRE creates a vector of Legendre-spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X_FIRST, X_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of Legendre-spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  call legendre_zeros ( n, x )

  x(1:n) = ( ( 1.0D+00 - x(1:n) ) * x_first  &
           + ( 1.0D+00 + x(1:n) ) * x_last ) &
           /   2.0D+00

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * a   &
             + real (     i - 1, kind = 8 ) * b ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_linspace2 ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE2 creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 2, 4, 6, 8, 10.
!
!    In other words, the interval is divided into N+1 even subintervals,
!    and the endpoints of internal intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = ( real ( n  - i + 1, kind = 8 ) * a &
           + real (      i,     kind = 8 ) * b ) &
           / real ( n      + 1, kind = 8 )
  end do

  return
end
function r8vec_lt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_LT == ( A1 < A2 ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>                              A1(1) < A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_lt
  integer ( kind = 4 ) i

  r8vec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      r8vec_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r8vec_lt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_mask_print ( n, a, mask_num, mask, title )

!*****************************************************************************80
!
!! R8VEC_MASK_PRINT prints a masked R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of masked elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the indices of the vector
!    to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mask(mask_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Masked vector printout:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, mask_num
    write ( *, '(2x,i8,a,1x,i8,2x,g14.6)' ) i, ':', mask(i), a(mask(i))
  end do

  return
end
subroutine r8vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine r8vec_max_abs_index ( n, a, max_index )

!*****************************************************************************80
!
!! R8VEC_MAX_ABS_INDEX: index of the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MAX_INDEX, the index of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_index

  if ( n <= 0 ) then

    max_index = -1

  else

    max_index = 1

    do i = 2, n
      if ( abs ( a(max_index) ) < abs ( a(i) ) ) then
        max_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_max_index ( n, a, max_index )

!*****************************************************************************80
!
!! R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MAX_INDEX, the index of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_index

  if ( n <= 0 ) then

    max_index = -1

  else

    max_index = 1

    do i = 2, n
      if ( a(max_index) < a(i) ) then
        max_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine r8vec_median ( n, a, median )

!*****************************************************************************80
!
!! R8VEC_MEDIAN returns the median of an unsorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, real ( kind = 8 ) MEDIAN, the value of the median of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) median

  k = ( n + 1 ) / 2

  call r8vec_frac ( n, a, k, median )

  return
end
subroutine r8vec_midspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_MIDSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This function divides the interval [a,b] into n subintervals, and then
!    returns the midpoints of those subintervals.
!
!  Example:
!
!    N = 5, A = 10, B = 20
!    X = [ 11, 13, 15, 17, 19 ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = 8 ) * a &
           + real (         2 * i - 1, kind = 8 ) * b ) &
           / real ( 2 * n,             kind = 8 )
  end do

  return
end
subroutine r8vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R8VEC_MIN returns the minimum value of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine r8vec_min_index ( n, a, min_index )

!*****************************************************************************80
!
!! R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MIN_INDEX, the index of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) min_index

  if ( n <= 0 ) then

    min_index = -1

  else

    min_index = 1

    do i = 2, n
      if ( a(i) < a(min_index) ) then
        min_index = i
      end if
    end do

  end if

  return
end
function r8vec_min_pos ( n, a )

!*****************************************************************************80
!
!! R8VEC_MIN_POS returns the minimum positive value of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_MIN_POS, the smallest positive entry,
!    or R8_HUGE if no entry is positive.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  real ( kind = 8 ) r8vec_min_pos
  real ( kind = 8 ) value

  value = r8_huge

  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      value = min ( value, a(i) )
    end if
  end do

  r8vec_min_pos = value

  return
end
subroutine r8vec_mirror_next ( n, a, done )

!*****************************************************************************80
!
!! R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In normal use, the user would set every element of A to be positive.
!    The routine will take the input value of A, and output a copy in
!    which the signs of one or more entries have been changed.  Repeatedly
!    calling the routine with the output from the previous call will generate
!    every distinct "variation" of A; that is, all possible sign variations.
!
!    When the output variable DONE is TRUE (or equal to 1), then the
!    output value of A_NEW is the last in the series.
!
!    Note that A may have some zero values.  The routine will essentially
!    ignore such entries; more exactly, it will not stupidly assume that -0
!    is a proper "variation" of 0!
!
!    Also, it is possible to call this routine with the signs of A set
!    in any way you like.  The routine will operate properly, but it
!    will nonethess terminate when it reaches the value of A in which
!    every nonzero entry has negative sign.
!
!    More efficient algorithms using the Gray code seem to require internal
!    memory in the routine, which is not one of MATLAB's strong points,
!    or the passing back and forth of a "memory array", or the use of
!    global variables, or unnatural demands on the user.  This form of
!    the routine is about as clean as I can make it.
!
!  Example:
!
!      Input         Output
!    ---------    --------------
!    A            A_NEW     DONE
!    ---------    --------  ----
!     1  2  3     -1  2  3  false
!    -1  2  3      1 -2  3  false
!     1 -2  3     -1 -2  3  false
!    -1 -2  3      1  2 -3  false
!     1  2 -3     -1  2 -3  false
!    -1  2 -3      1 -2 -3  false
!     1 -2 -3     -1 -2 -3  false
!    -1 -2 -3      1  2  3  true
!
!     1  0  3     -1  0  3  false
!    -1  0  3      1  0 -3  false
!     1  0 -3     -1  0 -3  false
!    -1  0 -3      1  0  3  true
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), a vector of real numbers.
!    On output, the signs of some entries have been changed.
!
!    Output, logical DONE, is TRUE if the input vector A was the last element
!    in the series (every entry was nonpositive); the output vector is reset
!    so that all entries are nonnegative, but presumably the ride is over!
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) positive
!
!  Seek the first strictly positive entry of A.
!
  positive = 0
  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      positive = i
      exit
    end if
  end do
!
!  If there is no strictly positive entry of A, there is no successor.
!
  if ( positive == 0 ) then
    a(1:n) = - a(1:n)
    done = .true.
    return
  end if
!
!  Otherwise, negate A up to the positive entry.
!
  a(1:positive) = - a(1:positive)
  done = .false.

  return
end
function r8vec_negative_strict ( n, a )

!*****************************************************************************80
!
!! R8VEC_NEGATIVE_STRICT: every element of an R8VEC is strictly negative.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N).
!
!    Output, logical R8VEC_NEGATIVE_STRICT, is TRUE every entry of the
!    vector is strictly negative.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_negative_strict

  r8vec_negative_strict = ( all ( a(1:n) < 0.0D+00 ) )

  return
end
subroutine r8vec_nint ( n, a )

!*****************************************************************************80
!
!! R8VEC_NINT rounds entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)

  a(1:n) = nint ( real ( a(1:n), kind = 8 ) )

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the vectors.
!
!    Input, real ( kind = 8 ) V0(N), the base vector.
!
!    Input, real ( kind = 8 ) V1(N), the vector whose affine norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end
function r8vec_norm_l0 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L0 returns the l0 "norm" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The l0 "norm" simply counts the number of nonzero entries in the vector.
!    It is not a true norm, but has some similarities to one.  It is useful
!    in the study of compressive sensing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, integer ( kind = 4 ) R8VEC_NORM_L0, the value of the norm.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r8vec_norm_l0
  integer ( kind = 4 ) value

  value = 0
  do i = 1, n
    if ( a(i) /= 0.0D+00 ) then
      value = value + 1
    end if
  end do

  r8vec_norm_l0 = value

  return
end
function r8vec_norm_l1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L1 returns the L1 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L1 norm is defined as:
!
!      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L1 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_l1

  r8vec_norm_l1 = sum ( abs ( a(1:n) ) )

  return
end
function r8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_l2

  r8vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r8vec_norm_li ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_LI returns the L-oo norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L-oo norm is defined as:
!
!      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L-oo norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_li

  r8vec_norm_li = maxval ( abs ( a(1:n) ) )

  return
end
function r8vec_norm_lp ( n, a, p )

!*****************************************************************************80
!
!! R8VEC_NORM_LP returns the LP norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector LP norm is defined as:
!
!      R8VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )^P )^(1/P).
!
!    Usually, the LP norms with
!      1 <= P <= oo
!    are of interest.  This routine allows
!      0 < P <= Huge ( P ).
!    If P = Huge ( P ), then the L-oo norm is returned, which is
!    simply the maximum of the absolute values of the vector components.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose LP norm is desired.
!
!    Input, real ( kind = 8 ) P, the index of the norm.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_LP, the LP norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) p
  real ( kind = 8 ) r8vec_norm_lp

  if ( p <= 0.0D+00 ) then
    r8vec_norm_lp = -1.0D+00
  else if ( p == huge ( p ) ) then
    r8vec_norm_lp = maxval ( abs ( a(1:n) ) )
  else if ( p == 1.0D+00 ) then
    r8vec_norm_lp = sum ( abs ( a(1:n) ) )
  else if ( p == 2.0D+00 ) then
    r8vec_norm_lp = sqrt ( sum ( a(1:n)**2 ) )
  else
    r8vec_norm_lp = ( sum ( ( abs ( a(1:n) ) )**p ) )**( 1.0D+00 / p )
  end if

  return
end
function r8vec_norm_squared ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_SQUARED returns the square of the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    R8VEC_NORM_SQUARED = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_SQUARED, the squared L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_squared

  r8vec_norm_squared = sum ( a(1:n)**2 )

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_normalize ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORMALIZE normalizes an R8VEC in the Euclidean norm.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The euclidean norm is also sometimes called the l2 or
!    least squares norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be normalized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) norm

  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0D+00 ) then
    a(1:n) = a(1:n) / norm
  end if

  return
end
subroutine r8vec_normalize_l1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORMALIZE_L1 normalizes an R8VEC to have unit sum.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be normalized.
!    On output, the entries of A should have unit sum.  However, if
!    the input vector has zero sum, the routine halts.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_sum

  a_sum = sum ( a(1:n) )

  if ( a_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_NORMALIZE_L1 - Fatal error!'
    write ( *, '(a)' ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
function r8vec_normsq ( n, v )

!*****************************************************************************80
!
!! R8VEC_NORMSQ returns the square of the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The square of the vector L2 norm is defined as:
!
!      R8VEC_NORMSQ = sum ( 1 <= I <= N ) V(I)^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the vector dimension.
!
!    Input, real ( kind = 8 ) V(N), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_NORMSQ, the squared L2 norm.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_normsq
  real ( kind = 8 ) v(n)

  r8vec_normsq = sum ( v(1:n)**2 )

  return
end
function r8vec_normsq_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORMSQ_AFFINE returns the affine squared norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine squared vector L2 norm is defined as:
!
!      R8VEC_NORMSQ_AFFINE(V0,V1)
!        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the vector dimension.
!
!    Input, real ( kind = 8 ) V0(N), the base vector.
!
!    Input, real ( kind = 8 ) V1(N), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_NORMSQ_AFFINE, the squared affine L2 norm.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_normsq_affine
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)

  r8vec_normsq_affine = sum ( ( v0(1:n) - v1(1:n) )**2 )

  return
end
subroutine r8vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R8VEC_ORDER_TYPE determines if R8VEC is (non)strictly ascending/descending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine r8vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine reorders the entries of A.  Using A(1) as the key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define
!    the three segments.  Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    KEY < A(I).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  real ( kind = 8 ) temp

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( key < a(l+1) ) then
      r = r - 1
      temp = a(r)
      a(r) = a(l+1)
      a(l+1) = temp
    else if ( a(l+1) == key ) then
      m = m + 1
      temp = a(m)
      a(m) = a(l+1)
      a(l+1) = temp
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine r8vec_permute ( n, p, a )

!*****************************************************************************80
!
!! R8VEC_PERMUTE permutes an R8VEC in place.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!    P(I) = J means that the I-th element of the output array should be
!    the J-th element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_temp
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine r8vec_permute_cyclic ( n, k, a )

!*****************************************************************************80
!
!! R8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For 0 <= K < N, this function cyclically permutes the input vector
!    to have the form
!
!     ( A(K+1), A(K+2), ..., A(N), A(1), ..., A(K) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) K, the increment used.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ipk
  integer ( kind = 4 ) k

  do i = 1, n
    ipk = i4_wrap ( i + k, 1, n )
    b(i) = a(ipk)
  end do

  a(1:n) = b(1:n)

  return
end
subroutine r8vec_permute_uniform ( n, a, seed )

!*****************************************************************************80
!
!! R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call perm_uniform ( n, base, seed, p )

  call r8vec_permute ( n, p, a )

  return
end
subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the vector to be polarized.
!
!    Input, real ( kind = 8 ) P(N), the polarizing direction.
!
!    Output, real ( kind = 8 ) A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_dot_p
  real ( kind = 8 ) a_normal(n)
  real ( kind = 8 ) a_parallel(n)
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
    return
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)

  return
end
function r8vec_positive_strict ( n, a )

!*****************************************************************************80
!
!! R8VEC_POSITIVE_STRICT: every element of an R8VEC is strictly positive.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N).
!
!    Output, logical R8VEC_POSITIVE_STRICT, is TRUE every entry of the
!    vector is strictly positive.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical r8vec_positive_strict

  r8vec_positive_strict = ( all ( 0.0D+00 < a(1:n) ) )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,a,1x,g14.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_print2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_PRINT2 prints out an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amin
  integer ( kind = 4 ) i
  character ( len = 11 ) iform
  logical integ
  integer ( kind = 4 ) lmax
  real ( kind = 8 ) r8_log_10
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, n

    if ( a(i) /= real ( int ( a(i) ), kind = 8 ) ) then
      integ = .false.
      exit
    end if

  end do
!
!  Find the range of the array.
!
  amax = maxval ( abs ( a(1:n) ) )
  amin = minval ( abs ( a(1:n) ) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real vectors.
!
  lmax = int ( r8_log_10 ( amax ) )

  if ( integ ) then
    write ( iform, '( ''(2x,i'', i2, '')'' )' ) lmax + 3
  else
    iform = ' '
  end if

  do i = 1, n

    if ( integ ) then
      write ( *, iform ) int ( a(i) )
    else
      write ( *, '(2x,g14.6)' ) a(i)
    end if

  end do

  return
end
function r8vec_product ( n, a )

!*****************************************************************************80
!
!! R8VEC_PRODUCT returns the product of the entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      R8VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      R8VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_product

  r8vec_product = product ( a(1:n) )

  return
end
subroutine r8vec_range ( n, x, xmin, xmax, y, ymin, ymax )

!*****************************************************************************80
!
!! R8VEC_RANGE finds the range of Y's within a restricted X range.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine is given a set of pairs of points (X,Y), and a range
!    XMIN to XMAX of valid X values.  Over this range, it seeks
!    YMIN and YMAX, the minimum and maximum values of Y for
!    valid X's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the X array.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the range of X values to check.
!
!    Input, real ( kind = 8 ) Y(N), the Y array.
!
!    Output, real ( kind = 8 ) YMIN, YMAX, the range of Y values whose
!    X value is within the X range.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  ymin =   huge ( ymin )
  ymax = - huge ( ymax )

  do i = 1, n

    if ( xmin <= x(i) .and. x(i) <= xmax ) then

      ymin = min ( ymin, y(i) )
      ymax = max ( ymax, y(i) )

    end if

  end do

  return
end
subroutine r8vec_range_2 ( n, a, amin, amax )

!*****************************************************************************80
!
!! R8VEC_RANGE_2 updates a range to include a new array.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Given a range AMIN to AMAX, and an array A, the routine will
!    decrease AMIN if necessary, or increase AMAX if necessary, so that
!    every entry of A is between AMIN and AMAX.
!
!    However, AMIN will not be increased, nor AMAX decreased.
!
!    This routine may be used to compute the maximum and minimum of a
!    collection of arrays one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Input/output, real ( kind = 8 ) AMIN, AMAX.  On input, the
!    current legal range of values for A.  On output, AMIN and AMAX
!    are either unchanged, or else "widened" so that all entries
!    of A are within the range.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) amin

  amax = max ( amax, maxval ( a(1:n) ) )
  amin = min ( amin, minval ( a(1:n) ) )

  return
end
subroutine r8vec_reverse ( n, a )

!*****************************************************************************80
!
!! R8VEC_REVERSE reverses the elements of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In FORTRAN90, calling R8VEC_REVERSE is equivalent to
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)

  a(1:n) = a(n:1:-1)

  return
end
function r8vec_rms ( n, a )

!*****************************************************************************80
!
!! R8VEC_RMS returns the RMS norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector RMS norm is defined as:
!
!      R8VEC_RMS = sqrt ( sum ( 1 <= I <= N ) A(I)^2 / N ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_RMS, the RMS norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_rms

  r8vec_rms = sqrt ( sum ( a(1:n)**2 ) / n )

  return
end
subroutine r8vec_rotate ( n, a, m )

!*****************************************************************************80
!
!! R8VEC_ROTATE "rotates" the entries of an R8VEC in place.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This routine rotates an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  real ( kind = 8 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy
      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

  end do

  return
end
function r8vec_scalar_triple_product ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.
!
!  Discussion:
!
!    STRIPLE = V1 dot ( V2 x V3 ).
!
!    STRIPLE is the volume of the parallelogram whose sides are
!    formed by V1, V2 and V3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the three vectors.
!
!    Output, real ( kind = 8 ) R8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
!    triple product.
!
  implicit none

  real ( kind = 8 ) r8vec_scalar_triple_product
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  r8vec_scalar_triple_product = &
      v1(1) * ( v2(2) * v3(3) - v2(3) * v3(2) ) &
    + v1(2) * ( v2(3) * v3(1) - v2(1) * v3(3) ) &
    + v1(3) * ( v2(1) * v3(2) - v2(2) * v3(1) )

  return
end
subroutine r8vec_search_binary_a ( n, a, aval, indx )

!*****************************************************************************80
!
!! R8VEC_SEARCH_BINARY_A searches an ascending sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be searched.  The array must
!    be sorted in ascending order.
!
!    Input, real ( kind = 8 ) AVAL, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, AVAL does not occur in the array.
!    I, A(I) = AVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = -1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == aval ) then
      indx = mid
      exit
    else if ( a(mid) < aval ) then
      low = mid + 1
    else if ( aval < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine r8vec_shift ( shift, n, x )

!*****************************************************************************80
!
!! R8VEC_SHIFT performs a shift on an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SHIFT, the amount by which each entry is to
!    be shifted.
!
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be shifted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) shift
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = x(1:n)

  x(1:n) = 0.0D+00

  ilo = max ( 1, 1 + shift )
  ihi = min ( n, n + shift )

  x(ilo:ihi) = y(ilo-shift:ihi-shift)

  return
end
subroutine r8vec_shift_circular ( shift, n, x )

!*****************************************************************************80
!
!! R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SHIFT, the amount by which each entry is to
!    be shifted.
!
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be shifted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) shift
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = x(1:n)

  do i = 1, n
    j = i4_wrap ( i - shift, 1, n )
    x(i) = y(j)
  end do

  return
end
subroutine r8vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(j) < a(i) ) then
        t    = a(i)
        a(i) = a(j)
        a(j) = t
      end if
    end do
  end do

  return
end
subroutine r8vec_sort_bubble_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(i) < a(j) ) then
        t    = a(i)
        a(i) = a(j)
        a(j) = t
      end if
    end do
  end do

  return
end
subroutine r8vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) n1
  real ( kind = 8 ) temp

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call r8vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  temp = a(1)
  a(1) = a(n)
  a(n) = temp
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call r8vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    temp = a(1)
    a(1) = a(n1)
    a(n1) = temp

  end do

  return
end
subroutine r8vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call r8vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call r8_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call r8vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call r8_swap ( a(1), a(n1) )

  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_heap_mask_a ( n, a, mask_num, mask, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An array A is given.  An array MASK of indices into A is given.
!    The routine produces a vector INDX, which is a permutation of the
!    entries of MASK, so that:
!
!      A(MASK(INDX(I)) <= A(MASK(INDX(J))
!
!    whenever
!
!      I <= J
!
!    In other words, only the elements of A that are indexed by MASK
!    are to be considered, and the only thing that happens is that
!    a rearrangment of the indices in MASK is returned that orders the
!    masked elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of mask elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the mask array.  This is
!    simply a list of indices of A.  The entries of MASK should
!    be unique, and each one should be between 1 and N.
!
!    Output, integer ( kind = 4 ) INDX(MASK_NUM), the sort index.  There are
!    MASK_NUM elements of A selected by MASK.  If we want to list those
!    elements in order, then the I-th element is A(MASK(INDX(I))).
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(mask_num)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mask(mask_num)

  if ( n < 1 ) then
    return
  end if

  if ( mask_num < 1 ) then
    return
  end if

  if ( mask_num == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( mask_num, indx )

  l = mask_num / 2 + 1
  ir = mask_num

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(mask(indxt))

    else

      indxt = indx(ir)
      aval = a(mask(indxt))
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(mask(indx(j))) < a(mask(indx(j+1))) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(mask(indx(j))) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine r8vec_sort_insert_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sorted indices.  The array
!    is sorted when listed from A(INDX(1)) through A(INDX(N)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(indx(j)) <= x ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine r8vec_sort_insert_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_INDEX_D descending index sorts an R8VEC using insertion.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sorted indices.  The array
!    is sorted when listed from A(INDX(1)) through A(INDX(N)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(indx(j)) ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine r8vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    Input:
!
!      N = 7
!      A = ( 6, 7, 3, 2, 9, 1, 8 )
!
!    Output:
!
!      A = ( 1, 2, 3, 6, 7, 8, 9 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 30
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r8vec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8vec_sort_shell_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) asave
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) maxpow

  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3**MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 3**maxpow < 2 * n + 1 )
    maxpow = maxpow + 1
  end do

  if ( 1 < maxpow ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3^IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc+k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine r8vec_sort2_a ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SORT2_A ascending sorts an R8VEC and adjusts an associated R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine sorts the elements of X, and whenever
!    an element of X is moved, the corresponding element of
!    Y is moved in the same way.  This action means that after
!    the sorting, every element of X is still paired to the
!    same Y value.
!
!    If you have more than one array associated with X, or
!    an integer array, or some other complication, you may want to
!    look at doing an "indexed sort" instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an unsorted array.
!    On output, X has been sorted.
!
!    Input/output, real ( kind = 8 ) Y(N), an array which is to be
!    shifted corresponding to the shifts made in X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call r8_swap ( x(i), x(j) )
      call r8_swap ( y(i), y(j) )

    else if ( indx < 0 ) then

      if ( x(i) <= x(j) ) then
        isgn = -1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec_sorted_merge_a ( na, a, nb, b, nc, c )

!*****************************************************************************80
!
!! R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(NA), the first sorted array.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, real ( kind = 8 ) B(NB), the second sorted array.
!
!    Output, integer ( kind = 4 ) NC, the number of elements in the output
!    array.  Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, real ( kind = 8 ) C(NC), the merged unique sorted array.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) c(na+nb)
  real ( kind = 8 ) d(na+nb)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) order

  na2 = na
  nb2 = nb

  ja = 0
  jb = 0
  nc = 0

  call r8vec_order_type ( na2, a, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORTED_MERGE_A - Fatal error!'
    write ( *, '(a)' ) '  The input array A is not ascending sorted!'
    stop
  end if

  call r8vec_order_type ( nb2, b, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORTED_MERGE_A - Fatal error!'
    write ( *, '(a)' ) '  The input array B is not ascending sorted!'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( na2 <= ja ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = b(jb)
        else if ( d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( nb2 <= jb ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = a(ja)
        else if ( d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = a(ja)
      else if ( d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = b(jb)
      else if ( d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
function r8vec_sorted_nearest ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), a sorted vector.
!
!    Input, real ( kind = 8 ) VALUE, the value whose nearest vector
!    entry is sought.
!
!    Output, integer ( kind = 4 ) R8VEC_SORTED_NEAREST, the index of the nearest
!    entry in the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) r8vec_sorted_nearest
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  real ( kind = 8 ) value

  if ( n < 1 ) then
    r8vec_sorted_nearest = -1
    return
  end if

  if ( n == 1 ) then
    r8vec_sorted_nearest = 1
    return
  end if

  if ( a(1) < a(n) ) then

    if ( value < a(1) ) then
      r8vec_sorted_nearest = 1
      return
    else if ( a(n) < value ) then
      r8vec_sorted_nearest = n
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = 1
    hi = n

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return
!
!  A descending sorted vector A.
!
  else

    if ( value < a(n) ) then
      r8vec_sorted_nearest = n
      return
    else if ( a(1) < value ) then
      r8vec_sorted_nearest = 1
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = n
    hi = 1

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return

  end if

  return
end
subroutine r8vec_sorted_range ( n, r, r_lo, r_hi, i_lo, i_hi )

!*****************************************************************************80
!
!! R8VEC_SORTED_RANGE searches a sorted vector for elements in a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!
!    Input, real ( kind = 8 ) R(N), the sorted vector.
!
!    Input, real ( kind = 8 ) R_LO, R_HI, the limits of the range.
!
!    Output, integer ( kind = 4 ) I_LO, I_HI, the range of indices
!    so that I_LO <= I <= I_HI => R_LO <= R(I) <= R_HI.  If no
!    values in R lie in the range, then I_HI < I_LO will be returned.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
!
!  Cases we can handle immediately.
!
  if ( r(n) < r_lo ) then
    i_lo = 0
    i_hi = - 1
    return
  end if

  if ( r_hi < r(1) ) then
    i_lo = 0
    i_hi = - 1
    return
  end if
!
!  Are there are least two intervals?
!
  if ( n == 1 ) then
    if ( r_lo <= r(1) .and. r(1) <= r_hi ) then
      i_lo = 1
      i_hi = 1
    else
      i_lo = 0
      i_hi = -1
    end if
    return
  end if
!
!  Bracket R_LO.
!
  if ( r_lo <= r(1) ) then

    i_lo = 1

  else
!
!  R_LO is in one of the intervals spanned by R(J1) to R(J2).
!  Examine the intermediate interval [R(I1), R(I1+1)].
!  Does R_LO lie here, or below or above?
!
    j1 = 1
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_lo < r(i1) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(i2) < r_lo ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_lo = i1
        exit
      end if

    end do

  end if
!
!  Bracket R_HI
!
  if ( r(n) <= r_hi ) then

    i_hi = n

  else

    j1 = i_lo
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_hi < r(i1) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(i2) < r_hi ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_hi = i2
        exit
      end if

    end do

  end if
!
!  We expect to have computed the largest I_LO and smallest I_HI such that
!    R(I_LO) <= R_LO <= R_HI <= R(I_HI)
!  but what we want is actually
!    R_LO <= R(I_LO) <= R(I_HI) <= R_HI
!  which we can usually get simply by incrementing I_LO and decrementing I_HI.
!
  if ( r(i_lo) < r_lo ) then
    i_lo = i_lo + 1
    if ( n < i_lo ) then
      i_hi = i_lo - 1
    end if
  end if

  if ( r_hi < r(i_hi) ) then
    i_hi = i_hi - 1
    if ( i_hi < 1 ) then
      i_lo = i_hi + 1
    end if
  end if

  return
end
subroutine r8vec_sorted_split ( n, a, split, i_lt, i_gt )

!*****************************************************************************80
!
!! R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Given a splitting value SPLIT, the routine seeks indices
!    I_LT and I_GT so that
!
!      A(I_LT) < SPLIT < A(I_GT),
!
!    and if there are intermediate index values between I_LT and
!    I_GT, then those entries of A are exactly equal to SPLIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), a sorted array.
!
!    Input, real ( kind = 8 ) SPLIT, a value to which the entries in A are
!    to be compared.
!
!    Output, integer ( kind = 4 ) I_LT:
!    0 if no entries are less than SPLIT;
!    N if all entries are less than SPLIT;
!    otherwise, the index of the last entry in A less than SPLIT.
!
!    Output, integer ( kind = 4 ) I_GT:
!    1 if all entries are greater than SPLIT;
!    N+1 if no entries are greater than SPLIT;
!    otherwise the index of the first entry in A greater than SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_gt
  integer ( kind = 4 ) i_lt
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  real ( kind = 8 ) split

  if ( n < 1 ) then
    i_lt = -1
    i_gt = -1
    return
  end if

  if ( split < a(1) ) then
    i_lt = 0
    i_gt = 1
    return
  end if

  if ( a(n) < split ) then
    i_lt = n
    i_gt = n + 1
    return
  end if

  lo = 1
  hi = n

  do

    if ( lo + 1 == hi ) then
      i_lt = lo
      exit
    end if

    mid = ( lo + hi ) / 2

    if ( split <= a(mid) ) then
      hi = mid
    else
      lo = mid
    end if

  end do

  do i = i_lt + 1, n
    if ( split < a(i) ) then
      i_gt = i
      return
    end if
  end do

  i_gt = n + 1

  return
end
subroutine r8vec_sorted_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNDEX returns unique sorted indexes for a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    Assuming X is already sorted, we examine the entries of X in order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector X could be
!    replaced by a compressed vector XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(I) = X(UNDX(I)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!    Here is an example of a vector X, the unique sort and
!    inverse unique sort vectors and the compressed unique sorted vector.
!
!      I      X      XU  Undx  Xdnu
!    ----+------+------+-----+-----+
!      1 | 11.0 |  11.0    1     1
!      2 | 11.0 |  22.0    5     1
!      3 | 11.0 |  33.0    8     1
!      4 | 11.0 |  55.0    9     1
!      5 | 22.0 |                2
!      6 | 22.0 |                2
!      7 | 22.0 |                2
!      8 | 33.0 |                3
!      9 | 55.0 |
!
!    INDX(2) = 3 means that sorted item(2) is X(3).
!    XDNI(2) = 5 means that X(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
!    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of data values.
!
!    Input, real ( kind = 8 ) X_VAL(X_NUM), the data values.
!
!    Input, integer ( kind = 4 ) X_UNIQUE_NUM, the number of unique values
!    in X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) x_unique_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(x_unique_num)
  real ( kind = 8 ) x_val(x_num)
  integer ( kind = 4 ) xdnu(x_num)
!
!  Walk through the sorted array X.
!
  i = 1

  j = 1
  undx(j) = i

  xdnu(i) = j

  do i = 2, x_num

    if ( tol < abs ( x_val(i) - x_val(undx(j)) ) ) then
      j = j + 1
      undx(j) = i
    end if

    xdnu(i) = j

  end do

  return
end
subroutine r8vec_sorted_unique ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE keeps the unique elements in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the sorted array of N elements;
!    On output, the sorted unique array of UNIQUE_NUM elements.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real ( kind = 8 ) tol

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( tol < abs ( a(i) - a(unique_num) ) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(i)
    end if

  end do

  return
end
subroutine r8vec_sorted_unique_count ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Because the array is sorted, this algorithm is O(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the sorted array to examine.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real ( kind = 8 ) tol

  if ( n < 1 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( tol < abs ( a(i-1) - a(i) ) ) then
      unique_num = unique_num + 1
    end if

  end do

  return
end
subroutine r8vec_sorted_unique_hist ( n, a, tol, maxuniq, unique_num, &
  auniq, acount )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array to examine.  The elements of A
!    should have been sorted.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Input, integer ( kind = 4 ) MAXUNIQ, the maximum number of unique elements
!    that can be handled.  If there are more than MAXUNIQ unique
!    elements in A, the excess will be ignored.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
!    Output, real ( kind = 8 ) AUNIQ(UNIQUE_NUM), the unique elements of A.
!
!    Output, integer ( kind = 4 ) ACOUNT(UNIQUE_NUM), the number of times
!    each element of AUNIQ occurs in A.
!
  implicit none

  integer ( kind = 4 ) maxuniq
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) acount(maxuniq)
  real ( kind = 8 ) auniq(maxuniq)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real ( kind = 8 ) tol
!
!  Start taking statistics.
!
  unique_num = 0

  do i = 1, n

    if ( i == 1 ) then

      unique_num = 1
      auniq(unique_num) = a(1)
      acount(unique_num) = 1

    else if ( abs ( a(i) - auniq(unique_num) ) <= tol ) then

      acount(unique_num) = acount(unique_num) + 1

    else if ( unique_num < maxuniq ) then

      unique_num = unique_num + 1
      auniq(unique_num) = a(i)
      acount(unique_num) = 1

    end if

  end do

  return
end
subroutine r8vec_split ( n, a, split, isplit )

!*****************************************************************************80
!
!! R8VEC_SPLIT "splits" an unsorted R8VEC based on a splitting value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:ISPLIT).
!
!    Input, real ( kind = 8 ) SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ( kind = 4 ) ISPLIT, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) isplit
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  real ( kind = 8 ) split
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n + 1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call r8_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  isplit = j1

  return
end
subroutine r8vec_std ( n, a, std )

!*****************************************************************************80
!
!! R8VEC_STD returns the standard deviation of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard deviation of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n - 1 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) std

  if ( n < 2 ) then

    std = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    std = sum ( ( a(1:n) - mean )**2 )

    std = sqrt ( std / real ( n - 1, kind = 8 ) )

  end if

  return
end
subroutine r8vec_step ( x0, n, x, fx )

!*****************************************************************************80
!
!! R8VEC_STEP evaluates a unit step function.
!
!  Discussion:
!
!    F(X) = 0 if X < X0
!           1 if     X0 <= X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, the location of the jump.
!
!    Input, integer ( kind = 4 ) N, the number of argument values.
!
!    Output, real ( kind = 8 ) X(N), the arguments.
!
!    Output, real ( kind = 8 ) FX(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0

  where ( x < x0 )
    fx = 0.0D+00
  else where
    fx = 1.0D+00
  end where

  return
end
subroutine r8vec_stutter ( n, a, m, am )

!*****************************************************************************80
!
!! R8VEC_STUTTER makes a "stuttering" copy of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
!    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input vector.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Input, integer ( kind = 4 ) M, the "stuttering factor".
!
!    Output, real ( kind = 8 ) AM(M*N), the stuttering vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) am(m*n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo

  do i = 1, n
    jlo = m * ( i - 1 ) + 1
    jhi = m *   i
    am(jlo:jhi) = a(i)
  end do

  return
end
function r8vec_sum ( n, a )

!*****************************************************************************80
!
!! R8VEC_SUM returns the sum of the entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      R8VEC_SUM ( N, A ) = SUM ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_sum

  r8vec_sum = sum ( a(1:n) )

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two R8VECs.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!        1.0    2.1    3.2    4.3    5.4
!        6.5    7.6    8.7    9.8   10.9
!       11.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
  end do

  return
end
subroutine r8vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8VEC_UNDEX returns unique sorted indexes for an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of X, which creates arrays INDX and XDNI.  (If all the entries
!    of X are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of X in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector X could be
!    replaced by a compressed vector XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(1:X_UNIQUE_NUM) = X(UNDX(1:X_UNIQUE_NUM)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    X   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11,     1     1
!      2 | 22.     3     5 |    22,     2     2
!      3 | 11.     6     2 |    33,     4     1
!      4 | 33.     9     8 |    55,     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is X(3).
!    XDNI(2) = 5 means that X(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
!    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of data values.
!
!    Input, real ( kind = 8 ) X_VAL(X_NUM), the data values.
!
!    Input, integer ( kind = 4 ) X_UNIQUE_NUM, the number of unique values
!    in X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) x_unique_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(x_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(x_unique_num)
  real ( kind = 8 ) x_val(x_num)
  integer ( kind = 4 ) xdnu(x_num)
!
!  Implicitly sort the array.
!
  call r8vec_sort_heap_index_a ( x_num, x_val, indx )
!
!  Walk through the implicitly sorted array X.
!
  i = 1

  j = 1
  undx(j) = indx(i)

  xdnu(indx(i)) = j

  do i = 2, x_num

    if ( tol < abs ( x_val(indx(i)) - x_val(undx(j)) ) ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_uniform_abvec ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Dimension I ranges from A(I) to B(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), B(N), the lower and upper limits
!    for each dimension.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_ABVEC - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a(i) + ( b(i) - a(i) ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_unique_count ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the unsorted array to examine.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_num
  real ( kind = 8 ) tol

  unique_num = 0

  do i = 1, n

    unique_num = unique_num + 1

    do j = 1, i - 1

      if ( abs ( a(i) - a(j) ) <= tol ) then
        unique_num = unique_num - 1
        exit
      end if

    end do

  end do

  return
end
subroutine r8vec_unique_index ( n, a, tol, unique_index )

!*****************************************************************************80
!
!! R8VEC_UNIQUE_INDEX indexes the unique occurrence of values in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
!    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
!    gathered in order, then
!
!      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the unique index.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do i = 1, n

    if ( unique_index(i) == -1 ) then

      unique_num = unique_num + 1
      unique_index(i) = unique_num

      do j = i + 1, n
        if ( abs ( a(i) - a(j) ) <= tol ) then
          unique_index(j) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine r8vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The variance of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      var ( X(1:n) ) = sum ( ( X(1:n) - mean )^2 ) / ( n - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    variance = sum ( ( a(1:n) - mean )**2 )

    variance = variance / real ( n - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_vector_triple_product ( v1, v2, v3, v )

!*****************************************************************************80
!
!! R8VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.
!
!  Discussion:
!
!    VTRIPLE = V1 x ( V2 x V3 )
!
!    VTRIPLE is a vector perpendicular to V1, lying in the plane
!    spanned by V2 and V3.  The norm of VTRIPLE is the product
!    of the norms of V1, V2 and V3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), three vectors.
!
!    Output, real ( kind = 8 ) V(3), the vector triple product.
!
  implicit none

  real ( kind = 8 ) v(3)
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) v4(3)

  call r8vec_cross_product_3d ( v2, v3, v4 )

  call r8vec_cross_product_3d ( v1, v4, v )

  return
end
subroutine r8vec_write ( n, r, output_file )

!*****************************************************************************80
!
!! R8VEC_WRITE writes an R8VEC to a file.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) R(N), the vector to be written.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) r(n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace' )

  do i = 1, n
    write ( output_unit, '(2x,g16.8)' ) r(i)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8vec_zero ( n, a )

!*****************************************************************************80
!
!! R8VEC_ZERO zeroes out an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, real ( kind = 8 ) A(N), the vector to be zeroed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)

  a(1:n) = 0.0D+00

  return
end
subroutine r8vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! R8VEC2_COMPARE compares two entries in an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    A1(I) A2(I)   A1(J) A2(J)  ISGN
!    -----------   -----------  ----
!    1.0   5.0  <  1.0   6.0     -1
!    1.0   5.0  <  2.0   8.0     -1
!    1.0   5.0  <  9.0   1.0     -1
!    1.0   5.0  =  1.0   5.0      0
!    1.0   5.0  >  0.0   2.0     +1
!    1.0   5.0  >  0.0   5.0     +1
!    1.0   5.0  >  1.0   3.0     +1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the two components of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vectors.
!
!    Input, real ( kind = 8 ) X1(N), X2(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title
  real ( kind = 8 ) x1(n)
  real ( kind = 8 ) x2(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '  ......  ..............  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end
subroutine r8vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_A ascending sorts an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_D descending sorts an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!  Reverse the value of ISGN to effect a descending sort.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec2_sort_heap_index_a ( n, x, y, indx )

!*****************************************************************************80
!
!! R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
!
!    * X(I) < X(J), or
!
!    * X(I) = X(J), and Y(I) < Y(J).
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      ( X(INDX(1:N)), Y(INDX(1:N) ), is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, indx, x )
!      call r8vec_permute ( n, indx, y )
!
!    after which ( X(1:N), Y(1:N) ), is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N),Y(N), pairs of X, Y coordinates of points.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array has coordinates ( X(INDX(I)), Y(INDX(I) ).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yval

  if ( n < 1 ) then
    return
  end if

  call i4vec_indicator ( n, indx )

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      xval = x(indxt)
      yval = y(indxt)

    else

      indxt = indx(ir)
      xval = x(indxt)
      yval = y(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        if ( x(indx(j)) < x(indx(j+1)) .or. &
          ( x(indx(j)) == x(indx(j+1)) .and. y(indx(j)) < y(indx(j+1)) ) ) then
          j = j + 1
        end if

      end if

      if ( xval < x(indx(j)) .or. &
          ( xval == x(indx(j)) .and. yval < y(indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! R8VEC2_SORTED_UNIQUE keeps unique elements in a sorted R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
subroutine r8vec2_sorted_unique_index ( n, a1, a2, unique_num, indx )

!*****************************************************************************80
!
!! R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it should be the
!    case that equal items are stored in adjacent vector locations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the array of N items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
!    Output, integer ( kind = 4 ) INDX(N), contains in entries 1 through
!    UNIQUE_NUM an index array of the unique items.  To build new arrays
!    with no repeated elements:
!      B1(1:UNIQUE_NUM) = A1(INDX(1:UNIQUE_NUM))
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  indx(1) = 1

  do itest = 2, n

    if ( a1(itest-1) /= a1(itest) .or. a2(itest-1) /= a2(itest) ) then

      unique_num = unique_num + 1

      indx(unique_num) = itest

    end if

  end do

  indx(unique_num+1:n) = 0

  return
end
subroutine r8vec2_sum_max_index ( n, a, b, sum_max_index )

!*****************************************************************************80
!
!! R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), B(N), two arrays whose sum
!    is to be examined.
!
!    Output, integer ( kind = 4 ) SUM_MAX_INDEX, the index of the largest
!    entry in A+B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) sum_max
  integer ( kind = 4 ) sum_max_index

  if ( n <= 0 ) then

    sum_max_index = -1

  else

    sum_max_index = 1
    sum_max = a(1) + b(1)

    do i = 2, n
      if ( sum_max < a(i) + b(i) ) then
        sum_max = a(i) + b(i)
        sum_max_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec3_print ( n, a1, a2, a3, title )

!*****************************************************************************80
!
!! R8VEC3_PRINT prints an R8VEC3.
!
!  Discussion:
!
!    An R8VEC3 is a dataset consisting of N triples of R8's, stored
!    as three separate vectors A1, A2, A3.
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
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i8,3g14.6)' ) i, a1(i), a2(i), a3(i)
  end do

  return
end
subroutine roots_to_r8poly ( n, x, c )

!*****************************************************************************80
!
!! ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of roots specified.
!
!    Input, real ( kind = 8 ) X(N), the roots.
!
!    Output, real ( kind = 8 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!    On return, if INDX is
!    * greater than 0,
!      > interchange items I and J;
!      > call again.
!    * less than 0,
!      > compare items I and J;
!      > set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      > call again.
!    * equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
