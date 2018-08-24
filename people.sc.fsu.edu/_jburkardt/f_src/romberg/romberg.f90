function romberg_trap ( a, b, f, tol, m )

!*****************************************************************************80
!
!! ROMBERG_TRAP approximates an integral by extrapolating the trapezoidal rule.
!
!  Discussion:
!
!    This routine computes a sequence of integral estimates involving the 
!    trapezoid rule.  Extrapolation of successive trapezoidal estimates 
!    produces an estimate with a higher rate of convergence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Input, real ( kind = 8 ) TOL, the tolerance.  When two successive integral
!    estimates differ by less than this value, the iteration will halt.
!
!    Output, integer ( kind = 4 ) M, the number of trapezoidal estimates
!    that were required.
!
!    Output, real ( kind = 8 ) ROMBERG_TRAP, the integral estimate.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) q_new
  real ( kind = 8 ) r
  real ( kind = 8 ) r_new
  real ( kind = 8 ) romberg_trap
  real ( kind = 8 ) tol
  real ( kind = 8 ) trapezoid_refine

  q_new = 0.0D+00
  r_new = 0.0D+00
  m = 1

  do

    q = q_new
    q_new = trapezoid_refine ( a, b, m, f, q )

    if ( m == 1 ) then
      r_new = q_new
    else
      r = r_new
      r_new = ( 4.0D+00 * q_new - q ) / 3.0D+00
      if ( abs ( r_new - r ) .lt. tol * ( 1.0 + abs ( r_new ) ) ) then
        exit
      end if
    end if

    if ( 20 <= m ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROMBERG_TRAP - Fatal error!'
      write ( *, '(a)' ) '  No convergence in 20 iterations.'
      write ( *, '(a)' ) '  The algorithm is halting.'
      stop
    end if

    m = m + 1

  end do

  romberg_trap = r_new

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
!    18 May 2013
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function trapezoid ( a, b, n, f )

!*****************************************************************************80
!
!! TRAPEZOID carries out the composite trapezoid rule.
!
!  Discussion:
!
!    N points are used to divide the interval [a,b] into equal subintervals.
!
!    Over each subinterval, the integral of the function F is estimated
!    by averaging its value at the subinterval endpoints and multiplying
!    by the length of the subinterval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the number of points to use.
!    N must be at least 1.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Output, real ( kind = 8 ) TRAPEZOID, the integral estimate.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) trapezoid
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAPEZOID - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    stop

  else if ( n == 1 ) then

    x = ( a + b ) / 2.0D+00
    value = ( b - a ) * f ( x )

  else

    value = 0.0D+00

    x = a
    value = value + 0.5D+00 * f ( x )

    do i = 2, n - 1
      x = ( real ( n - i,     kind = 8 ) * a   &
          + real (     i - 1, kind = 8 ) * b ) &
          / real ( n     - 1, kind = 8 )
      value = value + f ( x )
    end do

    x = b
    value = value + 0.5D+00 * f ( x )

    value = ( b - a ) * value / real ( n - 1, kind = 8 )

  end if

  trapezoid = value

  return
end
function trapezoid_refine ( a, b, m, f, q )

!*****************************************************************************80
!
!! TRAPEZOID_REFINE carries out a step of trapezoidal refinement.
!
!  Discussion:
!
!    This routine is designed to be an efficient way to carry out a 
!    sequence of integral estimates, using the trapezoidal rule
!    and a nested sequence of evaluation points, in such a way that
!    the function is not re-evaluated unnecessarily.
!
!    The user calls first with M = 1 and Q = 0 to get a 2 point
!    integral estimate.  On the second call, the user sets M = 2,
!    and the input value of Q should be the integral estimate returned
!    on the previous call.  By incrementing M and passing the previous
!    estimate, the user gets a sequence of increasingly improved
!    integral estimates:
!
!    q = 0.0
!    m = 1
!    do
!      q_new = trapezoid_refine ( a, b, m, f, q )
!      if ( satisfied ) then
!        exit
!      end if
!      q = q_new
!      m = m + 1
!    end do
!
!    The number of points used on each step of the iteration is:
!
!    M   N 
!    1   2
!    2   3
!    3   5
!    4   9
!    5  17
!    6  33
!    m   2^(m-1)+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) M, the current level of refinement.
!    On the first call, M should be 1.  On each subsequent call,
!    the user should increment M by 1.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Input, real ( kind = 8 ) Q, the integral estimate return on
!    the previous call.  But on the first call, set Q to zero.
!
!    Output, real ( kind = 8 ) TRAPEZOID_REFINE, the improved 
!    integral estimate.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) trapezoid_refine
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( m < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAPEZOID_REFINE - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of M.'
    stop

  else if ( m == 1 ) then

    value = ( b - a ) / 2.0D+00 * ( f ( a ) + f ( b ) )

  else if ( 2 <= m ) then

    k = 2 ** ( m - 2 )
    value = 0.0D+00
    do i = 1, k
      x = ( real ( 2 * k - 2 * i + 1, kind = 8 ) * a   &
          + real (         2 * i - 1, kind = 8 ) * b ) &
          / real ( 2 * k,             kind = 8 )
      value = value + f ( x )
    end do

    value = 0.5D+00 * q + ( b - a ) * value / real ( 2 * k, kind = 8 )

  end if

  trapezoid_refine = value

  return
end
function trapezoid_tol ( a, b, f, tol, m )

!*****************************************************************************80
!
!! TRAPEZOID_TOL repeatedly refines the trapezoid rule til a tolerance is met.
!
!  Discussion:
!
!    This routine computes a sequence of integral estimates involving the 
!    trapezoid rule, until two successive estimates differ by less than a 
!    user-designated tolerance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Input, real ( kind = 8 ) TOL, the tolerance.  When two successive integral
!    estimates differ by less than this value, the iteration will halt.
!
!    Output, integer ( kind = 4 ) M, the number of trapezoidal estimates
!    that were required.
!
!    Output, real ( kind = 8 ) TRAPEZOID_TOL, the integral estimate.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) q_new
  real ( kind = 8 ) tol
  real ( kind = 8 ) trapezoid_refine
  real ( kind = 8 ) trapezoid_tol

  q = 0.0D+00
  m = 1

  do

    q_new = trapezoid_refine ( a, b, m, f, q )

    if ( 1 < m .and. abs ( q - q_new ) .lt. tol ) then
      exit
    end if

    if ( 20 <= m ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRAPEZOID_TOL - Fatal error!'
      write ( *, '(a)' ) '  No convergence in 20 iterations.'
      write ( *, '(a)' ) '  The algorithm is halting.'
      stop
    end if

    q = q_new
    m = m + 1

  end do

  trapezoid_tol = q_new

  return
end

