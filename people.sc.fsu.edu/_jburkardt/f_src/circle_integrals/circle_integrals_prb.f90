program main

!*****************************************************************************80
!
!! MAIN is the main program for CIRCLE_INTEGRALS_PRB.
!
!  Discussion:
!
!    CIRCLE_INTEGRALS_PRB tests the CIRCLE_INTEGRALS library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_INTEGRALS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CIRCLE_INTEGRALS library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_INTEGRALS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses CIRCLE01_SAMPLE to compare exact and estimated monomial integrals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 4192

  real ( kind = 8 ) circle01_length
  integer ( kind = 4 ) e(m)
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use CIRCLE01_SAMPLE to compare exact and'
  write ( *, '(a)' ) '  estimated integrals along the circumference '
  write ( *, '(a)' ) '  of the unit circle in 2D.'
!
!  Get sample points.
!
  seed = 123456789
  call circle01_sample ( n, seed, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X, Y exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  If any exponent is odd, the integral is zero.'
  write ( *, '(a)' ) '  We restrict this test to randomly chosen even exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey     MC-Estimate           Exact      Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 5, seed, e )

    e(1:m) = e(1:m) * 2

    call monomial_value ( m, n, e, x, value )

    result = circle01_length ( ) * sum ( value(1:n) ) &
      / real ( n, kind = 8 )
    call circle01_monomial_integral ( e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end
