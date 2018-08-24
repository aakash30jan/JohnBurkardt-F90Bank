program main

!*****************************************************************************80
!
!! MAIN is the main program for E_VERSUS_D.
!
!  Discussion:
!
!    Real numerical constants in FORTRAN must have a precision.
!
!    A constant such as "4.0" will by default have single precision, sometimes
!    thought of as "8 digit precision".
!
!    A constant such as "1.2E+1" will explicitly request single precision, by
!    virtue of the use of the "E" exponent marker.  Note that this has the
!    often unexpected effect of TRUNCATING the numeric data to 8 digits,
!    before anything else is done with it.
!
!    A constant such as "3.4D+2" will explicitly request double precision,
!    by virtue of the occurrence of the "D" exponent marker, sometimes referred 
!    to as "16 digit precision".  Thus, as many as 16 digits may be usefully
!    specified in such a constant.
!
!    The point is that, while either format will allow you to list as many digits
!    as you like, the single precision format will essentially ignore all but
!    the first 8 significant digits.  This is true regardless of the type of
!    any arithmetic expression in which the constant occurs, or of the type of
!    any variable into which the constant is to be stored.
!
!    Thus, even if X and Y are double precision variables, the quantity
!      X = sin ( 3.14159265358979323846264338 * Y )
!    is computed as though the formula read
!      X = sin ( 3.14159265                   * Y ) ( 8 digits of PI)
!    and this is also true for the version
!      X = sin ( 3.14159265358979323846264338E+00 * Y )
!    while
!      X = sin ( 3.14159265358979323846264338D+00 * Y )
!    will be computed as though it read
!      X = sin ( 3.141592653589793 * Y ) (16 digits of PI )
!
!    This is a "feature" of FORTRAN that has existed since the beginning,
!    and one which occurs in no other language that I am aware of.
!
!    OK, now you've been warned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2016
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'E_VERSUS_D:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Demonstrating that FORTRAN constants are single precision'
  write ( *, '(a)' ) '  by default, no matter how many digits you type.'
  write ( *, '(a)' ) '  The only cure is to use the "D" exponent marker.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'E_VERSUS_D:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares storing constants using E and D exponents.
!
!  Discussion:
!
!    X1 and X2 are both double precision variables.
!
!    We initialize them with almost identical statements:
!
!      x1 = 1.2345678901234567890E+00
!      x2 = 1.2345678901234567890D+00
!
!    But X1 and X2 do not receive the same values, which we can
!    verify by printing them, or by computing their difference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2016
!
!  Author:
!
!    John Burkardt.
!
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  X1 and X2 are both double precision variables.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We initialize them with almost identical statements:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    x1 = 1.2345678901234567890E+00'
  write ( *, '(a)' ) '    x2 = 1.2345678901234567890D+00'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  But X1 and X2 do not receive the same values, which we can'
  write ( *, '(a)' ) '  verify by printing them, or by computing their difference.'

  x1 = 1.2345678901234567890E+00
  x2 = 1.2345678901234567890D+00

  x3 = x1 - x2

  write ( *, '(a)' ) ''
  write ( *, '(a,g30.20)' ) '  X1           = ', x1
  write ( *, '(a,g30.20)' ) '  X2           = ', x2
  write ( *, '(a,g30.20)' ) '  X3 = X1 - X2 = ', x3

  return
end

