program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS0_PRB.
!
!  Discussion:
!
!    BLAS0_PRB tests the BLAS0 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BLAS0_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS0 library.'

  call test01 ( )
  call test015 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BLAS0_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R4_ABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_abs
  real ( kind = 4 ) r4_absolute
  real ( kind = 4 ) r4_hi
  real ( kind = 4 ) r4_lo
  real ( kind = 4 ) r4_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  r4_hi = 5.0E+00
  r4_lo = -3.0E+00
  seed = 123456789
  test_num = 10

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R4_ABS returns the absolute value of an R4.'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    r4 = r4_uniform_ab ( r4_lo, r4_hi, seed )
    r4_absolute = r4_abs ( r4 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r4, r4_absolute
  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests R4_SIGN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 4 ) r4_sign
  integer ( kind = 4 ) test
  real ( kind = 4 ) x
  real ( kind = 4 ), save :: x_test(test_num) = (/ &
    -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  R4_SIGN returns the sign of a number.'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    x = x_test(test)
    write ( *, '(2x,f8.2,2x,f8.2)' ) x, r4_sign ( x )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8_ABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8
  real ( kind = 8 ) r8_abs
  real ( kind = 8 ) r8_absolute
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  r8_hi = 5.0D+00
  r8_lo = -3.0D+00
  seed = 123456789
  test_num = 10

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8_ABS returns the absolute value of an R8.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  X     R8_ABS(X)'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
    r8_absolute = r8_abs ( r8 )
    write ( *, '(2x,f10.6,2x,f1.06)' ) r8, r8_absolute
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R8_SIGN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: x_test(test_num) = (/ &
    -1.25D+00, -0.25D+00, 0.0D+00, +0.5D+00, +9.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R8_SIGN returns the sign of a number.'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    x = x_test(test)
    write ( *, '(2x,f8.4,2x,f8.4)' ) x, r8_sign ( x )
  end do

  return
end
