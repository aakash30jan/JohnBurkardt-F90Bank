program main

!*****************************************************************************80
!
!! MAIN is the main program for R4LIB_PRB.
!
!  Discussion:
!
!    R4LIB_PRB tests the R4LIB library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4LIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the R4LIB library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )

  call test010 ( )
  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test018 ( )
  call test019 ( )

  call test020 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test0235 ( )
  call test027 ( )
  call test028 ( )

  call test12555 ( )

  call test137 ( )

  call test150 ( )
  call test1504 ( )
  call test1505 ( )
  call test151 ( )
  call test1515 ( )
  call test153 ( )
  call test154 ( )
  call test155 ( )
  call test156 ( )
  call test157 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4LIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests R4_ABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 June 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_abs
  real ( kind = 4 ) r4_absolute
  real ( kind = 4 ) r4_uniform_ab
  real ( kind = 4 ) :: r4_hi = 5.0E+00
  real ( kind = 4 ) :: r4_lo = -3.0E+00
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  R4_ABS returns the absolute value of an R4.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r4 = r4_uniform_ab ( r4_lo, r4_hi, seed )
    r4_absolute = r4_abs ( r4 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r4, r4_absolute
  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests R4_ATAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 4 ) r4_atan
  integer ( kind = 4 ) test
  real ( kind = 4 ) x
  real ( kind = 4 ), dimension ( test_num ) :: xtest = (/ &
     1.0E+00,  1.0E+00,  0.0E+00, -1.0E+00, &
    -1.0E+00, -1.0E+00,  0.0E+00,  1.0E+00 /)
  real ( kind = 4 ) y
  real ( kind = 4 ), dimension ( test_num ) :: ytest = (/ &
     0.0E+00,  1.0E+00,  1.0E+00,  1.0E+00, &
     0.0E+00, -1.0E+00, -1.0E+00, -1.0E+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  R4_ATAN computes the arc-tangent given Y and X;'
  write ( *, '(a)' ) '  ATAN2 is the system version of this routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             Y          ATAN2(Y,X)    R4_ATAN(Y,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = xtest(test)
    y = ytest(test)
    write ( *, '(2x,4g14.6)' ) x, y, atan2 ( y, x ), r4_atan ( y, x )
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests R4_CAS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4_cas
  real ( kind = 4 ) r4_pi
  integer ( kind = 4 ), parameter :: test_num = 12
  integer ( kind = 4 ) test
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  R4_CAS evaluates the casine of a number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X           R8_CAS ( X )'
  write ( *, '(a)' ) ' '
  do test = 0, test_num
    x = r4_pi ( ) * real ( test, kind = 4 ) / real ( test_num, kind = 4 )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, r4_cas ( x )
  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests R4_CEILING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) r4_ceiling
  integer ( kind = 4 ) ival
  real ( kind = 4 ) rval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  R4_CEILING rounds a value up.'
  write ( *, '(a)' ) ' '

  do i = -6, 6
    rval = real ( i, kind = 4 ) / 5.0E+00
    ival = r4_ceiling ( rval )
    write ( *, '(2x,g14.6,2x,i8)' ) rval, ival
  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests R4_DIFF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 15

  integer ( kind = 4 ) ndig
  real ( kind = 4 ) r4_diff
  integer ( kind = 4 ) test
  real ( kind = 4 ) x
  real ( kind = 4 ), dimension ( test_num ) :: y_test = (/ &
    0.0625E+00, 0.125E+00, 0.25E+00, 0.50E+00,  0.874E+00, &
    0.876E+00,  0.90E+00,  0.95E+00, 0.99E+00,  1.0E+00, &
    1.01E+00,   1.05E+00,  1.10E+00, 3.0E+00,  10.0E+00 /)
  real ( kind = 4 ) y

  ndig = 3
  x = 1.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  R4_DIFF computes a difference X-Y to a given'
  write ( *, '(a)' ) '    number of binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  For this test, we use ', ndig, ' binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y       X-Y     R4_DIFF(X,Y)'
  write ( *, '(a)' ) ' '
  do test = 1, test_num
    y = y_test(test)
    write ( *, '(4f10.4)' ) x, y, x-y, r4_diff ( x, y, ndig )
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests R4_DIGIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdig = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) digit(-2:maxdig)
  integer ( kind = 4 ) idigit
  real ( kind = 4 ) r4_pi
  real ( kind = 4 ) x

  x = r4_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  R4_DIGIT extracts decimal digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Here, we get digits of ', x
  write ( *, '(a)' ) ' '

  do idigit = -2, maxdig
    call r4_digit ( x, idigit, digit(idigit) )
  end do

  write ( *, '(2x,25i3)' ) ( i, i = -2, maxdig )
  write ( *, '(2x,25i3)' ) digit(-2:maxdig)

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests R4_EPSILON
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4_epsilon
  real ( kind = 4 ) r
  real ( kind = 4 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  R4_EPSILON produces the R4 machine precision.'
  write ( *, '(a)' ) ' '

  r = r4_epsilon ( )
  write ( *, '(a,g14.6)' ) '  R = R4_EPSILON()   = ', r

  s = ( 1.0E+00 + r ) - 1.0E+00
  write ( *, '(a,g14.6)' ) '  ( 1 + R ) - 1     =  ', s

  s = ( 1.0E+00 + ( r / 2.0E+00 ) ) - 1.0E+00
  write ( *, '(a,g14.6)' ) '  ( 1 + (R/2) ) - 1 =  ', s

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests R4_FRACTION.
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
  implicit none

  real ( kind = 4 ) fraction
  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_fraction
  real ( kind = 4 ) r4_uniform_ab
  real ( kind = 4 ) :: r4_hi = 5.0E+00
  real ( kind = 4 ) :: r4_lo = -3.0E+00
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  R4_FRACTION returns the fraction part of an R4.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r4 = r4_uniform_ab ( r4_lo, r4_hi, seed )
    fraction = r4_fraction ( r4 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r4, fraction
  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests R4_HUGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4_huge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  R4_HUGE returns a "huge" R4;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '    R4_HUGE ( ) =      ', r4_huge ( )
  write ( *, '(a,g24.16)' ) '    HUGE ( 1.0E+00 ) = ', huge ( 1.0E+00 )

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests R4_LOG_2.
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
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 18

  real ( kind = 4 ) r4_log_2
  integer ( kind = 4 ) test
  real ( kind = 4 ) x
  real ( kind = 4 ), dimension(test_num) :: x_test = (/ &
    0.0E+00,  1.0E+00,  2.0E+00,   3.0E+00,  9.0E+00, &
   10.0E+00, 11.0E+00, 99.0E+00, 101.0E+00, -1.0E+00, &
   -2.0E+00, -3.0E+00, -9.0E+00,   0.5E+00,  0.33E+00, &
    0.25E+00, 0.20E+00, 0.01E+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  R4_LOG_2 computes the logarithm base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, R4_LOG_2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2g14.6 )' ) x, r4_log_2 ( x )
  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests R4_LOG_B.
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
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 4 ) b
  real ( kind = 4 ), dimension(test_num) :: b_test = (/ &
    2.0E+00, 3.0E+00, 4.0E+00, 5.0E+00, 6.0E+00, &
    7.0E+00, 8.0E+00, 16.0E+00, 32.0E+00, 256.0E+00 /)
  real ( kind = 4 ) r4_log_b
  integer ( kind = 4 ) test
  real ( kind = 4 ) x

  x = 16.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  R48_LOG_B computes the logarithm base B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, B, R4_LOG_B'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    b = b_test(test)

    write ( *, '( 2x,3g14.6, i12 )' ) x, b, r4_log_b ( x, b )

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests R8_MANT.
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
  implicit none

  integer ( kind = 4 ) is
  integer ( kind = 4 ) l
  real ( kind = 4 ) r
  real ( kind = 4 ) x

  x = -314.159E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  R4_MANT decomposes a value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number to be decomposed:'
  write ( *, '(2x,g14.6)' ) x

  call r4_mant ( x, is, r, l )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,g14.6,a,i8)' ) &
    '  R4_MANT: X = ', is, ' * ', r, ' * 2**', l

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests R4_MOD.
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
  implicit none

  real ( kind = 4 ) r4_mod
  real ( kind = 4 ) r4_uniform_ab
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x
  real ( kind = 4 ) :: x_hi =  10.0E+00
  real ( kind = 4 ) :: x_lo = -10.0E+00
  real ( kind = 4 ) y
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  R4_MOD returns the remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Y     MOD(X,Y)    R4_MOD(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    x = r4_uniform_ab ( x_lo, x_hi, seed )
    y = r4_uniform_ab ( x_lo, x_hi, seed )

    z1 =    mod ( x, y )
    z2 = r4_mod ( x, y )

    write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests R4_MODP.
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
  implicit none

  real ( kind = 4 ) r4_modp
  real ( kind = 4 ) r4_uniform_ab
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x
  real ( kind = 4 ) :: x_hi =  10.0E+00
  real ( kind = 4 ) :: x_lo = -10.0E+00
  real ( kind = 4 ) y
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  R4_MODP returns the remainder after division.'
  write ( *, '(a)' ) '  Unlike the FORTRAN MOD, R4_MODP ( X, Y ) is positive.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Y      MOD(X,Y)  R4_MODP(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    x = r4_uniform_ab ( x_lo, x_hi, seed )
    y = r4_uniform_ab ( x_lo, x_hi, seed )

    z1 =    mod  ( x, y )
    z2 = r4_modp ( x, y )

    write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests R4_NINT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) r4_nint
  real ( kind = 4 ) r4_uniform_ab
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  R4_NINT produces the nearest integer to an R4.'
  write ( *, '(a)' ) ' '

  b = -10.0E+00
  c = +10.0E+00

  do test = 1, test_num
    x = r4_uniform_ab ( b, c, seed )
    write ( *, '(2x,f10.4,2x,i8)' ) x, r4_nint ( x )
  end do

  return;
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests R4_NORMAL_01.
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
  implicit none

  real ( kind = 4 ) r4_normal_01
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  R4_NORMAL_01 generates normally distributed'
  write ( *, '(a)' ) '    random values.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = r4_normal_01 ( seed )
    write ( *, '(2x,g14.6)' ) x

  end do

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests R4_PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) four
  real ( kind = 4 ) one
  real ( kind = 4 ) r4_pi
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2

  four = real ( 4, kind = 4 )
  one = real ( 1, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  R4_PI returns the value of PI.'
  write ( *, '(a)' ) ' '
  v1 = r4_pi ( )
  write ( *, '(a,g24.16)' ) '  R4_PI =     ', v1
  v2 = four * atan ( one )
  write ( *, '(a,g24.16)' ) '  4*atan(1) = ', v2

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests R4_POWER.
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
  implicit none

  real ( kind = 4 ) r4_power
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p
  real ( kind = 4 ) r
  real ( kind = 4 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  R4_POWER computes R**P.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R          P       R**P'
  write ( *, '(a)' ) ' '

  do i = -5, 5

    r = 2.0E+00
    p = i
    value = r4_power ( r, p )
    write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, value

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests R4_POWER_FAST.
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
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) mults
  integer ( kind = 4 ) p
  real ( kind = 4 ) r
  real ( kind = 4 ) rp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  R4_POWER_FAST computes R**P, economizing on'
  write ( *, '(a)' ) '    multiplications.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R          P       R**P       Mults'
  write ( *, '(a)' ) ' '

  do i = -10, 40

    r = 2.0E+00
    p = i
    call r4_power_fast ( r, p, rp, mults )
    write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, rp, mults

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests R4_ROUND2.
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
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  real ( kind = 4 ) r4_pi
  real ( kind = 4 ) x
  real ( kind = 4 ) xround

  x = r4_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  R4_ROUND2 rounds a number to a'
  write ( *, '(a)' ) '    specified number of base 2 digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    nplace = i
    call r4_round2 ( nplace, x, xround )
    write ( *, '(2x,i8,g14.6)' ) i, xround
  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests R4_ROUNDB.
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
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  real ( kind = 4 ) r4_pi
  real ( kind = 4 ) x
  real ( kind = 4 ) xround

  base = 3
  x = r4_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  R4_ROUNDB rounds a number to a '
  write ( *, '(a)' ) '    specified number of base BASE digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here, we will use BASE = ',base
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    nplace = i
    call r4_roundb ( base, nplace, x, xround )
    write ( *, '(2x,i8,g14.6)' ) i, xround
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Try with a negative base:'
  x = 121.0E+00
  base = -3
  nplace = 3
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Input quantity is X = ', x
  write ( *, '(a,i8)' ) '  to be rounded in base ', base

  do nplace = 1, 5

    call r4_roundb ( base, nplace, x, xround )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,g14.6)' ) '  Output value to ', nplace, &
      ' places is ', xround

  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests R4_ROUNDX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r4_pi
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x
  real ( kind = 4 ) xround

  seed = 123456789
  x = r4_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  R4_ROUNDX rounds a number to a '
  write ( *, '(a)' ) '    specified number of decimal digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g16.10)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    nplace = i
    call r4_roundx ( nplace, x, xround )
    write ( *, '(2x,i8,g16.10)' ) i, xround
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on random values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  X     XROUND'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    x = r4_uniform_01 ( seed )

    write ( *, '(a)' ) ' '

    do nplace = 0, 10, 2
      call r4_roundx ( nplace, x, xround )
      write ( *, '(2x,i8,2x,g16.10,2x,g16.10)' ) nplace, x, xround
    end do

  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests R4_SIGN and R4_SIGN3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_sign
  real ( kind = 4 ) r4_sign3 
  real ( kind = 4 ), parameter, dimension ( test_num ) :: r4_test = (/ &
    -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 /)
  real ( kind = 4 ) s1
  real ( kind = 4 ) s2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  R4_SIGN returns the sign of an R4.'
  write ( *, '(a)' ) '  R4_SIGN3 returns the three-way sign of an R4.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      R4    R4_SIGN(R4)  R4_SIGN3(R4)'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    r4 = r4_test(test)
    s1 = r4_sign ( r4 )
    s2 = r4_sign3 ( r4 )
    write ( *, '(2x,f8.4,2x,f8.0,2x,f8.0)' ) r4, s1, s2
  end do

  return
end
subroutine test0235 ( )

!*****************************************************************************80
!
!! TEST0235 tests R4_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) x
  real ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0235'
  write ( *, '(a)' ) '  R4_SWAP swaps two reals.'

  x = 1.0E+00
  y = 3.141592653589793E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before swapping:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    X = ', x
  write ( *, '(a,g14.6)' ) '    Y = ', y

  call r4_swap ( x, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After swapping:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    X = ', x
  write ( *, '(a,g14.6)' ) '    Y = ', y

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests R4_UNIFORM_01.
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
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_old
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  R4_UNIFORM_01 produces a sequence of random values.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random seed ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED   R4_UNIFORM_01(SEED)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    seed_old = seed
    x = r4_uniform_01 ( seed )
    write ( *, '(2x,i12,2x,g14.6)' ) seed, x
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the sequence can be restarted.'
  write ( *, '(a)' ) '  Set the seed back to its original value, and see that'
  write ( *, '(a)' ) '  we generate the same sequence.'

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED   R4_UNIFORM_01(SEED)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_old = seed
    x = r4_uniform_01 ( seed )
    write ( *, '(2x,i12,2x,g14.6)' ) seed, x
  end do

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests R4_UNIFORM_01
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
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) mean
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01
  real ( kind = 4 ) variance
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  R4_UNIFORM_01 samples a uniform random'
  write ( *, '(a)' ) '  distribution in [0,1].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Starting with seed = ', seed

  do i = 1, n
    x(i) = r4_uniform_01 ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few values:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  call r4vec_mean ( n, x, mean )

  call r4vec_variance ( n, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,g14.6)' ) '  Minimum value was ', minval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Maximum value was ', maxval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  return
end
subroutine test12555 ( )

!*****************************************************************************80
!
!! TEST12555 tests R4VEC_INDICATOR0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 4 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12555'
  write ( *, '(a)' ) '  R4VEC_INDICATOR0 returns an indicator vector.'

  call r4vec_indicator0 ( n, a )

  call r4vec_print ( n, a, '  The indicator0 vector:' )

  return
end
subroutine test137 ( )

!*****************************************************************************80
!
!! TEST137 tests R4VEC_SORT_BUBBLE_A.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 4 ) a(n)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137'
  write ( *, '(a)' ) '  R4VEC_SORT_BUBBLE_A ascending sorts a R4VEC.'

  b = 0.0E+00
  c = 3.0E+00 * real ( n, kind = 4 )
  seed = 123456789

  call r4vec_uniform_ab ( n, b, c, seed, a )

  call r4vec_print_some ( n, a, 1, 10, '  Original array:' )

  call r4vec_sort_bubble_a ( n, a )

  call r4vec_print_some ( n, a, 1, 10, '  Ascending sorted array:' )

  return
end
subroutine test150 ( )

!*****************************************************************************80
!
!! TEST150 tests R4VEC_SORTED_UNIQUE_HIST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: unique_max = 30
  integer ( kind = 4 ), parameter :: n = 30

  real ( kind = 4 ) a(n)
  integer ( kind = 4 ) acount(unique_max)
  real ( kind = 4 ) auniq(unique_max)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed
  real ( kind = 4 ), parameter :: tol = 0.25E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150'
  write ( *, '(a)' ) '  R4VEC_SORTED_UNIQUE_HIST makes a historgram of'
  write ( *, '(a)' ) '  the unique entries in a real vector.'

  b = 0.0E+00
  c = real ( n, kind = 4 )
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random number seed ', seed

  call r4vec_uniform_ab ( n, b, c, seed, a )

  a(1:n) = real ( int ( a(1:n) ), kind = 4 ) + 0.5E+00

  call r4vec_print ( n, a, '  Unsorted array:' )

  call r4vec_sort_bubble_a ( n, a )

  call r4vec_print ( n, a, '  Ascending sorted array:' )

  call r4vec_sorted_unique_hist ( n, a, tol, unique_max, unique_num, &
    auniq, acount )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  R4VEC_SORTED_UNIQUE_HIST counts ' , unique_num, &
    ' unique entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Value  Multiplicity'
  write ( *, '(a)' ) ' '
  do i = 1, unique_num
    write ( *, '(2x,i8,2x,g14.6,2x,i8)' ) i, auniq(i), acount(i)
  end do

  return
end
subroutine test1504 ( )

!*****************************************************************************80
!
!! TEST1504 tests R4VEC_TRANSPOSE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  integer ( kind = 4 ) seed
  real ( kind = 4 ) x(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1504'
  write ( *, '(a)' ) '  R4VEC_TRANSPOSE_PRINT prints an R4VEC "tranposed",'
  write ( *, '(a)' ) '  that is, placing multiple entries on a line.'

  call r4vec_uniform_01 ( n, seed, x )

  call r4vec_transpose_print ( n, x, '  The vector X:' )

  return
end
subroutine test1505 ( )

!*****************************************************************************80
!
!! TEST1505 tests R4VEC_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: x_num = 9

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_epsilon
  real ( kind = 4 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) x_unique_num
  real ( kind = 4 ), dimension ( x_num ) :: x_val = (/ &
    33.0, 55.0, 11.0, 11.0, 55.0, 33.0, 22.0, 22.0, 11.0 /)
  integer ( kind = 4 ) xdnu(x_num)
  real ( kind = 4 ), allocatable, dimension ( : ) :: xu_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1505'
  write ( *, '(a)' ) '  R4VEC_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique elements of an (unsorted) R4VEC,'
  write ( *, '(a)' ) '  and a map from the original vector to the (implicit)'
  write ( *, '(a)' ) '  vector of sorted unique elements.'

  call r4vec_print ( x_num, x_val, '  The vector X:' )

  tol = r4_epsilon ( )
  call r4vec_unique_count ( x_num, x_val, tol, x_unique_num )

  allocate ( undx(1:x_unique_num) )
  allocate ( xu_val(1:x_unique_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tolerance for equality is ', tol
  write ( *, '(a,i8)' ) '  Number of unique entries in X is ', x_unique_num

  call r4vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to list the unique elements of X'
  write ( *, '(a)' ) '  in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX   X(UNDX)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), x_val(undx(i))
  end do

  xu_val(1:x_unique_num) = x_val(undx(1:x_unique_num))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to created XU, a copy of X'
  write ( *, '(a)' ) '  containing only the unique elements, in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX     XU(I)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), xu_val(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU can be used to match each element of X with one of the'
  write ( *, '(a)' ) '  unique elements'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU    X(I)       XU(XDNU(I))'
  write ( *, '(a)' ) ' '

  do i = 1, x_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1,2x,f12.1)' ) i, xdnu(i), x_val(i), xu_val(xdnu(i))
  end do

  deallocate ( undx )
  deallocate ( xu_val )

  return
end
subroutine test151 ( )

!*****************************************************************************80
!
!! TEST151 tests R4VEC_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 4 ), parameter :: b = 10.0E+00
  real ( kind = 4 ), parameter :: c = 20.0E+00
  real ( kind = 4 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151'
  write ( *, '(a)' ) '  R4VEC_UNIFORM_AB returns a random R4VEC '
  write ( *, '(a)' ) '  with entries in a given range [ B, C ]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this problem:'
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,g14.6)' ) '  C = ', c

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Input SEED = ', seed

    call r4vec_uniform_ab ( n, b, c, seed, r )

    call r4vec_print_some ( n, r, 1, 10, '  Random vector:' )

  end do

  return
end
subroutine test1515 ( )

!*****************************************************************************80
!
!! TEST1515 tests R4VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 4 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1515'
  write ( *, '(a)' ) '  R4VEC_UNIFORM_01 returns a random R4VEC '
  write ( *, '(a)' ) '  with entries in [0,1].'

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Input SEED = ', seed

    call r4vec_uniform_01 ( n, seed, r )

    call r4vec_print_some ( n, r, 1, 10, '  Random vector:' )

  end do

  return
end
subroutine test153 ( )

!*****************************************************************************80
!
!! TEST153 tests R4VEC2_SORT_A and R4VEC2_SORT_D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  For a pair of R4VEC''s:'
  write ( *, '(a)' ) '  R4VEC2_SORT_A ascending sorts;'
  write ( *, '(a)' ) '  R4VEC2_SORT_D descending sorts;'

  b = 1.0E+00
  c = 3.0E+00
  seed = 123456789

  call r4vec_uniform_ab ( n, b,  c, seed, a1 )

  b = 5.0E+00
  c = 10.0E+00

  call r4vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r4vec2_print ( n, a1, a2, '  The pair of arrays:' )

  call r4vec2_sort_a ( n, a1, a2 )

  call r4vec2_print ( n, a1, a2, '  Arrays after ascending sort:' )

  call r4vec2_sort_d ( n, a1, a2 )

  call r4vec2_print ( n, a1, a2, '  Arrays after descending sort:' )

  return
end
subroutine test154 ( )

!*****************************************************************************80
!
!! TEST154 tests R4VEC2_SORT_HEAP_INDEX_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed
  real ( kind = 4 ) temp
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  R4VEC2_SORT_HEAP_INDEX_A creates a sort index'
  write ( *, '(a)' ) '  for an (X,Y) array.'

  do i = 1, n

    x(i) = real ( i4_uniform_ab ( 0, n, seed ), kind = 4 ) / real ( n, kind = 4 )
    y(i) = real ( i4_uniform_ab ( 0, n, seed ), kind = 4 ) / real ( n, kind = 4 )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The unsorted array:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,6x,2g14.6)' ) i, x(i), y(i)
  end do

  call r4vec2_sort_heap_index_a ( n, x, y, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After sorting:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,2g14.6)' ) i, indx(i), x(i), y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), X(INDX(I)), Y(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,2g14.6)' ) i, indx(i), x(indx(i)), y(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R4VEC_PERMUTE carries out the permutation.'

  call r4vec_permute ( n, indx, x )
  call r4vec_permute ( n, indx, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,6x,2g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test155 ( )

!*****************************************************************************80
!
!! TEST155 tests R4VEC2_SORTED_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  For a pair of R4VEC''s:'
  write ( *, '(a)' ) '  R4VEC2_SORTED_UNIQUE counts unique entries.'

  b = 1.0E+00
  c = 3.0E+00
  seed = 123456789

  call r4vec_uniform_ab ( n, b, c, seed, a1 )

  b = 5.0E+00
  c = 10.0E+00

  call r4vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r4vec2_print ( n, a1, a2, '  The pair of arrays:' )

  call r4vec2_sort_a ( n, a1, a2 )

  call r4vec2_print ( n, a1, a2, '  Arrays after ascending sort:' )

  call r4vec2_sorted_unique ( n, a1, a2, unique_num )

  call r4vec2_print ( unique_num, a1, a2, '  UNIQed array:' )

  return
end
subroutine test156 ( )

!*****************************************************************************80
!
!! TEST156 tests R4VEC2_SORTED_UNIQUE_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  For a pair of R4VEC''s:'
  write ( *, '(a)' ) '  R4VEC2_SORTED_UNIQUE_INDEX indexes unique entries.'

  b = 1.0E+00
  c = 3.0E+00
  seed = 123456789

  call r4vec_uniform_ab ( n, b, c, seed, a1 )

  b = 5.0E+00
  c = 10.0E+00

  call r4vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r4vec2_sort_a ( n, a1, a2 )

  call r4vec2_print ( n, a1, a2, '  Sorted arrays:' )

  call r4vec2_sorted_unique_index ( n, a1, a2, unique_num, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of unique elements is ', unique_num

  call i4vec_print ( unique_num, indx, '  Index of Unique Elements:' )

  call r4vec_index_order ( unique_num, a1, indx )
  call r4vec_index_order ( unique_num, a2, indx )

  call r4vec2_print ( unique_num, a1, a2, '  After Indexed Nonunique Deletion.' )

  return
end
subroutine test157 ( )

!*****************************************************************************80
!
!! TEST157 tests R4VEC2_SUM_MAX_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157'
  write ( *, '(a)' ) '  For a pair of R4VEC''s:'
  write ( *, '(a)' ) '  R4VEC2_SUM_MAX_INDEX: index of the sum vector'
  write ( *, '(a)' ) '  with maximum value.'

  b = 0.0E+00
  c = 10.0E+00
  seed = 123456789

  call r4vec_uniform_ab ( n, b, c, seed, a1 )

  b = 0.0E+00
  c = 5.0E+00

  call r4vec_uniform_ab ( n, b, c, seed, a2 )

  call r4vec2_print ( n, a1, a2, '  The pair of vectors:' )

  call r4vec2_sum_max_index ( n, a1, a2, ival )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Index of maximum in A+B: ', ival

  return
end

