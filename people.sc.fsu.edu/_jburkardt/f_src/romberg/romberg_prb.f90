program main

!*****************************************************************************80
!
!! MAIN is the main program for ROMBERG_PRB.
!
!  Discussion:
!
!    ROMBERG_PRB tests the ROMBERG library.
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
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROMBERG_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ROMBERG library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROMBERG_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of TRAPEZOID.
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
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) e01
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: f01
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) trapezoid

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TRAPEZOID carries out an iterative procedure to estimate'
  write ( *, '(a)' ) '  an integral using the composite trapezoid rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N     Integral         Error'
  write ( *, '(a)' ) '             estimate'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 2.0D+00
  e = e01 ( )
  n = 2
 
  do m = 1, 12
    q = trapezoid ( a, b, n, f01 )
    err = abs ( q - e )
    write ( *, '(2x,i6,2x,g14.6,2x,g14.6)' ) n, q, err
    n = 2 * n - 1
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of TRAPEZOID_REFINE.
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
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) e01
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: f01
  integer ( kind = 8 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) q_new
  real ( kind = 8 ) trapezoid_refine

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  TRAPEZOID_REFINE carries out an iterative procedure'
  write ( *, '(a)' ) '  to estimate an integral using the composite '
  write ( *, '(a)' ) '  trapezoid rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     M       N     Integral         Error'
  write ( *, '(a)' ) '                   estimate'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 2.0D+00
  e = e01 ( )
  n = 2
  q = 0.0D+00

  do m = 1, 12
    q_new = trapezoid_refine ( a, b, m, f01, q )
    err = abs ( q_new - e )
    write ( *, '(2x,i4,2x,i6,2x,g14.6,2x,g14.6)' ) m, n, q_new, err
    q = q_new
    n = 2 * n - 1
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of TRAPEZOID_TOL.
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
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) e01
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: f01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) tol
  real ( kind = 8 ) trapezoid_tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  TRAPEZOID_TOL estimates an integral by refining a'
  write ( *, '(a)' ) '  trapezoid estimate until a tolerance is met. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TOL            Integral        Error        M'
  write ( *, '(a)' ) '                      estimate'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 2.0D+00
  e = e01 ( )
  tol = 1.0D+00

  do i = 1, 8
    q = trapezoid_tol ( a, b, f01, tol, m )
    err = abs ( q - e )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,i2)' ) tol, q, err, m
    tol = tol / 10.0D+00
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates the use of ROMBERG_TRAP.
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
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) e01
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: f01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) tol
  real ( kind = 8 ) romberg_trap

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  ROMBERG_TRAP estimates an integral by Romberg'
  write ( *, '(a)' ) '  extrapolation from trapezoid estimates. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TOL            Integral        Error      M'
  write ( *, '(a)' ) '                      estimate'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 2.0D+00
  e = e01 ( )
  tol = 1.0D+00

  do i = 1, 8
    q = romberg_trap ( a, b, f01, tol, m )
    err = abs ( q - e )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,i2)' ) tol, q, err, m
    tol = tol / 10.0D+00
  end do

  return
end
function e01 ( )

!*****************************************************************************80
!
!! E01 returns the value of the first test integral.
!
!  Discussion:
!
!    Integral ( 0 <= x <= 2 ) x^4 log ( x + sqrt ( x^2+1 ) ) dx 
!      = ( 8 - 40 Sqrt ( 5 ) + 480 ArcSinh(2) ) / 75
!      = 8.1533641198111650205...
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
!    Output, real E01, the value of the integral.
!
  implicit none

  real ( kind = 8 ) e01

  e01 = ( 8.0D+00 - 40.0D+00 * sqrt ( 5.0D+00 ) &
    + 480.0D+00 * asinh ( 2.0D+00 ) ) / 75.0D+00

  return
end
function f01 ( x )

!*****************************************************************************80
!
!! F01 evaluates the first test integrand.
!
!  Discussion:
!
!    Integral ( 0 <= x <= 2 ) x^4 log ( x + sqrt ( x^2+1 ) ) dx 
!      = ( 8 - 40 Sqrt ( 5 ) + 480 ArcSinh(2) ) / 75
!      = 8.1533641198111650205...
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
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F01, the value of the integrand.
!
  implicit none

  real ( kind = 8 ) f01
  real ( kind = 8 ) x

  f01 = x ** 4 * log ( x + sqrt ( x * x + 1.0D+00 ) )

  return
end


