program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS726_PRB.
!
!  Discussion:
!
!    TOMS726_PRB tests the TOMS726 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS726_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS726 library.'

  n = 1
  call test01 ( n )
  n = 4
  call test01 ( n )
  n = 16
  call test01 ( n )
  n = 64
  call test01 ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS726_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n )

!*****************************************************************************80
!
!! TEST01 calls QLAG_R8 to compute a Gauss-Laguerre quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule to be generated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xj(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  QLAG_R8 computes points and weights for'
  write ( *, '(a)' ) '  a Gauss-Laguerre quadrature rule.'
!
!  Compute the rule.
!
  allocate ( x(1:n) )
  allocate ( xj(1:n) )
  allocate ( w(1:n) )

  call qlag_r8 ( n, x, w, ierr )
!
!  Print the rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I      X(i)                      W(i)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i6,2x,g25.18,2x,g25.18)' ) i, x(i), w(i)
  end do
!
!  Test the rule on a few monomials.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       J      Integral(X^J)             Quad(X^J)'
  write ( *, '(a)' ) ' '
  do j = 0, 5
    exact = r8_factorial ( j )
    xj(1:n) = x(1:n) ** j
    quad = dot_product ( w, xj )
    write ( *, '(2x,i6,2x,g25.18,2x,g25.18)' ) j, exact, quad
  end do

  deallocate ( w )
  deallocate ( x )
  deallocate ( xj )

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
subroutine qcheb_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHEB_R8 returns a Gauss-Chebyshev rule.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  real ( kind = 8 ) om2
  real ( kind = 8 ) pi
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  common / common08 / om2

  pi = 4.0D+00 * atan ( 1.0D+00 )

  do k = 1, n
    x(k) = cos ( real ( 2 * k - 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 ) )
    w(k) = pi / ( real ( n, kind = 8 ) * sqrt ( 1.0D+00 - om2 * x(k)**2 ) )
  end do

  return
end
subroutine qjac_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QJAC_R8 returns a Gauss-Jacobi rule.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(41)
  real ( kind = 8 ) b(41)
  real ( kind = 8 ) e(41)
  real ( kind = 8 ) epsma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  common / common05 / a, b, e, epsma

  call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

  w(1:n) = w(1:n) / b(1)

  return
end
subroutine qlag_r8 ( n, x, w, ierr )

!*****************************************************************************80
!
!! QLAG_R8 returns a Gauss-Laguerre rule.
!
!  Modified:
!
!    24 April 2013
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) epsma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
!
!  Get the recursion coefficients.
!
  ipoly = 7
  alpha = 0.0D+00
  beta = 0.0D+00

  call recur_r8 ( n, ipoly, alpha, beta, a, b, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Fatal error!'
    write ( *, '(a,i6)' ) '  RECUR_R8 returned IERR = ', ierr
    stop
  end if
!
!  Generate the Gaussian quadrature formula.
!
  epsma = epsilon ( epsma )
  call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

  if ( ierr == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Warning'
    write ( *, '(a,i6)' ) '  The accuracy request was not achieved.'
  else if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Fatal error!'
    write ( *, '(a,i6)' ) '  GAUSS_R8 returned IERR = ', ierr
    stop
  end if

  return
end
subroutine qchle_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHLE_R8 returns Gauss-Chebyshev or Gauss-Legendre rules.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(81)
  real ( kind = 8 ) b(81)
  real ( kind = 8 ) c
  real ( kind = 8 ) e(81)
  real ( kind = 8 ) epsma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  real ( kind = 8 ) pi
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  pi = 4.0D+00 * atan ( 1.0D+00 )
!
!  Gauss-Chebyshev rule.
!
  if ( i == 1 ) then

    do k = 1, n
      x(k) = cos ( real ( 2 * k - 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 ) )
      w(k) = pi / real ( k, kind = 8 )
    end do
!
!  Gauss-Legendre rule.
!
  else if ( i == 2 ) then

    call recur_r8 ( n, 1, 0.0D+00, 0.0D+00, a, b, ierr )

    call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

    w(1:n) = c * w(1:n)

  end if

  return
end
