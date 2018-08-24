program main

!*****************************************************************************80
!
!! MAIN is the main program for QUADPACK_PRB.
!
!  Discussion:
!
!    QUADPACK_PRB tests the QUADPACK library.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QUADPACK library.'

  call qag_test ( )
  call qagi_test ( )
  call qagp_test ( )
  call qags_test ( )
  call qawc_test ( )
  call qawf_test ( )
  call qawo_test ( )
  call qaws_test ( )
  call qk15_test ( )
  call qk21_test ( )
  call qk31_test ( )
  call qk41_test ( )
  call qk51_test ( )
  call qk61_test ( )
  call qng_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )
 
  stop
end
subroutine qag_test ( )

!*****************************************************************************80
!
!! QAG_TEST tests QAG.
!
!  Discussion:
!
!    QAG is an adaptive automatic integrator using a Gauss-Kronrod rule.
!
!    integrate cos(100*sin(x)) from 0 to pi.
!
!    The exact answer is pi * j0(100), or roughly 0.06278740.
!
!    KEY chooses the order of the integration rule, from 1 to 6.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f02
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: key = 6
  integer ( kind = 4 ) neval
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = 0.06278740E+00

  b = r4_pi

  call qag ( f02, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAG_TEST'
  write ( *, '(a)' ) '  Test QAG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is COS(100*SIN(X))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qagi_test ( )

!*****************************************************************************80
!
!! QAGI_TEST tests QAGI.
!
!  Discussion:
!
!    QAGI is an adaptive quadrature routine for infinite intervals.
!
!    integrate log(x)/(1+100*x*x) from 0 to infinity.
!
!    The exact answer is -pi*log(10)/20 = -0.3616892
!
!    give the type of infinity
!
!    inf=1 means a to infinity
!       -1      -infinity to a
!        2      -infinity to infinity
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f05
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: inf = 1
  integer ( kind = 4 ) neval
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) result
  real ( kind = 4 ) true

  call qagi ( f05, a, inf, epsabs, epsrel, result, abserr, neval, ier )

  true = - r4_pi * log ( 10.0E+00 ) / 20.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAGI_TEST'
  write ( *, '(a)' ) '  Test QAGI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+100*x*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =    Infinity'
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qagp_test ( )

!*****************************************************************************80
!
!! QAGP_TEST tests QAGP.
!
!  Discussion:
!
!    QAGP is an adaptive integrator that can handle singularities
!    of the integrand at user specified points,
!
!    Integrate 
!
!      x**3 * log(abs( (x*x-1)*(x*x-2) )) 
!
!    from 0 to 3.
!
!    The exact answer is 61*log(2)+77*log(7)/4 - 27.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: npts = 2
  integer ( kind = 4 ), parameter :: npts2 = 2 * npts

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 3.0E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f04
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) neval
  real ( kind = 4 ) points(npts2)
  real ( kind = 4 ) result
  real ( kind = 4 ) true
!
!  Singularity points:
!
  points(1) = 1.0E+00
  points(2) = sqrt ( 2.0E+00 )

  call qagp ( f04, a, b, npts2, points, epsabs, epsrel, result, abserr, &
    neval, ier )

  true = 61.0E+00 * log ( 2.0E+00 ) &
    + 77.0E+00 * log ( 7.0E+00 ) / 4.0E+00 - 27.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAGP_TEST'
  write ( *, '(a)' ) '  Test QAGP'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is x**3 * log(abs((x*x-1)*(x*x-2)))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qags_test ( )

!*****************************************************************************80
!
!! QAGS_TEST tests QAGS.
!
!  Discussion:
!
!    QAGS is an adaptive integrator for endpoint singularities.
!
!    integrate log(x)/sqrt(x) from 0 to 1.
!
!    The exact answer is -4.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f03
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) neval
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = -4.0E+00

  call qags ( f03, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAGS_TEST'
  write ( *, '(a)' ) '  Test QAGS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is LOG(X)/SQRT(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qawc_test ( )

!*****************************************************************************80
!
!! QAWC_TEST tests QAWC.
!
!  Discussion:
!
!    QAWC is an adaptive integrator for finding the Cauchy
!    principal value of the integral of f(x)*w(x) over (a,b)
!    where w(x)=1/(x-c), c between a and b.
!
!    Integrate 1/(x*(5*x*x*x+6)) from -1 to 5
!
!    The exact answer is log(125/631) / 18 = -0.08994401
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = -1.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 5.0E+00
  real ( kind = 4 ), parameter :: c = 0.0E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f09
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) neval
  real ( kind = 4 ) result
  real ( kind = 4 ) true

  call qawc ( f09, a, b, c, epsabs, epsrel, result, abserr, neval, ier )

  true = log ( 125.0E+00 / 631.0E+00 ) / 18.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAWC_TEST'
  write ( *, '(a)' ) '  Test QAWC'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is 1/(x*(5*x**3+6)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Point of singularity c =      ', c
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qawf_test ( )

!*****************************************************************************80
!
!! QAWF_TEST tests QAWF.
!
!  Discussion:
!
!    QAWF handles fourier integration of f(x)*w(x) from
!    a to infinity, with w(x)=cos(omega*x) or sine(omega*x)
!
!    integrate cos(pi*x/2) /sqrt(x) from 0 to infinity.
!
!    The exact answer is 1.0
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: epsabs = 0.001E+00
  real ( kind = 4 ), external :: f07
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: integr = 1
  integer ( kind = 4 ) neval
  real ( kind = 4 ) omega
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = 1.0E+00
!
!  set argument of sine or cosine
!  set integr=1 for cosine, 2 for sine
!
  omega = 0.5E+00 * r4_pi

  call qawf ( f07, a, omega, integr, epsabs, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAWF_TEST'
  write ( *, '(a)' ) '  Test QAWF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is cos(pi*x/2)/sqrt(x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qawo_test ( )

!*****************************************************************************80
!
!! QAWO_TEST tests QAWO.
!
!  Discussion:
!
!    QAWO integrates functions of the form 
!      f(x) * sin(omega*x)
!    or 
!      f(x) * cos(omega*x)
!
!    Here, we estimate 
!
!      Integral ( 0 <= x <= 1 ) log(x) sin(10*pi*x) dx
!
!    The exact answer is
!
!      exact = - ( gamma + log(10*pi) - ci(10*pi) ) / (10*pi)
!            = - 0.1281316...
!
!    Here Gamma is Euler's constant.
!
!    ci is the cosine integral:
!
!      ci(x) = integral ( x <= v < +oo ) - cos ( v ) / v dv.
!
!    We specify 
!      * INTEGR=1 for integrands with a cosine factor;
!      * INTEGR=2 for integrands with a sine factor.
!
!    Thanks to William Gandler for pointing out errors in the documentation
!    and text of this example, 29 October 2010.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), parameter :: ci = -0.001007E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f06
  real ( kind = 4 ), parameter :: gamma = 0.5772156649E+00
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: integr = 2
  integer ( kind = 4 ) neval
  real ( kind = 4 ) omega
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) result
  real ( kind = 4 ) true

  omega = 10.0E+00 * r4_pi

  call qawo ( f06, a, b, omega, integr, epsabs, epsrel, result, abserr, &
    neval, ier )

  true = - ( gamma + log ( 10.0E+00 * r4_pi ) - ci ) / ( 10.0E+00 * r4_pi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAWO_TEST'
  write ( *, '(a)' ) '  Test QAWO'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)*sin(10*pi*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qaws_test ( )

!*****************************************************************************80
!
!! QAWS_TEST tests QAWS.
!
!  Discussion:
!
!    QAWS is an adaptive integrator for integrands with
!    algebraic or logarithmic singularities at the endpoints.
!
!    Here we estimate:
!
!      Integral ( 0 <= x <= 1 ) log(x) / ( 1 + log(x)^2 )^2 dx
!
!    The exact answer is 
!
!      exact = 0.5 * ( ci(1.0) * sin(1.0) - si(G&R)(1.0) * cos(1.0) - 1.0 )
!
!    Numerically:
!
!      exact = -0.18927518788209332118...
!
!    ci is the cosine integral:
!
!      ci(x) = integral ( x <= v < oo ) - cos(v) / v dv
!
!    si (according to Mathematica, for instance) is the sine integral:
!
!      si(x) = integral ( 0 <= v <= x )    sin(v) / v dv
!
!    Note that when Gradshteyn and Rizhik refer to si(x), they
!    mean, instead the complementary form:
!
!      si(G&R)(x) = integral ( x <= v < oo ) sin(v) / v dv
!                 = ( pi / 2 ) - si(x)
!
!    ci(1.0)      = 0.33740392290096813466
!    si(1.0)      = 0.94608307036718301494
!    si(G&R)(1.0) = 0.62471325642771360429
!
!    Thanks to William Gandler for questioning a previous version of
!    the documentation and text of this example, which led to the 
!    clarification of the difference between the Gradshteyn & Rizhik
!    convention versus the Mathematica convention for the sine integral
!    function, 02 November 2010.
!
!    Note that the original QUADPACK documentation lists the answer as
!      (Ci(1)sin(1)+(pi/2-Si(1))*cos(1))/pi
!    which has an incorrect final divisor of pi (it should be 2) and
!    which uses the Mathematica convention for Si.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: alfa = 0.0E+00
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), parameter :: beta = 0.0E+00
  real ( kind = 4 ), parameter :: ci = 0.33740392290096813466E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f08
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) integr
  integer ( kind = 4 ) neval
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: si = 0.94608307036718301494E+00
  real ( kind = 4 ) true
!
!  INTEGR = 2 means the weight function is:
!
!    (x-a)**alfa * (b-x)**beta * log ( x - a )
!
  integr = 2

  call qaws ( f08, a, b, alfa, beta, integr, epsabs, epsrel, result, &
    abserr, neval, ier )

  true = 0.5E+00 * ( ci * sin ( 1.0E+00 ) &
    + ( r4_pi/ 2.0E+00 - si ) * cos ( 1.0E+00 ) - 1.0E+00 ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAWS_TEST'
  write ( *, '(a)' ) '  Test QAWS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+(log(x))**2)**2'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine qk15_test ( )

!*****************************************************************************80
!
!! QK15_TEST tests QK15.
!
!  Discussion:
!
!    QK15 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk15 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK15_TEST'
  write ( *, '(a)' ) '  Test QK15'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qk21_test ( )

!*****************************************************************************80
!
!! QK21_TEST tests QK21.
!
!  Discussion:
!
!    QK21 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk21 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK21_TEST'
  write ( *, '(a)' ) '  Test QK21'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qk31_test ( )

!*****************************************************************************80
!
!! QK31_TEST tests QK31.
!
!  Discussion:
!
!    QK31 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk31 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK31_TEST'
  write ( *, '(a)' ) '  Test QK31'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qk41_test ( )

!*****************************************************************************80
!
!! QK41_TEST tests QK41.
!
!  Discussion:
!
!    QK41 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk41 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK41_TEST'
  write ( *, '(a)' ) '  Test QK41'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qk51_test ( )

!*****************************************************************************80
!
!! QK51_TEST tests QK51.
!
!  Discussion:
!
!    QK51 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk51 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK51_TEST'
  write ( *, '(a)' ) '  Test QK51'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qk61_test ( )

!*****************************************************************************80
!
!! QK61_TEST tests QK61.
!
!  Discussion:
!
!    QK61 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), external :: f01
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qk61 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QK61_TEST'
  write ( *, '(a)' ) '  Test QK61'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine qng_test ( )

!*****************************************************************************80
!
!! QNG_TEST tests QNG.
!
!  Discussion:
!
!    QNG is a nonadaptive automatic integrator using a Gauss-Kronrod or 
!    Patterson rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), parameter :: a = 0.0E+00
  real ( kind = 4 ) abserr
  real ( kind = 4 ), parameter :: b = 1.0E+00
  real ( kind = 4 ), parameter :: epsabs = 0.0E+00
  real ( kind = 4 ), parameter :: epsrel = 0.001E+00
  real ( kind = 4 ), external :: f01
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) neval
  real ( kind = 4 ) result
  real ( kind = 4 ), parameter :: true = - 4.0E+00 / 9.0E+00

  call qng ( f01, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QNG_TEST'
  write ( *, '(a)' ) '  Test QNG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
function f01 ( x )

!*****************************************************************************80
!
!! F01 is the integrand function SQRT(X) * LOG(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f01
  real ( kind = 4 ) x

  if ( x <= 0.0E+00 )then
    f01 = 0.0E+00
  else
    f01 = sqrt ( x ) * log ( x )
  end if

  return
end
function f02 ( x )

!*****************************************************************************80
!
!! F02 is the integrand function COS(100*SIN(X)).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f02
  real ( kind = 4 ) x

  f02 = cos ( 100.0E+00 * sin ( x ) )

  return
end
function f03 ( x )

!*****************************************************************************80
!
!! F03 is the integrand function LOG(X)/SQRT(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f03
  real ( kind = 4 ) x

  if ( x <= 0.0E+00 ) then
    f03 = 0.0E+00
  else
    f03 = log ( x ) / sqrt ( x )
  end if

  return
end
function f04 ( x )

!*****************************************************************************80
!
!! F04 is the integrand function X^3 LOG((X^2-1)*(X^2-2))
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f04
  real ( kind = 4 ) x

  f04 = x**3 * log ( abs ( ( x**2 - 1.0E+00 ) * ( x**2 - 2.0E+00 ) ) )

  return
end
function f05 ( x )

!*****************************************************************************80
!
!! F05 is the integrand function LOG(X)/(1+100X^2).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f05
  real ( kind = 4 ) x

  f05 = log ( x ) / ( 1.0E+00 + 100.0E+00 * x**2 )

  return
end
function f06 ( x )

!*****************************************************************************80
!
!! F06 is the integrand function LOG(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f06
  real ( kind = 4 ) x

  if ( x <= 0.0E+00 ) then
    f06 = 0.0E+00
  else
    f06 = log ( x )
  end if

  return
end
function f07 ( x )

!*****************************************************************************80
!
!! F07 is the integrand function 1/SQRT(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f07
  real ( kind = 4 ) x

  if ( x <= 0.0E+00 ) then
    f07 = 0.0E+00
  else
    f07 = 1.0E+00 / sqrt ( x )
  end if

  return
end
function f08 ( x )

!*****************************************************************************80
!
!! F08 is the integrand function 1/(1+LOG(X)**2)**2
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f08
  real ( kind = 4 ) x

  if ( 0.0E+00 < x ) then
    f08 = 1.0E+00 / ( 1.0E+00 + log ( x )**2 )**2
  else
    f08 = 0.0E+00
  end if

  return
end
function f09 ( x )

!*****************************************************************************80
!
!! F09 is the integrand function 1 / ( 5 X^3 + 6 ).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) f09
  real ( kind = 4 ) x

  f09 = 1.0E+00 / ( 5.0E+00 * x**3 + 6.0E+00 )

  return
end
