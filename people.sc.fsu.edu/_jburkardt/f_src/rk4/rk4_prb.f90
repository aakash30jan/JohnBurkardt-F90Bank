program main

!*****************************************************************************80
!
!! MAIN is the main program for RK4_PRB.
!
!  Discussion:
!
!    RK4_PRB tests the RK4 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the RK4 library.'

  call rk4_test ( )
  call rk4vec_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine rk4_test ( )

!*****************************************************************************80
!
!! RK4_TEST tests the RK4 routine for a scalar ODE.
!
!  Discussion:
!
!    RK4_PRB tests the RK4 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, real ( kind = 8 ) DT, the time step.
!
!    Local, real ( kind = 8 ) T0, the time at which the solution is known.
!
!    Local, real ( kind = 8 ) TMAX, the maximum time at which a solution 
!    is desired.
!
!    Local, real ( kind = 8 ) U0, the estimated solution at time T0.
!
  implicit none

  real ( kind = 8 ), parameter :: dt = 0.1D+00
  external rk4_test_f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ), parameter :: tmax = 12.0D+00 * pi
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4_TEST'
  write ( *, '(a)' ) '  RK4 takes a Runge Kutta step for a scalar ODE.'
  write ( *, '(a)' ) ' '

  t0 = 0.0D+00
  u0 = 0.5D+00

  do
!
!  Print (T0,U0).
!
    write ( *, '(2x,g14.6,2x,g14.6)' ) t0, u0
!
!  Stop if we've exceeded TMAX.
!
    if ( tmax <= t0 ) then
      exit
    end if
!
!  Otherwise, advance to time T1, and have RK4 estimate 
!  the solution U1 there.
!
    t1 = t0 + dt
    call rk4 ( t0, u0, dt, rk4_test_f, u1 )
!
!  Shift the data to prepare for another step.
!
    t0 = t1
    u0 = u1

  end do

  return
end
subroutine rk4_test_f ( t, u, uprime )

!*****************************************************************************80
!
!! RK4_TEST_F evaluates the right hand side of a particular ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) U, the current solution value.
!
!    Output, real ( kind = 8 ) UPRIME, the value of the derivative, dU/dT.
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) uprime
  
  uprime = u * cos ( t )
 
  return
end
subroutine rk4vec_test ( )

!*****************************************************************************80
!
!! RK4VEC_TEST tests RK4VEC for a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), parameter :: dt = 0.1D+00
  external rk4vec_test_f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ), parameter :: tmax = 12.0D+00 * pi
  real ( kind = 8 ) u0(n)
  real ( kind = 8 ) u1(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4VEC_TEST'
  write ( *, '(a)' ) '  RK4VEC takes a Runge Kutta step for a vector ODE.'
  write ( *, '(a)' ) ' '

  t0 = 0.0D+00
  u0(1) = 0.0D+00
  u0(2) = 1.0D+00

  do
!
!  Print (T0,U0).
!
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) t0, u0(1), u0(2)
!
!  Stop if we've exceeded TMAX.
!
    if ( tmax <= t0 ) then
      exit
    end if
!
!  Otherwise, advance to time T1, and have RK4 estimate 
!  the solution U1 there.
!
    t1 = t0 + dt
    call rk4vec ( t0, n, u0, dt, rk4vec_test_f, u1 )
!
!  Shift the data to prepare for another step.
!
    t0 = t1
    u0(1:n) = u1(1:n)

  end do

  return
end
subroutine rk4vec_test_f ( t, n, u, uprime )

!*****************************************************************************80
!
!! RK4VEC_TEST_F evaluates the right hand side of a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, integer ( kind = 4 ) N, the dimension of the system.
!
!    Input, real ( kind = 8 ) U(N), the current solution value.
!
!    Output, real ( kind = 8 ) UPRIME(N), the value of the derivative, dU/dT.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) t
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uprime(n)
  
  uprime(1) = u(2)
  uprime(2) = - u(1)
 
  return
end
