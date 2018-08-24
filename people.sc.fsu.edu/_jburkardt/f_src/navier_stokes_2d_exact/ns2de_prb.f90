program main

!*****************************************************************************80
!
!! NS2DE_TEST tests the NS2DE library.
!
!  Location:
!
!    http://people.sc.fsu.edu/~jburkardt/f_src/navier_stokes_2d_exact/ns2de_prb.f90
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NS2DE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NS2DE library.'
!
!  Lucas flow.
!
  call uvp_lucas_test ( )
  call uvp_lucas_test2 ( )
  call rhs_lucas_test ( )
  call resid_lucas_test ( )
  call gnuplot_lucas_test ( )
!
!  Poiseuille flow.
!
  call uvp_poiseuille_test ( )
  call uvp_poiseuille_test2 ( )
  call rhs_poiseuille_test ( )
  call resid_poiseuille_test ( )
  call gnuplot_poiseuille_test ( )
  call parameter_poiseuille_test ( )
!
!  Spiral flow.
!
  call uvp_spiral_test ( )
  call uvp_spiral_test2 ( )
  call rhs_spiral_test ( )
  call resid_spiral_test ( )
  call gnuplot_spiral_test ( )
  call parameter_spiral_test ( )
!
!  Taylor flow.
!
  call uvp_taylor_test ( )
  call uvp_taylor_test2 ( )
  call rhs_taylor_test ( )
  call resid_taylor_test ( )
  call gnuplot_taylor_test ( )
  call parameter_taylor_test ( )
!
!  Vortex flow.
!
  call uvp_vortex_test ( )
  call uvp_vortex_test2 ( )
  call rhs_vortex_test ( )
  call resid_vortex_test ( )
  call gnuplot_vortex_test ( )
  call parameter_vortex_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NS2DE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  return
end
subroutine uvp_lucas_test ( )

!*****************************************************************************80
!
!! UVP_LUCAS_TEST samples Lucas Bystricky's solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_LUCAS_TEST'
  write ( *, '(a)' ) '  Lucas Bystricky''s flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = +0.0D+00
  r8_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call uvp_lucas ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_lucas_test2 ( )

!*****************************************************************************80
!
!! UVP_LUCAS_TEST2 samples Lucas Bystricky's solution on the bounadry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  r8_lo = 0.0D+00
  r8_hi = 1.0D+00

  nu = 1.0D+00
  rho = 1.0D+00
  t = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_LUCAS_TEST2'
  write ( *, '(a)' ) '  Lucas Bystricky''s flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  along the boundary,'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 400
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call r8vec_linspace ( 100, r8_lo, r8_hi, x(1:100) )
  y(1:100) = r8_lo

  x(101:200) = r8_hi
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(101:200) )

  call r8vec_linspace ( 100, r8_hi, r8_lo, x(201:300) )
  y(201:300) = r8_hi

  x(301:400) = r8_lo
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(301:400) )

  call uvp_lucas ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine rhs_lucas_test ( )

!*****************************************************************************80
!
!! RHS_LUCAS_TEST samples Lucas Bystricky's right hand side at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ), allocatable :: g(:)
  real ( kind = 8 ), allocatable :: h(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RHS_LUCAS_TEST'
  write ( *, '(a)' ) '  Lucas Bystricky''s flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes right hand sides'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( h(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = +0.0D+00
  r8_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call rhs_lucas ( nu, rho, n, x, y, t, f, g, h )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Uf:  ', minval ( f ), maxval ( f )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vf:  ', minval ( g ), maxval ( g )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pf:  ', minval ( h ), maxval ( h )

  deallocate ( f )
  deallocate ( g )
  deallocate ( h )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine resid_lucas_test ( )

!*****************************************************************************80
!
!! RESID_LUCAS_TEST samples Lucas Bystricky's residual at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: ur(:)
  real ( kind = 8 ), allocatable :: vr(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RESID_LUCAS_TEST'
  write ( *, '(a)' ) '  Lucas Bystricky''s flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes residuals'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( pr(1:n) )
  allocate ( ur(1:n) )
  allocate ( vr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = +0.0D+00
  r8_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call resid_lucas ( nu, rho, n, x, y, t, ur, vr, pr )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Ur:  ', minval ( abs ( ur ) ), maxval ( abs ( ur ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Vr:  ', minval ( abs ( vr ) ), maxval ( abs ( vr ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( ur )
  deallocate ( vr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine gnuplot_lucas_test ( )

!*****************************************************************************80
!
!! GNUPLOT_LUCAS_TEST plots a Lucas Bystricky velocity field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) header
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable :: y(:,:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo
  integer ( kind = 4 ) y_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GNUPLOT_LUCAS_TEST:'
  write ( *, '(a)' ) '  Lucas Bystricky flow'
  write ( *, '(a)' ) '  Generate a velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.0D+00
  x_hi = 1.0D+00
  x_num = 21

  y_lo = 0.0D+00
  y_hi = 1.0D+00
  y_num = 21

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  nu = 1.0D+00
  rho = 1.0D+00
  n = x_num * y_num
  t = 0.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )
  allocate ( p(1:x_num,1:y_num) )

  call uvp_lucas ( nu, rho, n, x, y, t, u, v, p )

  header = 'lucas'
  s = 0.25D+00
  call ns2de_gnuplot ( header, n, x, y, u, v, p, s )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_poiseuille_test ( )

!*****************************************************************************80
!
!! UVP_POISEUILLE_TEST samples Poiseuille solution at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo

  nu = 1.0D+00
  rho = 1.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_POISEUILLE_TEST'
  write ( *, '(a)' ) '  Poiseuille flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  at the initial time T = 0, using a channel region.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  x_lo = +0.0D+00
  x_hi = +6.0D+00
  y_lo = -1.0D+00
  y_hi = +1.0D+00

  seed = 123456789

  call r8vec_uniform_ab ( n, x_lo, x_hi, seed, x )
  call r8vec_uniform_ab ( n, y_lo, y_hi, seed, y )
  t = 0.0D+00

  call uvp_poiseuille ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_poiseuille_test2 ( )

!*****************************************************************************80
!
!! UVP_POISEUILLE_TEST2 samples Poiseuille flow on the bounadry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo

  nu = 1.0D+00
  rho = 1.0D+00
  t = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_POISEUILLE_TEST2'
  write ( *, '(a)' ) '  Poiseuille flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  along the boundary,'
  write ( *, '(a)' ) '  at the initial time T = 0, using a channel region.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 400
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  x_lo = 0.0D+00
  x_hi = 6.0D+00
  y_lo = -1.0D+00
  y_hi = +1.0D+00

  call r8vec_linspace ( 100, x_lo, x_hi, x(1:100) )
  y(1:100) = y_lo

  x(101:200) = x_hi
  call r8vec_linspace ( 100, y_lo, y_hi, y(101:200) )

  call r8vec_linspace ( 100, x_hi, x_lo, x(201:300) )
  y(201:300) = y_hi

  x(301:400) = x_lo
  call r8vec_linspace ( 100, y_hi, y_lo, y(301:400) )

  call uvp_poiseuille ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine rhs_poiseuille_test ( )

!*****************************************************************************80
!
!! RHS_POISEUILLE_TEST samples Poiseuille flow right hand side at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ), allocatable :: g(:)
  real ( kind = 8 ), allocatable :: h(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RHS_POISEUILLE_TEST'
  write ( *, '(a)' ) '  Poiseuille flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes right hand sides'
  write ( *, '(a)' ) '  at the initial time T = 0, using a channel region.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( h(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  x_lo = +0.0D+00
  x_hi = +6.0D+00
  y_lo = -1.0D+00
  y_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, x_lo, x_hi, seed, x )

  call r8vec_uniform_ab ( n, y_lo, y_hi, seed, y )

  t = 0.0D+00

  call rhs_poiseuille ( nu, rho, n, x, y, t, f, g, h )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Uf:  ', minval ( f ), maxval ( f )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vf:  ', minval ( g ), maxval ( g )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pf:  ', minval ( h ), maxval ( h )

  deallocate ( f )
  deallocate ( g )
  deallocate ( h )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine resid_poiseuille_test ( )

!*****************************************************************************80
!
!! RESID_POISEUILLE_TEST samples Poiseuille flow residuals at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: ur(:)
  real ( kind = 8 ), allocatable :: vr(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RESID_POISEUILLE_TEST'
  write ( *, '(a)' ) '  Poiseuille flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes residuals'
  write ( *, '(a)' ) '  at the initial time T = 0, using a channel region.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( pr(1:n) )
  allocate ( ur(1:n) )
  allocate ( vr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  x_lo = +0.0D+00
  x_hi = +6.0D+00
  y_lo = -1.0D+00
  y_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, x_lo, x_hi, seed, x )

  call r8vec_uniform_ab ( n, y_lo, y_hi, seed, y )

  t = 0.0D+00

  call resid_poiseuille ( nu, rho, n, x, y, t, ur, vr, pr )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Ur:  ', minval ( abs ( ur ) ), maxval ( abs ( ur ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Vr:  ', minval ( abs ( vr ) ), maxval ( abs ( vr ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( ur )
  deallocate ( vr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine gnuplot_poiseuille_test ( )

!*****************************************************************************80
!
!! GNUPLOT_POISEUILLE_TEST plots a Poiseuille velocity field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) header
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable :: y(:,:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo
  integer ( kind = 4 ) y_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GNUPLOT_POISEUILLE_TEST:'
  write ( *, '(a)' ) '  Poiseuille flow'
  write ( *, '(a)' ) '  Generate a velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.0D+00
  x_hi = 6.0D+00
  x_num = 61

  y_lo = -1.0D+00
  y_hi =  1.0D+00
  y_num = 21

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  nu = 1.0D+00
  rho = 1.0D+00
  n = x_num * y_num
  t = 0.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )
  allocate ( p(1:x_num,1:y_num) )
  call uvp_poiseuille ( nu, rho, n, x, y, t, u, v, p )

  header = 'poiseuille'
  s = 0.25D+00
  call ns2de_gnuplot ( header, n, x, y, u, v, p, s )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine parameter_poiseuille_test ( )

!*****************************************************************************80
!
!! PARAMETER_POISEUILLE_TEST monitors Poiseuille solution norms for various NU, RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) p_norm
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ) u_norm
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) v_norm
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARAMETER_POISEUILLE_TEST'
  write ( *, '(a)' ) '  Poiseuille Flow'
  write ( *, '(a)' ) '  Monitor solution norms over time for various'
  write ( *, '(a)' ) '  values of NU, RHO.'

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  x_lo =  0.0D+00
  x_hi =  6.0D+00
  y_lo = -1.0D+00
  y_hi =  1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, x_lo, x_hi, seed, x )
  call r8vec_uniform_ab ( n, y_lo, y_hi, seed, y )
!
!  Vary RHO.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Let RHO vary:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00

  do j = 1, 3

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_poiseuille ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''
    rho = rho / 100.0D+00

  end do
!
!  Vary NU.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Let NU vary:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00
  
  do i = 1, 4

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_poiseuille ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''

    nu = nu / 10.0D+00

  end do

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_spiral_test ( )

!*****************************************************************************80
!
!! UVP_SPIRAL_TEST samples the Spiral solution at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_SPIRAL_TEST'
  write ( *, '(a)' ) '  Spiral flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  at the initial time T = 0, over the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.0D+00
  r8_hi = +1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call uvp_spiral ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_spiral_test2 ( )

!*****************************************************************************80
!
!! UVP_SPIRAL_TEST2 samples the Spiral solution at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  r8_lo = 0.0D+00
  r8_hi = 1.0D+00

  nu = 1.0D+00
  rho = 1.0D+00
  t = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_SPIRAL_TEST2'
  write ( *, '(a)' ) '  Spiral flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  along the boundary,'
  write ( *, '(a)' ) '  at the initial time T = 0, over the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 400
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call r8vec_linspace ( 100, r8_lo, r8_hi, x(1:100) )
  y(1:100) = r8_lo

  x(101:200) = r8_hi
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(101:200) )

  call r8vec_linspace ( 100, r8_hi, r8_lo, x(201:300) )
  y(201:300) = r8_hi

  x(301:400) = r8_lo
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(301:400) )

  call uvp_spiral ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine rhs_spiral_test ( )

!*****************************************************************************80
!
!! RHS_SPIRAL_TEST samples the Spiral right hand side at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ), allocatable :: g(:)
  real ( kind = 8 ), allocatable :: h(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RHS_SPIRAL_TEST'
  write ( *, '(a)' ) '  Spiral Flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes right hand sides'
  write ( *, '(a)' ) '  at the initial time T = 0, over the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( h(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.0D+00
  r8_hi = 1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call rhs_spiral ( nu, rho, n, x, y, t, f, g, h )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Uf:  ', minval ( f ), maxval ( f )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vf:  ', minval ( g ), maxval ( g )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pf:  ', minval ( h ), maxval ( h )

  deallocate ( f )
  deallocate ( g )
  deallocate ( h )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine resid_spiral_test ( )

!*****************************************************************************80
!
!! RESID_SPIRAL_TEST samples the Spiral residual at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: ur(:)
  real ( kind = 8 ), allocatable :: vr(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RESID_SPIRAL_TEST'
  write ( *, '(a)' ) '  Spiral Flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes residuals'
  write ( *, '(a)' ) '  at the initial time T = 0, over the unit square.'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( pr(1:n) )
  allocate ( ur(1:n) )
  allocate ( vr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.0D+00
  r8_hi = 1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call resid_spiral ( nu, rho, n, x, y, t, ur, vr, pr )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Ur:  ', minval ( abs ( ur ) ), maxval ( abs ( ur ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Vr:  ', minval ( abs ( vr ) ), maxval ( abs ( vr ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) &
    '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( ur )
  deallocate ( vr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine gnuplot_spiral_test ( )

!*****************************************************************************80
!
!! GNUPLOT_SPIRAL_TEST plots a Spiral velocity field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) header
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable :: y(:,:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo
  integer ( kind = 4 ) y_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GNUPLOT_SPIRAL_TEST:'
  write ( *, '(a)' ) '  Spiral Flow'
  write ( *, '(a)' ) '  Generate a velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.0D+00
  x_hi = 1.0D+00
  x_num = 21

  y_lo = 0.0D+00
  y_hi = 1.0D+00
  y_num = 21

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  nu = 1.0D+00
  rho = 1.0D+00
  n = x_num * y_num
  t = 0.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )
  allocate ( p(1:x_num,1:y_num) )
  call uvp_spiral ( nu, rho, n, x, y, t, u, v, p )

  header = 'spiral'
  s = 5.0D+00
  call ns2de_gnuplot ( header, n, x, y, u, v, p, s )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine parameter_spiral_test ( )

!*****************************************************************************80
!
!! PARAMETER_SPIRAL_TEST monitors Spiral solution norms for various NU, RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) p_norm
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ) u_norm
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) v_norm
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARAMETER_SPIRAL_TEST'
  write ( *, '(a)' ) '  Spiral Flow'
  write ( *, '(a)' ) '  Monitor solution norms over time for various'
  write ( *, '(a)' ) '  values of NU, RHO.'

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.0D+00
  r8_hi = 1.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
!
!  Vary RHO.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Let RHO vary:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00

  do j = 1, 3

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_spiral ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''
    rho = rho / 100.0D+00

  end do
!
!  Vary NU.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Let NU vary:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00
  
  do i = 1, 4

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_spiral ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''

    nu = nu / 10.0D+00

  end do

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_taylor_test ( )

!*****************************************************************************80
!
!! UVP_TAYLOR_TEST samples the Taylor solution at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_TAYLOR_TEST'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +2.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call uvp_taylor ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_taylor_test2 ( )

!*****************************************************************************80
!
!! UVP_TAYLOR_TEST2 samples the Taylor solution on the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  r8_lo = 0.5D+00
  r8_hi = 2.5D+00

  nu = 1.0D+00
  rho = 1.0D+00
  t = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_TAYLOR_TEST2'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  along the boundary,'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 400
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call r8vec_linspace ( 100, r8_lo, r8_hi, x(1:100) )
  y(1:100) = r8_lo

  x(101:200) = r8_hi
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(101:200) )

  call r8vec_linspace ( 100, r8_hi, r8_lo, x(201:300) )
  y(201:300) = r8_hi

  x(301:400) = r8_lo
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(301:400) )

  call uvp_taylor ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine rhs_taylor_test ( )

!*****************************************************************************80
!
!! RHS_TAYLOR_TEST samples the Taylor right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ), allocatable :: g(:)
  real ( kind = 8 ), allocatable :: h(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RHS_TAYLOR_TEST'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes right hand sides'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( h(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +2.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call rhs_taylor ( nu, rho, n, x, y, t, f, g, h )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Uf:  ', minval ( f ), maxval ( f )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vf:  ', minval ( g ), maxval ( g )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pf:  ', minval ( h ), maxval ( h )

  deallocate ( f )
  deallocate ( g )
  deallocate ( h )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine resid_taylor_test ( )

!*****************************************************************************80
!
!! RESID_TAYLOR_TEST samples the Taylor residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: ur(:)
  real ( kind = 8 ), allocatable :: vr(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RESID_TAYLOR_TEST'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes residuals'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( pr(1:n) )
  allocate ( ur(1:n) )
  allocate ( vr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +2.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call resid_taylor ( nu, rho, n, x, y, t, ur, vr, pr )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Ur:  ', minval ( abs ( ur ) ), maxval ( abs ( ur ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vr:  ', minval ( abs ( vr ) ), maxval ( abs ( vr ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( ur )
  deallocate ( vr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine gnuplot_taylor_test ( )

!*****************************************************************************80
!
!! GNUPLOT_TAYLOR_TEST plots a Taylor velocity field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) header
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable :: y(:,:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo
  integer ( kind = 4 ) y_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GNUPLOT_TAYLOR_TEST:'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Generate a velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.5D+00
  x_hi = 2.5D+00
  x_num = 21

  y_lo = 0.5D+00
  y_hi = 2.5D+00
  y_num = 21

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  nu = 1.0D+00
  rho = 1.0D+00
  n = x_num * y_num
  t = 0.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )
  allocate ( p(1:x_num,1:y_num) )
  call uvp_taylor ( nu, rho, n, x, y, t, u, v, p )

  header = 'taylor'
  s = 0.10D+00
  call ns2de_gnuplot ( header, n, x, y, u, v, p, s )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine parameter_taylor_test ( )

!*****************************************************************************80
!
!! PARAMETER_TAYLOR_TEST monitors Taylor solution norms for various NU, RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) p_norm
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ) u_norm
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) v_norm
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARAMETER_TAYLOR_TEST'
  write ( *, '(a)' ) '  Taylor flow'
  write ( *, '(a)' ) '  Monitor solution norms over time for various'
  write ( *, '(a)' ) '  values of NU, RHO.'

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +2.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
!
!  Vary RHO.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  RHO affects the pressure scaling.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00

  do j = 1, 3

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_taylor ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''
    rho = rho / 100.0D+00

  end do
!
!  Vary NU.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  NU affects the time scaling.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00
  
  do i = 1, 4

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_taylor ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''

    nu = nu / 10.0D+00

  end do

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_vortex_test ( )

!*****************************************************************************80
!
!! UVP_VORTEX_TEST samples the Vortex solution at the initial time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_VORTEX_TEST'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +1.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call uvp_vortex ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine uvp_vortex_test2 ( )

!*****************************************************************************80
!
!! UVP_VORTEX_TEST2 samples the Vortex solution on the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  r8_lo = 0.5D+00
  r8_hi = 1.5D+00

  nu = 1.0D+00
  rho = 1.0D+00
  t = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UVP_VORTEX_TEST2'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Estimate the range of velocity and pressure'
  write ( *, '(a)' ) '  along the boundary,'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 400
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call r8vec_linspace ( 100, r8_lo, r8_hi, x(1:100) )
  y(1:100) = r8_lo

  x(101:200) = r8_hi
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(101:200) )

  call r8vec_linspace ( 100, r8_hi, r8_lo, x(201:300) )
  y(201:300) = r8_hi

  x(301:400) = r8_lo
  call r8vec_linspace ( 100, r8_lo, r8_hi, y(301:400) )

  call uvp_vortex ( nu, rho, n, x, y, t, u, v, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  U:  ', minval ( u ), maxval ( u )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  V:  ', minval ( v ), maxval ( v )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  P:  ', minval ( p ), maxval ( p )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine rhs_vortex_test ( )

!*****************************************************************************80
!
!! RHS_VORTEX_TEST samples the Vortex right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ), allocatable :: g(:)
  real ( kind = 8 ), allocatable :: h(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RHS_VORTEX_TEST'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes right hand sides'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( h(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +1.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call rhs_vortex ( nu, rho, n, x, y, t, f, g, h )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Uf:  ', minval ( f ), maxval ( f )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vf:  ', minval ( g ), maxval ( g )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pf:  ', minval ( h ), maxval ( h )

  deallocate ( f )
  deallocate ( g )
  deallocate ( h )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine resid_vortex_test ( )

!*****************************************************************************80
!
!! RESID_VORTEX_TEST samples the Vortex residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: ur(:)
  real ( kind = 8 ), allocatable :: vr(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  nu = 1.0D+00
  rho = 1.0D00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RESID_VORTEX_TEST'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Sample the Navier-Stokes residuals'
  write ( *, '(a)' ) '  at the initial time T = 0, using a region that is'
  write ( *, '(a)' ) '  the square centered at (1.5,1.5) with "radius" 1.0,'
  write ( *, '(a,g14.6)' ) '  Kinematic viscosity NU = ', nu
  write ( *, '(a,g14.6)' ) '  Fluid density RHO = ', rho

  n = 1000
  allocate ( pr(1:n) )
  allocate ( ur(1:n) )
  allocate ( vr(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +1.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  t = 0.0D+00

  call resid_vortex ( nu, rho, n, x, y, t, ur, vr, pr )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '           Minimum       Maximum'
  write ( *, '(a)' ) ''
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Ur:  ', minval ( abs ( ur ) ), maxval ( abs ( ur ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Vr:  ', minval ( abs ( vr ) ), maxval ( abs ( vr ) )
  write ( *, '(a,2x,g14.6,2x,g14.6)' ) '  Pr:  ', minval ( abs ( pr ) ), maxval ( abs ( pr ) )

  deallocate ( pr )
  deallocate ( ur )
  deallocate ( vr )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine gnuplot_vortex_test ( )

!*****************************************************************************80
!
!! GNUPLOT_VORTEX_TEST plots a Vortex velocity field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) header
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable :: y(:,:)
  real ( kind = 8 ) y_hi
  real ( kind = 8 ) y_lo
  integer ( kind = 4 ) y_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GNUPLOT_VORTEX_TEST:'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Generate a velocity field on a regular grid.'
  write ( *, '(a)' ) '  Store in GNUPLOT data and command files.'

  x_lo = 0.5D+00
  x_hi = 1.5D+00
  x_num = 21

  y_lo = 0.5D+00
  y_hi = 1.5D+00
  y_num = 21

  allocate ( x(1:x_num,1:y_num) )
  allocate ( y(1:x_num,1:y_num) )

  call grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

  nu = 1.0D+00
  rho = 1.0D+00
  n = x_num * y_num
  t = 0.0D+00

  allocate ( u(1:x_num,1:y_num) )
  allocate ( v(1:x_num,1:y_num) )
  allocate ( p(1:x_num,1:y_num) )
  call uvp_vortex ( nu, rho, n, x, y, t, u, v, p )

  header = 'vortex'
  s = 0.10D+00
  call ns2de_gnuplot ( header, n, x, y, u, v, p, s )

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine parameter_vortex_test ( )

!*****************************************************************************80
!
!! PARAMETER_VORTEX_TEST monitors Vortex solution norms for various NU, RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) p_norm
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ) u_norm
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) v_norm
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PARAMETER_VORTEX_TEST'
  write ( *, '(a)' ) '  Vortex flow'
  write ( *, '(a)' ) '  Monitor solution norms over time for various'
  write ( *, '(a)' ) '  values of NU, RHO.'

  n = 1000
  allocate ( p(1:n) )
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  r8_lo = 0.5D+00
  r8_hi = +2.5D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
!
!  Vary RHO.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  RHO affects the pressure scaling.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00

  do j = 1, 3

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_vortex ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''
    rho = rho / 100.0D+00

  end do
!
!  Vary NU.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  NU affects the time scaling.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     RHO         NU           T     ||U||       ||V||       ||P||'
  write ( *, '(a)' ) ''

  nu = 1.0D+00
  rho = 1.0D+00
  
  do i = 1, 4

    do k = 0, 5

      t = real ( k, kind = 8 ) / 5.0D+00

      call uvp_vortex ( nu, rho, n, x, y, t, u, v, p )

      u_norm = r8vec_norm_l2 ( n, u ) / real ( n, kind = 8 )
      v_norm = r8vec_norm_l2 ( n, v ) / real ( n, kind = 8 )
      p_norm = r8vec_norm_l2 ( n, p ) / real ( n, kind = 8 )

      write ( *, '(2x,g10.4,2x,g10.4,2x,f8.4,2x,g10.4,2x,g10.4,2x,g10.4)' ) &
        rho, nu, t, u_norm, v_norm, p_norm

    end do

    write ( *, '(a)' ) ''

    nu = nu / 10.0D+00

  end do

  deallocate ( p )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )

  return
end

