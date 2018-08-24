program main

!*****************************************************************************80
!
!! MAIN is the main program for R8PBU_PRB.
!
!  Discussion:
!
!    R8PBU_PRB tests the R8PBU library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8PBU library.'

  call r8pbu_cg_test ( )
  call r8pbu_det_test ( )
  call r8pbu_dif2_test ( )
  call r8pbu_fa_test ( )
  call r8pbu_indicator_test ( )
  call r8pbu_ml_test ( )
  call r8pbu_mv_test ( )
  call r8pbu_print_test ( )
  call r8pbu_print_some_test ( )
  call r8pbu_random_test ( )
  call r8pbu_res_test ( )
  call r8pbu_sl_test ( )
  call r8pbu_sor_test ( )
  call r8pbu_to_r8ge_test ( )
  call r8pbu_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8pbu_cg_test ( )

!*****************************************************************************80
!
!! R8PBU_CG_TEST tests R8PBU_CG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) err
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_CG_TEST'
  write ( *, '(a)' ) '  R8PBU_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '  to a symmetric positive definite banded '
  write ( *, '(a)' ) '  linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix values.
!
  a(2,1:n) = 2.0D+00
  a(1,2:n) = -1.0D+00

  call r8pbu_print_some ( n, mu, a, 1, 1, 10, 10, &
    'The symmetric banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  call r8pbu_mv ( n, n, mu, a, x, b )
!
!  Set the approximate solution.
!
  x(1:n) = 1.0D+00
!
!  Call the conjugate gradient method.
!
  call r8pbu_cg ( n, mu, a, b, x )
!
!  Compute the residual, A*x-b
!
  call r8pbu_mv ( n, n, mu, a, x, r )
 
  err = maxval ( abs ( r(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, 10, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual = ', err
 
  return
end
subroutine r8pbu_det_test ( )

!*****************************************************************************80
!
!! R8PBU_DET_TEST tests R8PBU_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_DET_TEST'
  write ( *, '(a)' ) '  R8PBU_DET, determinant of a positive definite'
  write ( *, '(a)' ) '  symmetric banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8pbu_random ( n, mu, seed, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8pbu_to_r8ge ( n, mu, a, a2 )
!
!  Factor the matrix.
!
  call r8pbu_fa ( n, mu, a, info )

  call r8pbu_print ( n, mu, a, '  The R8PBU factored matrix:' )
!
!  Compute the determinant.
!
  call r8pbu_det ( n, mu, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8PBU_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call r8ge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call r8ge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  R8GE_DET computes the determinant =  ', det

  return
end
subroutine r8pbu_dif2_test ( )

!*****************************************************************************80
!
!! R8PBU_DIF2_TEST tests R8PBU_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_DIF2_TEST'
  write ( *, '(a)' ) '  R8PBU_DIF2 sets up an R8PBU second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  m = n
  call r8pbu_dif2 ( m, n, mu, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU second difference matrix:' )

  return
end
subroutine r8pbu_fa_test ( )

!*****************************************************************************80
!
!! R8PBU_FA_TEST tests R8PBU_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_FA_TEST'
  write ( *, '(a)' ) '  R8PBU_FA factors a banded positive definite symmetric matrix,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  call r8pbu_random ( n, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  call r8pbu_mv ( n, n, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8pbu_fa ( n, mu, a, info )
!
!  Solve the linear system.
!
  call r8pbu_sl ( n, mu, a, b )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
 
  return
end
subroutine r8pbu_indicator_test ( )

!*****************************************************************************80
!
!! R8PBU_INDICATOR_TEST tests R8PBU_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8PBU_INDICATOR sets up an R8PBU indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call r8pbu_indicator ( n, mu, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU indicator matrix:' )

  return
end
subroutine r8pbu_ml_test ( )

!*****************************************************************************80
!
!! R8PBU_ML_TEST tests R8PBU_ML.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_ML_TEST'
  write ( *, '(a)' ) '  R8PBU_ML computes A*x '
  write ( *, '(a)' ) '  where A has been factored by R8PBU_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8pbu_random ( n, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8pbu_mv ( n, n, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8pbu_fa ( n, mu, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_ML_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8PBU_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
  call r8pbu_ml ( n, mu, a, x, b2 )

  call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )

  return
end
subroutine r8pbu_mv_test ( )

!*****************************************************************************80
!
!! R8PBU_MV_TEST tests R8PBU_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_MV_TEST'
  write ( *, '(a)' ) '  R8PBU_MV computes A*x '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8pbu_random ( n, mu, seed, a )
  call r8pbu_print ( n, mu, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )
!
!  Compute the corresponding right hand side.
!
  call r8pbu_mv ( n, n, mu, a, x, b )
  call r8vec_print ( n, b, '  Product b=A*x' )

  return
end
subroutine r8pbu_print_test ( )

!*****************************************************************************80
!
!! R8PBU_PRINT_TEST tests R8PBU_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_PRINT_TEST'
  write ( *, '(a)' ) '  R8PBU_PRINT prints an R8PBU matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call r8pbu_indicator ( n, mu, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU matrix:' )

  return
end
subroutine r8pbu_print_some_test ( )

!*****************************************************************************80
!
!! R8PBU_PRINT_SOME_TEST tests R8PBU_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 4

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8PBU_PRINT_SOME prints some of an R8PBU matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call r8pbu_indicator ( n, mu, a )

  call r8pbu_print_some ( n, mu, a, 4, 5, 8, 9, '  Row(4:8), Col(5:9):' )

  return
end
subroutine r8pbu_random_test ( )

!*****************************************************************************80
!
!! R8PBU_RANDOM_TEST tests R8PBU_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_RANDOM_TEST'
  write ( *, '(a)' ) '  R8PBU_RANDOM randomizes an R8PBU matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8pbu_random ( n, mu, seed, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU matrix:' )
 
  return
end
subroutine r8pbu_res_test ( )

!*****************************************************************************80
!
!! R8PBU_RES_TEST tests R8PBU_RES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) m
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n)

  m = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_RES_TEST'
  write ( *, '(a)' ) '  R8PBU_RES returns the residual b-A*x where A is'
  write ( *, '(a)' ) '  a positive definite symmetric band matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  call r8pbu_random ( n, mu, seed, a )
  call r8pbu_print ( n, mu, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Exact solution x:' )
!
!  Compute the right hand side.
!
  call r8pbu_mv ( m, n, mu, a, x, b )
  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Jostle the solution.
!
  call r8vec_uniform_01 ( n, seed, e )
  x2(1:n) = x(1:n) + 0.01D+00 * e(1:n)
  call r8vec_print ( n, x2, '  Approximate solution x2:' )
!
!  Compute the residual.
!
  call r8pbu_res ( m, n, mu, a, x2, b, r )
 
  call r8vec_print ( n, r, '  Residual r = b-A*x2:' )
 
  return
end
subroutine r8pbu_sl_test ( )

!*****************************************************************************80
!
!! R8PBU_SL_TEST tests R8PBU_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_SL_TEST'
  write ( *, '(a)' ) '  R8PBU_SL solves a linear system that was factored'
  write ( *, '(a)' ) '  by R8PBU_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  call r8pbu_random ( n, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  call r8pbu_mv ( n, n, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8pbu_fa ( n, mu, a, info )
!
!  Solve the linear system.
!
  call r8pbu_sl ( n, mu, a, b )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
 
  return
end
subroutine r8pbu_sor_test ( )

!*****************************************************************************80
!
!! R8PBU_SOR_TEST tests R8PBU_SOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itchk
  integer ( kind = 4 ) itknt
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) k
  real ( kind = 8 ) omega
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_SOR_TEST'
  write ( *, '(a)' ) '  R8PBU_SOR, SOR routine for iterative'
  write ( *, '(a)' ) '  solution of A*x=b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do k = 1, 3
 
    if ( k == 1 ) then
      omega = 0.25D+00
    else if ( k == 2 ) then
      omega = 0.75D+00
    else
      omega = 1.00D+00
    end if
!
!  Set matrix values.
!
    a(2,1:n) = 2.0D+00

    a(1,1) = 0.0D+00
    a(1,2:n) = -1.0D+00
!
!  Set the desired solution.
!
    do i = 1, n
      t = pi * real ( i - 1, kind = 8 ) / real ( n - 1, kind = 8 )
      x(i) = sin ( t )
    end do
!
!  Compute the right hand side.
!
    call r8pbu_mv ( n, n, mu, a, x, b ) 
!
!  Set the initial solution estimate.
!
    x(1:n) = 1.0D+00
 
    itchk = 1
    itmax = 8000
    eps = 0.0001D+00

    call r8pbu_sor ( n, mu, a, b, eps, itchk, itknt, itmax, omega, x )
!
!  Compute residual, A*x-b
!
    call r8pbu_mv ( n, n, mu, a, x, b2 )
 
    err = 0.0D+00
    do i = 1, n
      err = max ( err, abs ( b2(i) - b(i) ) )
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SOR iteration.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Relaxation factor OMEGA = ', omega
    write ( *, '(a,i8)' ) '  Iterations taken = ', itknt

    call r8vec_print_some ( n, x, 10, '  Solution:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Maximum error = ', err
 
  end do
 
  return
end
subroutine r8pbu_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8PBU_TO_R8GE_TEST tests R8PBU_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8PBU_TO_R8GE converts an R8PBU matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call r8pbu_indicator ( n, mu, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU matrix:' )

  call r8pbu_to_r8ge ( n, mu, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8pbu_zeros_test ( )

!*****************************************************************************80
!
!! R8PBU_ZEROS_TEST tests R8PBU_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_ZEROS_TEST'
  write ( *, '(a)' ) '  R8PBU_ZEROS sets up an R8PBU zero matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call r8pbu_zeros ( n, mu, a )

  call r8pbu_print ( n, mu, a, '  The R8PBU zero matrix:' )

  return
end

