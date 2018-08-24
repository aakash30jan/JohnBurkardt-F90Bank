program main

!*****************************************************************************80
!
!! MAIN is the main program for R83T_PRB.
!
!  Discussion:
!
!    R83T_PRB tests the R83T library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R83T library.'

  call r83t_cg_test ( )
  call r83t_dif2_test ( )
  call r83t_gs_sl_test ( )
  call r83t_indicator_test ( )
  call r83t_jac_sl_test ( )
  call r83t_mtv_test ( )
  call r83t_mv_test ( )
  call r83t_print_test ( )
  call r83t_print_some_test ( )
  call r83t_random_test ( )
  call r83t_res_test ( )
  call r83t_to_r8ge_test ( )
  call r83t_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r83t_cg_test ( )

!*****************************************************************************80
!
!! R83T_CG_TEST tests R83T_CG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(m,3)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_CG_TEST'
  write ( *, '(a)' ) '  R83T_CG solves an R83T linear system using'
  write ( *, '(a)' ) '  the conjugate gradient method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The right hand side B:' )

  x(1:n) = 0.0D+00
  call r83t_cg ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  return
end
subroutine r83t_dif2_test ( )

!*****************************************************************************80
!
!! R83T_DIF2_TEST tests R83T_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_DIF2_TEST'
  write ( *, '(a)' ) '  R83T_DIF2 sets up an R83T second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T second difference matrix:' )

  return
end
subroutine r83t_gs_sl_test ( )

!*****************************************************************************80
!
!! R83T_GS_SL_TEST tests R83T_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n,3)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: it_max = 25
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_GS_SL_TEST'
  write ( *, '(a)' ) '  R83T_GS_SL solves a linear system using'
  write ( *, '(a)' ) '  Gauss-Seidel iteration, with R83T matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83t_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83t_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83t_gs_sl ( n, a, b, x, it_max )

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83t_indicator_test ( )

!*****************************************************************************80
!
!! R83T_INDICATOR_TEST tests R83T_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83T_INDICATOR sets up an R83T indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T indicator matrix:' )

  return
end
subroutine r83t_jac_sl_test ( )

!*****************************************************************************80
!
!! R83T_JAC_SL_TEST tests R83T_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n,3)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: it_max = 25
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83T_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '  Jacobi iteration, with R83T matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83t_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83t_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83t_jac_sl ( n, a, b, x, it_max )

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83t_mtv_test ( )

!*****************************************************************************80
!
!! R83T_MTV_TEST tests R83T_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) a(m,3)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_MTV_TEST'
  write ( *, '(a)' ) '  R83T_MTV multiplies an R83T matrix transposed times a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r83t_mtv ( m, n, a, x, b )

  call r8vec_print ( n, b, '  The product b = A''*x:' )

  return
end
subroutine r83t_mv_test ( )

!*****************************************************************************80
!
!! R83T_MV_TEST tests R83T_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) a(m,3)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_MV_TEST'
  write ( *, '(a)' ) '  R83T_MV multiplies an R83T matrix times a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The product b = A*x:' )

  return
end
subroutine r83t_print_test ( )

!*****************************************************************************80
!
!! R83T_PRINT_TEST tests R83T_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRINT_TEST'
  write ( *, '(a)' ) '  R83T_PRINT prints an R83T matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix:' )

  return
end
subroutine r83t_print_some_test ( )

!*****************************************************************************80
!
!! R83T_PRINT_SOME_TEST tests R83T_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 9

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n

  n = 9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83T_PRINT_SOME prints some of an R83T matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print_some ( m, n, a, 3, 5, 6, 8, '  Rows 3:6, Cols 5:8:' )

  return
end
subroutine r83t_random_test ( )

!*****************************************************************************80
!
!! R83T_RANDOM_TEST tests R83T_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  n = 5
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_RANDOM_TEST'
  write ( *, '(a)' ) '  R83T_RANDOM sets up an R83T random matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_random ( m, n, seed, a )

  call r83t_print ( m, n, a, '  The R83T random matrix:' )

  return
end
subroutine r83t_res_test ( )

!*****************************************************************************80
!
!! R83T_RES_TEST tests R83T_RES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(m,3)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_RES_TEST'
  write ( *, '(a)' ) '  R83T_RES evaluates the residual given an approximate'
  write ( *, '(a)' ) '  solution of a linear system A*x=b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The right hand side B:' )

  x(1:n) = 0.0D+00
  call r83t_cg ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  call r83t_res ( m, n, a, x, b, r )

  call r8vec_print ( m, r, '  The residual b-A*x:' )

  return
end
subroutine r83t_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83T_TO_R8GE_TEST tests R83T_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a_r83t(m,3)
  real ( kind = 8 ) a_r8ge(m,n)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83T_TO_R8GE converts an R83T matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a_r83t )

  call r83t_print ( m, n, a_r83t, '  The R83T indicator matrix:' )

  call r83t_to_r8ge ( m, n, a_r83t, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE format matrix:' )

  return
end
subroutine r83t_zeros_test ( )

!*****************************************************************************80
!
!! R83T_ZEROS_TEST tests R83T_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ) a(m,3)
  integer ( kind = 4 ) n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_ZEROS_TEST'
  write ( *, '(a)' ) '  R83T_ZEROS sets up an R83T zero matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_zeros ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T zero matrix:' )

  return
end
