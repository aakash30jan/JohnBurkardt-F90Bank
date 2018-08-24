program main

!*****************************************************************************80
!
!! MAIN is the main program for R83S_PRB.
!
!  Discussion:
!
!    R83S_PRB tests the R83S library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R83S library.'

  call r83s_cg_test ( )
  call r83s_dif2_test ( )
  call r83s_gs_sl_test ( )
  call r83s_indicator_test ( )
  call r83s_jac_sl_test ( )
  call r83s_mtv_test ( )
  call r83s_mv_test ( )
  call r83s_print_test ( )
  call r83s_print_some_test ( )
  call r83s_random_test ( )
  call r83s_res_test ( )
  call r83s_to_r8ge_test ( )
  call r83s_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r83s_cg_test ( )

!*****************************************************************************80
!
!! R83S_CG_TEST tests R83S_CG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) e_norm
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ) r_norm
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) r8vec_norm
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_CG_TEST'
  write ( *, '(a)' ) '  R83S_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '  to solve a linear system with an R83S matrix.'

  seed = 123456789
  n = 10
!
!  Let A be the -1 2 -1 matrix.
!
  call r83s_dif2 ( n, n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83s_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r83s_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83s_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83s_dif2_test ( )

!*****************************************************************************80
!
!! R83S_DIF2_TEST tests R83S_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_DIF2_TEST'
  write ( *, '(a)' ) '  R83S_DIF2 sets an R83S matrix to the second difference.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    call r83s_dif2 ( m, n, a )

    call r83s_print ( m, n, a, '  Second difference in R83S format:' )

  end do

  return
end
subroutine r83s_gs_sl_test ( )

!*****************************************************************************80
!
!! R83S_GS_SL_TEST tests R83S_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: it_max = 25
  real ( kind = 8 ) :: tol = 0.000001D+00
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_GS_SL_TEST'
  write ( *, '(a)' ) '  R83S_GS_SL solves a linear system using'
  write ( *, '(a)' ) '  Gauss-Seidel iteration, with R83S matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83s_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83s_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83s_gs_sl ( n, a, b, x, tol, it_max, it, diff )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
    write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83s_indicator_test ( )

!*****************************************************************************80
!
!! R83S_INDICATOR_TEST tests R83S_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83S_INDICATOR sets an R83S matrix to an indicator matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    call r83s_indicator ( m, n, a )

    call r83s_print ( m, n, a, '  R83S indicator matrix:' )

  end do

  return
end
subroutine r83s_jac_sl_test ( )

!*****************************************************************************80
!
!! R83S_JAC_SL_TEST tests R83S_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: it_max = 25
  real ( kind = 8 ) :: tol = 0.000001D+00
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83S_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '  Jacobi iteration, with R83S matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83s_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83s_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83s_jac_sl ( n, a, b, x, tol, it_max, it, diff )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
    write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83s_mtv_test ( )

!*****************************************************************************80
!
!! R83S_MTV_TEST tests R83S_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a_83s(3)
  real ( kind = 8 ), allocatable :: a_ge(:,:)
  real ( kind = 8 ), allocatable :: ax_83s(:)
  real ( kind = 8 ), allocatable :: ax_ge(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_MTV_TEST'
  write ( *, '(a)' ) '  R83S_MTV computes b=A''*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_ge(m,n) )
    allocate ( ax_83s(n) )
    allocate ( ax_ge(n) )
    allocate ( x(m) )

    seed = 123456789
    call r83s_random ( m, n, seed, a_83s )
    call r8vec_indicator1 ( m, x )
    call r83s_mtv ( m, n, a_83s, x, ax_83s )
    call r83s_to_r8ge ( m, n, a_83s, a_ge )
    call r8ge_mtv ( m, n, a_ge, x, ax_ge )
    call r8vec2_print ( n, ax_83s, ax_ge, '  Product comparison:' )

    deallocate ( a_ge )
    deallocate ( ax_83s )
    deallocate ( ax_ge )
    deallocate ( x )

  end do

  return
end
subroutine r83s_mv_test ( )

!*****************************************************************************80
!
!! R83S_MV_TEST tests R83S_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a_83s(3)
  real ( kind = 8 ), allocatable :: a_ge(:,:)
  real ( kind = 8 ), allocatable :: ax_83s(:)
  real ( kind = 8 ), allocatable :: ax_ge(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_MV_TEST'
  write ( *, '(a)' ) '  R83S_MV computes b=A*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_ge(m,n) )
    allocate ( ax_83s(m) )
    allocate ( ax_ge(m) )
    allocate ( x(n) )

    seed = 123456789
    call r83s_random ( m, n, seed, a_83s )
    call r8vec_indicator1 ( n, x )
    call r83s_mv ( m, n, a_83s, x, ax_83s )
    call r83s_to_r8ge ( m, n, a_83s, a_ge )
    call r8ge_mv ( m, n, a_ge, x, ax_ge )
    call r8vec2_print ( m, ax_83s, ax_ge, '  Product comparison:' )

    deallocate ( a_ge )
    deallocate ( ax_83s )
    deallocate ( ax_ge )
    deallocate ( x )

  end do

  return
end
subroutine r83s_print_test ( )

!*****************************************************************************80
!
!! R83S_PRINT_TEST tests R83S_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_PRINT_TEST'
  write ( *, '(a)' ) '  R83S_PRINT prints an R83S matrix.'

  m = 5
  n = 4

  call r83s_indicator ( m, n, a )

  call r83s_print ( m, n, a, '  R83S indicator matrix:' )

  return
end
subroutine r83s_print_some_test ( )

!*****************************************************************************80
!
!! R83S_PRINT_SOME_TEST tests R83S_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83S_PRINT_SOME prints some of an R83S matrix.'
!
!  Set the matrix.
!
  call r83s_indicator ( m, n, a )

  call r83s_print_some ( m, n, a, 2, 2, 5, 4, '  Rows 2-5, Cols 2-4:' )

  return
end
subroutine r83s_random_test ( )

!*****************************************************************************80
!
!! R83S_RANDOM_TEST tests R83S_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_RANDOM_TEST'
  write ( *, '(a)' ) '  R83S_RANDOM randomizes an R83S matrix.'

  m = 5
  n = 4
  seed = 123456789

  call r83s_random ( m, n, seed, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

  return
end
subroutine r83s_res_test ( )

!*****************************************************************************80
!
!! R83S_RES_TEST tests R83S_RES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_RES_TEST'
  write ( *, '(a)' ) '  R83S_RES computes b-A*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( b(1:m) )
    allocate ( r(1:m) )
    allocate ( x(1:n) )

    seed = 123456789
    call r83s_random ( m, n, seed, a )
    call r8vec_indicator1 ( n, x )
    call r83s_mv ( m, n, a, x, b )
    call r83s_res ( m, n, a, x, b, r )
    call r8vec_print ( m, r, '  Residual A*x-b:' )

    deallocate ( b )
    deallocate ( r )
    deallocate ( x )

  end do

  return
end
subroutine r83s_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83S_TO_R8GE_TEST tests R83S_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) a_ge(m,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83S_TO_R8GE converts an R83S matrix to R8GE format.'

  seed = 123456789

  call r83s_random ( m, n, seed, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

  call r83s_to_r8ge ( m, n, a, a_ge )

  call r8ge_print ( m, n, a_ge, '  R8GE matrix:' )

  return
end
subroutine r83s_zeros_test ( )

!*****************************************************************************80
!
!! R83S_ZEROS_TEST tests R83S_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_ZEROS_TEST'
  write ( *, '(a)' ) '  R83S_ZEROS zeros an R83S matrix.'

  m = 5
  n = 4

  call r83s_zeros ( m, n, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

  return
end

