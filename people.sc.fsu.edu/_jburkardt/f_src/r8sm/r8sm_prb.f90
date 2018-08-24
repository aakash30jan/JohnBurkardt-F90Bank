program main

!*****************************************************************************80
!
!! MAIN is the main program for R8SM_PRB.
!
!  Discussion:
!
!    R8SM_PRB tests the R8SM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8SM library.'

  call r8sm_indicator_test ( )
  call r8sm_ml_test ( )
  call r8sm_mtv_test ( )
  call r8sm_mv_test ( )
  call r8sm_print_test ( )
  call r8sm_print_some_test ( )
  call r8sm_random_test ( )
  call r8sm_sl_test ( )
  call r8sm_slt_test ( )
  call r8sm_to_r8ge_test ( )
  call r8sm_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8sm_indicator_test ( )

!*****************************************************************************80
!
!! R8SM_INDICATOR_TEST tests R8SM_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8SM_INDICATOR sets up an R8SM indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM indicator matrix:' )

  return
end
subroutine r8sm_ml_test ( )

!*****************************************************************************80
!
!! R8SM_ML_TEST tests R8SM_ML.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = m

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_ML_TEST'
  write ( *, '(a)' ) '  R8SM_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  where A is a Sherman Morrison matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r8sm_random ( m, n, seed, a, u, v )

    call r8sm_print ( m, n, a, u, v, '  The Sherman Morrison matrix:' )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8sm_mv ( m, n, a, u, v, x, b )
    else
      call r8sm_mtv ( m, n, a, u, v, x, b )
    end if
!
!  Factor the matrix.
!
    call r8ge_fa ( n, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8SM_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8sm_ml ( n, a, u, v, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8sm_mtv_test ( )

!*****************************************************************************80
!
!! R8SM_MTV_TEST tests R8SM_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_MTV_TEST'
  write ( *, '(a)' ) '  R8SM_MTV computes A''*x=b, where A is an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r8sm_mtv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8sm_mv_test ( )

!*****************************************************************************80
!
!! R8SM_MV_TEST tests R8SM_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_MV_TEST'
  write ( *, '(a)' ) '  R8SM_MV computes A*x=b, where A is an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sm_mv ( m, n, a, u, v, x, b )

  call r8vec_print ( m, b, '  The product b=A*x:' )

  return
end
subroutine r8sm_print_test ( )

!*****************************************************************************80
!
!! R8SM_PRINT_TEST tests R8SM_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRINT_TEST'
  write ( *, '(a)' ) '  R8SM_PRINT prints an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  return
end
subroutine r8sm_print_some_test ( )

!*****************************************************************************80
!
!! R8SM_PRINT_SOME_TEST tests R8SM_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 9
  integer ( kind = 4 ), parameter :: n = 9

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8SM_PRINT_SOME prints some of an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print_some ( m, n, a, u, v, 2, 3, 5, 7, '  Rows 2-5, Cols 3-7:' )

  return
end
subroutine r8sm_random_test ( )

!*****************************************************************************80
!
!! R8SM_RANDOM_TEST tests R8SM_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_RANDOM_TEST'
  write ( *, '(a)' ) '  R8SM_RANDOM sets up a random R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  seed = 123456789
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  return
end
subroutine r8sm_sl_test ( )

!*****************************************************************************80
!
!! R8SM_SL_TEST tests R8SM_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = m

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_SL_TEST'
  write ( *, '(a)' ) '  R8SM_SL solves B*x=b, where B is an R8SM matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8sm_mv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
  call r8ge_fa ( n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8sm_sl ( n, a, u, v, b, ierror, pivot )
 
  call r8vec_print ( n, b, '  Solution to A * X = B:' )

  return
end
subroutine r8sm_slt_test ( )

!*****************************************************************************80
!
!! R8SM_SLT_TEST tests R8SM_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = m

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_SLT_TEST'
  write ( *, '(a)' ) '  R8SM_SLT solves B''*x=b, where B is an R8SM matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8sm_mtv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
  call r8ge_fa ( n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SLT_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8sm_slt ( n, a, u, v, b, ierror, pivot )

  call r8vec_print ( n, b, '  Solution to A'' * X = B:' )

  return
end
subroutine r8sm_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8SM_TO_R8GE_TEST tests R8SM_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_r8ge(m,n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8SM_TO_R8GE converts an R8SM matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8sm_to_r8ge ( m, n, a, u, v, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8sm_zeros_test ( )

!*****************************************************************************80
!
!! R8SM_ZEROS_TEST tests R8SM_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_ZEROS_TEST'
  write ( *, '(a)' ) '  R8SM_ZEROS sets up a zero R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_zeros ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM zero matrix:' )

  return
end
