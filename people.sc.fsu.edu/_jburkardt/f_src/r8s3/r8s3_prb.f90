program main

!*****************************************************************************80
!
!! MAIN is the main program for R8S3_PRB.
!
!  Discussion:
!
!    R8S3_PRB tests the R8S3 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8S3 library.'

  call r8s3_diagonal_test ( )
  call r8s3_dif2_test ( )
  call r8s3_indicator_test ( )
  call r8s3_jac_sl_test ( )
  call r8s3_mtv_test ( )
  call r8s3_mv_test ( )
  call r8s3_print_test ( )
  call r8s3_print_some_test ( )
  call r8s3_random_test ( )
  call r8s3_write_test ( )
  call r8s3_read_test ( )
  call r8s3_res_test ( )
  call r8s3_to_r8ge_test ( )
  call r8s3_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8s3_diagonal_test ( )

!*****************************************************************************80
!
!! R8S3_DIAGONAL_TEST tests R8S3_DIAGONAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_num = 20

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    5, 6, 2, 2, 3, 4, 4, 5, 1, 6, &
    4, 6, 5, 1, 6, 3, 1, 2, 1, 3 /)
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 3, 4, 6, 5, 2, 6, 3, 1, 2, &
    4, 6, 5, 4, 4, 3, 6, 2, 3, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_DIAGONAL_TEST'
  write ( *, '(a)' ) '  R8S3_DIAGONAL rearranges an R8S3 matrix'
  write ( *, '(a)' ) '  so that the diagonal is listed first.'

  m = 6
  n = 6
  sym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_indicator ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before rearrangement:'
  write ( *, '(a)' ) '       K  ROW(K)  COL(K)      A(K)'
  write ( *, '(a)' ) ' '
  do k = 1, nz_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) k, row(k), col(k), a(k)
  end do

  call r8s3_diagonal ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After rearrangement:'
  write ( *, '(a)' ) '       K  ROW(K)  COL(K)      A(K)'
  write ( *, '(a)' ) ' '
  do k = 1, nz_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) k, row(k), col(k), a(k)
  end do

  return
end
subroutine r8s3_dif2_test ( )

!*****************************************************************************80
!
!! R8SP_DIF2_TEST tests R8S3_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10
  integer ( kind = 4 ), parameter :: nz_max = 3 * n_max - 2

  real ( kind = 8 ) a(nz_max)
  integer ( kind = 4 ) col(nz_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row(nz_max)
  integer ( kind = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_DIF2_TEST'
  write ( *, '(a)' ) '  R8S3_DIF2 sets an R8S3 matrix to the second difference.'

  m = 5
  n = 5
  nz_num = 3 * n - 2
  sym = 0

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix A:' )

  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix A:' )

  return
end
subroutine r8s3_indicator_test ( )

!*****************************************************************************80
!
!! R8S3_INDICATOR_TEST tests R8S3_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_num = 20

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    5, 6, 2, 2, 3, 4, 4, 5, 1, 6, &
    4, 6, 5, 1, 6, 3, 1, 2, 1, 3 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) sym
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 3, 4, 6, 5, 2, 6, 3, 1, 2, &
    4, 6, 5, 4, 4, 3, 6, 2, 3, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_INDICATOR'
  write ( *, '(a)' ) '  R8S3_INDICATOR sets up an R8S3 indicator matrix.'

  m = 6
  n = 6
  sym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_indicator ( m, n, nz_num, sym, row, col, a )

  call r8s3_print ( m, n, nz_num, sym, row, col, a, &
    '  The R8S3 indicator matrix:' )

  return
end
subroutine r8s3_jac_sl_test ( )

!*****************************************************************************80
!
!! R8S3_JAC_SL_TEST tests R8S3_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_max = 30

  real ( kind = 8 ) a(nz_max)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ) col(nz_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row(nz_max)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_JAC_SL_TEST'
  write ( *, '(a)' ) '  R8S3_JAC_SL uses Jacobi iteration to solve a linear system'
  write ( *, '(a)' ) '  with an R8S3 matrix.'

  m = 10
  n = 10
  sym = 0
  it_max = 25

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym
  write ( *, '(a,i8)' ) '  Iterations per call    = ', it_max
!
!  Set the matrix values.
!
  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solving A * x = b.'
  write ( *, '(a)' ) ' '
!
!  Set the desired solution.
!
  allocate ( x(1:n) )

  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )

  call r8s3_mv ( n, n, nz_num, sym, row, col, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r8s3_jac_sl ( n, nz_num, sym, row, col, a, b, x, it_max )

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  deallocate ( b )
  deallocate ( x )

  return
end
subroutine r8s3_mtv_test ( )

!*****************************************************************************80
!
!! R8S3_MTV_TEST tests R8S3_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: col(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable :: row(:)
  integer ( kind = 4 ) sym
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_MTV_TEST'
  write ( *, '(a)' ) '  R8S3_MTV computes b=A''*x, where A is an R8S3 matrix.'

  m = 5
  n = 4
  if ( m == n ) then
    nz_num = 3 * n - 2
  else
    nz_num = 3 * n - 1
  end if
  sym = 0

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  allocate ( x(1:m) )
  call r8vec_indicator1 ( m, x )
  call r8vec_print ( m, x, '  x:' )

  allocate ( b(1:n) )
  call r8s3_mtv ( m, n, nz_num, sym, row, col, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( row )
  deallocate ( x )
!
!  Try symmetric option.
!
  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  allocate ( x(1:m) )
  call r8vec_indicator1 ( m, x )
  call r8vec_print ( m, x, '  x:' )

  allocate ( b(1:n) )
  call r8s3_mtv ( m, n, nz_num, sym, row, col, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( row )
  deallocate ( x )

  return
end
subroutine r8s3_mv_test ( )

!*****************************************************************************80
!
!! R8S3_MV_TEST tests R8S3_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: col(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable :: row(:)
  integer ( kind = 4 ) sym
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_MV_TEST'
  write ( *, '(a)' ) '  R8S3_MV computes b=A*x, where A is an R8S3 matrix.'

  m = 5
  n = 4
  if ( m == n ) then
    nz_num = 3 * n - 2
  else
    nz_num = 3 * n - 1
  end if
  sym = 0

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  allocate ( b(1:m) )
  call r8s3_mv ( m, n, nz_num, sym, row, col, a, x, b )

  call r8vec_print ( m, b, '  b=A*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( row )
  deallocate ( x )
!
!  Try symmetric option.
!
  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  allocate ( b(1:m) )
  call r8s3_mv ( m, n, nz_num, sym, row, col, a, x, b )

  call r8vec_print ( m, b, '  b=A*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( row )
  deallocate ( x )

  return
end
subroutine r8s3_print_test ( )

!*****************************************************************************80
!
!! R8S3_PRINT_TEST tests R8S3_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_max = 20

  real ( kind = 8 ) a(nz_max)
  integer ( kind = 4 ) col(nz_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row(nz_max)
  integer ( kind = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_PRINT_TEST'
  write ( *, '(a)' ) '  R8S3_PRINT prints an R8S3 matrix.'

  m = 5
  n = 5
  nz_num = 3 * n - 2
  sym = 0

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix A:' )

  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix A:' )

  return
end
subroutine r8s3_print_some_test ( )

!*****************************************************************************80
!
!! R8S3_PRINT_SOME_TEST tests R8S3_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_max = 20

  real ( kind = 8 ) a(nz_max)
  integer ( kind = 4 ) col(nz_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row(nz_max)
  integer ( kind = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8S3_PRINT_SOME prints some of an R8S3 matrix.'

  m = 5
  n = 5
  nz_num = 3 * n - 2
  sym = 0

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_print_some ( m, n, nz_num, sym, row, col, a, 2, 3, 4, 5, '  Rows 2:4, Cols 3:5:' )

  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  call r8s3_print_some ( m, n, nz_num, sym, row, col, a, 2, 3, 4, 5, '  Rows 2:4, Cols 3:5:' )

  return
end
subroutine r8s3_random_test ( )

!*****************************************************************************80
!
!! R8S3_RANDOM_TEST tests R8S3_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_num = 20

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    5, 6, 2, 2, 3, 4, 4, 5, 1, 6, &
    4, 6, 5, 1, 6, 3, 1, 2, 1, 3 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sym
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 3, 4, 6, 5, 2, 6, 3, 1, 2, &
    4, 6, 5, 4, 4, 3, 6, 2, 3, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_RANDOM'
  write ( *, '(a)' ) '  R8S3_RANDOM randomizes an R8S3 matrix.'

  m = 6
  n = 6
  sym = 0
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_random ( m, n, nz_num, sym, row, col, seed, a )

  call r8s3_print ( m, n, nz_num, sym, row, col, a, &
    '  The R8S3 indicator matrix:' )

  return
end
subroutine r8s3_read_test ( )

!*****************************************************************************80
!
!! R8S3_READ_TEST tests R8S3_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  character ( len = 80 ) input_file
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row

  input_file = 'r8s3_matrix.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_READ_TEST'
  write ( *, '(a)' ) '  R8S3_READ_SIZE reads the size of an R8S3 matrix.'
  write ( *, '(a)' ) '  R8S3_READ reads an R8S3 matrix from a file.'

  call r8s3_read_size ( input_file, m, n, nz_num )

  sym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_read ( input_file, m, n, nz_num, row, col, a )

  sym = 0

  call r8s3_print_some ( m, n, nz_num, sym, row, col, a, 1, 1, &
    10, 10, '  Initial 10x10 block of recovered R8S3 matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Deleting the matrix data file "' &
    // trim ( input_file ) // '".'

  call file_delete ( input_file )

  deallocate ( row )
  deallocate ( col )
  deallocate ( a )

  return
end
subroutine r8s3_res_test ( )

!*****************************************************************************80
!
!! R8S3_RES_TEST tests R8S3_RES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: col(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  real ( kind = 8 ), allocatable :: r(:)
  integer ( kind = 4 ), allocatable :: row(:)
  integer ( kind = 4 ) sym
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_RES_TEST'
  write ( *, '(a)' ) '  R8S3_RES computes r=b-A*x, where A is an R8S3 matrix.'

  m = 5
  n = 4
  if ( m == n ) then
    nz_num = 3 * n - 2
  else
    nz_num = 3 * n - 1
  end if
  sym = 0

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  allocate ( b(1:m) )
  call r8s3_mv ( m, n, nz_num, sym, row, col, a, x, b )

  allocate ( r(1:m) )
  call r8s3_res ( m, n, nz_num, sym, row, col, a, x, b, r )

  call r8vec_print ( m, r, '  r=b-A*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( r )
  deallocate ( row )
  deallocate ( x )
!
!  Try symmetric option.
!
  m = 5
  n = 5
  nz_num = 2 * n - 1
  sym = 1

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  allocate ( b(1:m) )
  call r8s3_mv ( m, n, nz_num, sym, row, col, a, x, b )

  allocate ( r(1:m) )
  call r8s3_res ( m, n, nz_num, sym, row, col, a, x, b, r )

  call r8vec_print ( m, r, '  r=b-A*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( r )
  deallocate ( row )
  deallocate ( x )

  return
end
subroutine r8s3_write_test ( )

!*****************************************************************************80
!
!! R8S3_WRITE_TEST tests R8S3_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 100
  integer ( kind = 4 ), parameter :: nz_num = 3 * n_max - 2

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 80 ) output_file
  integer ( kind = 4 ) row(nz_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_WRITE_TEST'
  write ( *, '(a)' ) '  R8S3_WRITE writes an R8S3 matrix to a file.'

  m = 100
  n = 100
  sym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym
!
!  Set the matrix values.
!
  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a )
!
!  Print some of the matrix.
!
  call r8s3_print_some ( m, n, nz_num, sym, row, col, a, 1, 1, &
    10, 10, '  Initial 10x10 block of R8S3 matrix:' )
!
!  Write the matrix to a file.
!
  output_file = 'r8s3_matrix.txt'

  call r8s3_write ( m, n, nz_num, sym, row, col, a, output_file )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8S3_WRITE wrote the matrix data to "' &
    // trim ( output_file ) // '".'

  return
end
subroutine r8s3_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8S3_TO_R8GE_TEST tests R8S3_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a_r8s3(:)
  real ( kind = 8 ), allocatable :: a_r8ge(:,:)
  integer ( kind = 4 ), allocatable ::  col(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable ::  row(:)
  integer ( kind = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8S3_TO_R8GE converts an R8S3 matrix to R8GE format.'

  m = 5
  n = 5
  nz_num = 3 * n - 2
  sym = 0
  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a_r8s3(1:nz_num) )

  call r8s3_dif2 ( m, n, nz_num, sym, row, col, a_r8s3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a_r8s3, '  R8S3 matrix A:' )

  allocate ( a_r8ge(1:m,1:n) )

  call r8s3_to_r8ge ( m, n, nz_num, sym, row, col, a_r8s3, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  R8GE matrix A:' )

  deallocate ( a_r8ge )
  deallocate ( a_r8s3 )
  deallocate ( col )
  deallocate ( row )

  return
end
subroutine r8s3_zeros_test ( )

!*****************************************************************************80
!
!! R8S3_ZEROS_TEST tests R8S3_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nz_max = 30

  real ( kind = 8 ) a(nz_max)
  integer ( kind = 4 ) col(nz_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row(nz_max)
  integer ( kind = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8S3_ZEROS_TEST'
  write ( *, '(a)' ) '  R8S3_ZEROS sets an R8S3 matrix to zeros.'

  m = 5
  n = 5
  nz_num = 13
  sym = 0
  row(1:nz_num) = (/ 1, 2, 3, 4, 5, 1, 2, 3, 3, 4, 4, 5, 5 /)
  col(1:nz_num) = (/ 1, 2, 3, 4, 5, 2, 4, 2, 5, 1, 5, 1, 3 /)

  call r8s3_zeros ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM =    ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix:' )

  m = 5
  n = 5
  nz_num = 10
  sym = 1
  row(1:nz_num) = (/ 1, 2, 3, 4, 5, 3, 4, 5, 5 /)
  col(1:nz_num) = (/ 1, 2, 3, 4, 5, 2, 1, 1, 3 /)

  call r8s3_zeros ( m, n, nz_num, sym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M         = ', m
  write ( *, '(a,i8)' ) '  Matrix order N         = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a,i8)' ) '  Symmetry option SYM    = ', sym

  call r8s3_print ( m, n, nz_num, sym, row, col, a, '  R8S3 matrix:' )

  return
end
