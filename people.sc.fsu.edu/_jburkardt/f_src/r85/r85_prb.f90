program main

!*****************************************************************************80
!
!! MAIN is the main program for R85_PRB.
!
!  Discussion:
!
!    R85_PRB tests the R85 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R85 library.'

  call r85_dif2_test ( )
  call r85_indicator_test ( )
  call r85_mtv_test ( )
  call r85_mv_test ( )
  call r85_np_fs_test ( )
  call r85_print_test ( )
  call r85_print_some_test ( )
  call r85_random_test ( )
  call r85_to_r8ge_test ( )
  call r85_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R85_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine r85_dif2_test ( )

!*****************************************************************************80
!
!! R85_DIF2_TEST tests R85_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_DIF2_TEST'
  write ( *, '(a)' ) '  R85_DIF2 sets up an R85 second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_dif2 ( n, a )

  call r85_print ( n, a, '  The R85 second difference matrix:' )

  return
end
subroutine r85_indicator_test ( )

!*****************************************************************************80
!
!! R85_INDICATOR_TEST tests R85_INDICATOR.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_INDICATOR_TEST'
  write ( *, '(a)' ) '  R85_INDICATOR sets up an R85 indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 indicator matrix:' )

  return
end
subroutine r85_mtv_test ( )

!*****************************************************************************80
!
!! R85_MTV_TEST tests R85_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_MTV_TEST'
  write ( *, '(a)' ) '  R85_MTV computes b=A''*x, where A is an R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 indicator matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r85_mtv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x' )

  return
end
subroutine r85_mv_test ( )

!*****************************************************************************80
!
!! R85_MV_TEST tests R85_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_MV_TEST'
  write ( *, '(a)' ) '  R85_MV computes b=A*x, where A is an R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 indicator matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r85_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x' )

  return
end
subroutine r85_np_fs_test ( )

!*****************************************************************************80
!
!! R85_NP_FS_TEST tests R85_NP_FS.
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

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_NP_FS_TEST'
  write ( *, '(a)' ) '  R85_NP_FS factors and solves an R85'
  write ( *, '(a)' ) '  linear system, with no pivoting.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix to a random value.
!
  call r85_random ( n, seed, a )

  call r85_print ( n, a, '  The pentadiagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute b = A * x.
!
  call r85_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the system.
!
  call r85_np_fs ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine r85_print_test ( )

!*****************************************************************************80
!
!! R85_PRINT_TEST tests R85_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_PRINT_TEST'
  write ( *, '(a)' ) '  R85_PRINT prints an R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 matrix:' )

  return
end
subroutine r85_print_some_test ( )

!*****************************************************************************80
!
!! R85_PRINT_SOME_TEST tests R85_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R85_PRINT_SOME prints some of an R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print_some ( n, a, 2, 3, 6, 5, '  Rows 2-6, Cols 3-5:' )

  return
end
subroutine r85_random_test ( )

!*****************************************************************************80
!
!! R85_RANDOM_TEST tests R85_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_RANDOM_TEST'
  write ( *, '(a)' ) '  R85_RANDOM sets up a random R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call r85_random ( n, seed, a )

  call r85_print ( n, a, '  The random R85 matrix:' )

  return
end
subroutine r85_to_r8ge_test ( )

!*****************************************************************************80
!
!! R85_TO_R8GE_TEST tests R85_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R85_TO_R8GE converts an R85 matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 matrix:' )

  call r85_to_r8ge ( n, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, "  The R8GE format" )

  return
end
subroutine r85_zeros_test ( )

!*****************************************************************************80
!
!! R85_ZEROS_TEST tests R85_ZEROS.
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

  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R85_ZEROS_TEST'
  write ( *, '(a)' ) '  R85_ZEROS zeros an R85 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_zeros ( n, a )

  call r85_print ( n, a, '  The zero R85 matrix:' )

  return
end

