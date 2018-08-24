program main

!*****************************************************************************80
!
!! MAIN is the main program for R8TO_PRB.
!
!  Discussion:
!
!    R8TO_PRB tests the R8TO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8TO library.'

  call r8to_dif2_test ( )
  call r8to_indicator_test ( )
  call r8to_mtv_test ( )
  call r8to_mv_test ( )
  call r8to_print_test ( )
  call r8to_print_some_test ( )
  call r8to_random_test ( )
  call r8to_sl_test ( )
  call r8to_slt_test ( )
  call r8to_to_r8ge_test ( )
  call r8to_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8to_dif2_test ( )

!*****************************************************************************80
!
!! R8TO_DIF2_TEST tests R8TO_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_DIF2_TEST'
  write ( *, '(a)' ) '  R8TO_DIF2 sets the second difference as an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8to_dif2 ( n, a )

  call r8to_print ( n, a, '  The R8TO matrix:' )

  return
end
subroutine r8to_indicator_test ( )

!*****************************************************************************80
!
!! R8TO_INDICATOR_TEST tests R8TO_INDICATOR.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8TO_INDICATOR sets up an R8TO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8to_indicator ( n, a )

  call r8to_print ( n, a, '  The R8TO indicator matrix:' )

  return
end
subroutine r8to_mtv_test ( )

!*****************************************************************************80
!
!! R8TO_MTV_TEST tests R8TO_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_MTV_TEST'
  write ( *, '(a)' ) '  R8TO_MTV computes b=A''*x=b, where A is an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8to_indicator ( n, a )
  call r8to_print ( n, a, '  The Toeplitz matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )
!
!  Compute b=A'*x.
!
  call r8to_mtv ( n, a, x, b )
  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8to_mv_test ( )

!*****************************************************************************80
!
!! R8TO_MV_TEST tests R8TO_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_MV_TEST'
  write ( *, '(a)' ) '  R8TO_MV computes b=A*x=b, where A is an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8to_indicator ( n, a )
  call r8to_print ( n, a, '  The Toeplitz matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )
!
!  Compute b=A*x.
!
  call r8to_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8to_print_test ( )

!*****************************************************************************80
!
!! R8TO_PRINT_TEST tests R8TO_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_PRINT_TEST'
  write ( *, '(a)' ) '  R8TO_PRINT prints an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8to_indicator ( n, a )

  call r8to_print ( n, a, '  The R8TO indicator matrix:' )

  return
end
subroutine r8to_print_some_test ( )

!*****************************************************************************80
!
!! R8TO_PRINT_SOME_TEST tests R8TO_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8TO_PRINT_SOME prints some of an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8to_indicator ( n, a )

  call r8to_print_some ( n, a, 2, 1, 5, 3, '  Rows2:5, Cols 1:3' )

  return
end
subroutine r8to_random_test ( )

!*****************************************************************************80
!
!! R8TO_RANDOM_TEST tests R8TO_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_RANDOM_TEST'
  write ( *, '(a)' ) '  R8TO_RANDOM randomizes an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call r8to_random ( n, seed, a )

  call r8to_print ( n, a, '  The R8TO matrix:' )

  return
end
subroutine r8to_sl_test ( )

!*****************************************************************************80
!
!! R8TO_SL_TEST tests R8TO_SL.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_SL_TEST'
  write ( *, '(a)' ) '  R8TO_SL solves A*x=b, where A is an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8to_random ( n, seed, a )

  call r8to_print ( n, a, '  The Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8to_mv ( n, a, x, b )
!
!  Solve the linear system.
!
  call r8to_sl ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )
 
  return
end
subroutine r8to_slt_test ( )

!*****************************************************************************80
!
!! R8TO_SLT_TEST tests R8TO_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_SLT_TEST'
  write ( *, '(a)' ) '  R8TO_SLT solves A''*x=b, where A is an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8to_random ( n, seed, a )

  call r8to_print ( n, a, '  The Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8to_mtv ( n, a, x, b )
!
!  Solve the linear system.
!
  call r8to_slt ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution to transposed system:' )
 
  return
end
subroutine r8to_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8TO_TO_TO_R8GE_TEST tests R8TO_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a_r8ge(n,n)
  real ( kind = 8 ) a_r8to(2*n-1)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8TO_TO_R8GE converts a matrix from R8TO to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8to_random ( n, seed, a_r8to )

  call r8to_print ( n, a_r8to, '  The matrix in R8TO format:' )

  call r8to_to_r8ge ( n, a_r8to, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The matrix in R8GE format:' )
 
  return
end
subroutine r8to_zeros_test ( )

!*****************************************************************************80
!
!! R8TO_ZEROS_TEST tests R8TO_ZEROS.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8TO_ZEROS_TEST'
  write ( *, '(a)' ) '  R8TO_ZEROS zeros an R8TO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8to_zeros ( n, a )

  call r8to_print ( n, a, '  The R8TO matrix:' )

  return
end

