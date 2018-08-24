program main

!*****************************************************************************80
!
!! R8LTT_TEST tests the R8LTT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the R8LTT library.'

  call r8ltt_det_test ( )
  call r8ltt_indicator_test ( )
  call r8ltt_inverse_test ( )
  call r8ltt_mm_test ( )
  call r8ltt_mtm_test ( )
  call r8ltt_mtv_test ( )
  call r8ltt_mv_test ( )
  call r8ltt_print_test ( )
  call r8ltt_print_some_test ( )
  call r8ltt_random_test ( )
  call r8ltt_sl_test ( )
  call r8ltt_slt_test ( )
  call r8ltt_to_r8ge_test ( )
  call r8ltt_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine r8ltt_det_test ( )

!*****************************************************************************80
!
!! R8LTT_DET_TEST tests R8LTT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_DET_TEST'
  write ( *, '(a)' ) '  R8LTT_DET computes the determinant of an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call  r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  The matrix:' )
!
!  Compute the determinant.
!
  call r8ltt_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  The determinant = ', det

  return
end
subroutine r8ltt_indicator_test ( )

!*****************************************************************************80
!
!! R8LTT_INDICATOR_TEST tests R8LTT_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8LTT_INDICATOR sets up an indicator matrix in R8LTT format'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_indicator ( n, a )

  call r8ltt_print ( n, a, '  The indicator matrix:' )

  return
end
subroutine r8ltt_inverse_test ( )

!*****************************************************************************80
!
!! R8LTT_INVERSE_TEST tests R8LTT_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_INVERSE_TEST'
  write ( *, '(a)' ) '  R8LTT_INVERSE computes the inverse of an R8LTT matrix.'

  call r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  The matrix A:' )
!
!  Compute the inverse matrix B.
!
  call r8ltt_inverse ( n, a, b )

  call r8ltt_print ( n, b, '  The inverse matrix B:' )
!
!  Check
!
  call r8ltt_mm ( n, a, b, c )

  call r8ltt_print ( n, c, '  The product A * B:' )

  return
end
subroutine r8ltt_mm_test ( )

!*****************************************************************************80
!
!! R8LTT_MM_TEST tests R8LTT_MM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_ge(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b_ge(n,n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) c_ge(n,n)
  integer ( kind = 4 ) seed
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_MM_TEST'
  write ( *, '(a)' ) '  R8LTT_MM computes C = A * B for R8LTT matrices.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ltt_random ( n, seed, a )
  call r8ltt_print ( n, a, '  Factor A:' )
  call r8ltt_random ( n, seed, b )
  call r8ltt_print ( n, b, '  Factor B:' )
  call r8ltt_mm ( n, a, b, c )
  call r8ltt_print ( n, c, '  The product C = A * B' )

  call r8ltt_to_r8ge ( n, a, a_ge )
  call r8ltt_to_r8ge ( n, b, b_ge )
  c_ge = matmul ( a_ge(1:n,1:n), b_ge(1:n,1:n) )
  call r8ge_print ( n, n, c_ge, '  The R8GE product C:' )

  return
end
subroutine r8ltt_mtm_test ( )

!*****************************************************************************80
!
!! R8LTT_MTM_TEST tests R8LTT_MTM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_ge(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b_ge(n,n)
  real ( kind = 8 ) c_ge(n,n)
  integer ( kind = 4 ) seed

  seed = 123456789
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_MTM_TEST'
  write ( *, '(a)' ) '  R8LTT_MTM computes C = A'' * B for R8LTT matrices.'
 
  call r8ltt_random ( n, seed, a )
  call r8ltt_print ( n, a, '  The matrix A:' )
  call r8ltt_random ( n, seed, b )
  call r8ltt_print ( n, b, '  The matrix B:' )

  call r8ltt_mtm ( n, a, b, c_ge )
  call r8ge_print ( n, n, c_ge, '  The product C = A'' * B:' )
!
!  Compare with an R8GE calculation.
!
  call r8ltt_to_r8ge ( n, a, a_ge )
  call r8ltt_to_r8ge ( n, b, b_ge )
  c_ge = matmul ( transpose ( a_ge ), b_ge )
  call r8ge_print ( n, n, c_ge, '  The R8GE product C = A'' * B:' )

  return
end
subroutine r8ltt_mtv_test ( )

!*****************************************************************************80
!
!! R8LTT_MTV_TEST tests R8LTT_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_MTV_TEST'
  write ( *, '(a)' ) '  R8LTT_MTV computes a matrix product b=A''*x for an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_indicator ( n, a )
  call r8ltt_print ( n, a, '  The matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector X:' )

  call r8ltt_mtv ( n, a, x, b )
  call r8vec_print ( n, b, '  The vector b=A''*x:' )

  return
end
subroutine r8ltt_mv_test ( )

!*****************************************************************************80
!
!! R8LTT_MV_TEST tests R8LTT_MV
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_MV_TEST'
  write ( *, '(a)' ) '  R8LTT_MV computes a product b=A*x for an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_indicator ( n, a )
  call r8ltt_print ( n, a, '  The R8LTT matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )

  call r8ltt_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Vector b = A*x:' )

  return
end
subroutine r8ltt_print_some_test ( )

!*****************************************************************************80
!
!! R8LTT_PRINT_SOME_TEST tests R8LTT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8LTT_PRINT_SOME prints some of an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_indicator ( n, a )

  call r8ltt_print_some ( n, a, 2, 1, 5, 4, '  Some of the matrix:' )

  return
end
subroutine r8ltt_print_test ( )

!*****************************************************************************80
!
!! R8LTT_PRINT_TEST tests R8LTT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_PRINT_TEST'
  write ( *, '(a)' ) '  R8LTT_PRINT prints an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_indicator ( n, a )

  call r8ltt_print ( n, a, '  The matrix:' )

  return
end
subroutine r8ltt_random_test ( )

!*****************************************************************************80
!
!! R8LTT_RANDOM_TEST tests R8LTT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_RANDOM_TEST'
  write ( *, '(a)' ) '  R8LTT_RANDOM randomizes an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  Matrix A:' )

  return
end
subroutine r8ltt_sl_test ( )

!*****************************************************************************80
!
!! R8LTT_SL_TEST tests R8LTT_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_SL_TEST'
  write ( *, '(a)' ) '  R8LTT_SL solves a linear system A*x=b with R8LTT matrix'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8ltt_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Solve the linear system.
!
  call r8ltt_sl ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution x:' )

  return
end
subroutine r8ltt_slt_test ( )

!*****************************************************************************80
!
!! R8LTT_SLT_TEST tests R8LTT_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_SLT_TEST'
  write ( *, '(a)' ) '  R8LTT_SLT solves a linear system A''x=b with R8LTT matrix'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8ltt_mtv ( n, a, x, b )

  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Solve the linear system.
!
  call r8ltt_slt ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution x:' )

  return
end
subroutine r8ltt_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8LTT_TO_R8GE_TEST tests R8LTT_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_ge(n,n)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8LTT_TO_R8GE converts an R8LTT matrix to R8GE format.'

  call r8ltt_random ( n, seed, a )

  call r8ltt_print ( n, a, '  The random R8LTT matrix:' )
 
  call r8ltt_to_r8ge ( n, a, a_ge )

  call r8ge_print ( n, n, a_ge, '  The R8GE matrix:' )

  return
end
subroutine r8ltt_zeros_test ( )

!*****************************************************************************80
!
!! R8LTT_ZEROS_TEST tests R8LTT_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8LTT_ZEROS_TEST'
  write ( *, '(a)' ) '  R8LTT_ZEROS zeros out space for an R8LTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8ltt_zeros ( n, a )

  call r8ltt_print ( n, a, '  The matrix:' )

  return
end
