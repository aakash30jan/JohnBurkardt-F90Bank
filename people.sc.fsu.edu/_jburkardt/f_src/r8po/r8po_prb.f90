program main

!*****************************************************************************80
!
!! MAIN is the main program for R8PO_PRB.
!
!  Discussion:
!
!    R8PO_PRB tests the R8PO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8PO library.'

  call r8ge_to_r8po_test ( )

  call r8po_det_test ( )
  call r8po_dif2_test ( )
  call r8po_fa_test ( )
  call r8po_indicator_test ( )
  call r8po_inverse_test ( )
  call r8po_ml_test ( )
  call r8po_mm_test ( )
  call r8po_mv_test ( )
  call r8po_print_test ( )
  call r8po_print_some_test ( )
  call r8po_random_test ( )
  call r8po_sl_test ( )
  call r8po_to_r8ge_test ( )
  call r8po_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8ge_to_r8po_test ( )

!*****************************************************************************80
!
!! R8GE_TO_R8PO_TEST tests R8GE_TO_R8PO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_TO_R8PO_TEST'
  write ( *, '(a)' ) '  R8GE_TO_R8PO converts an R8GE matrix to R8PO format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ge_random ( n, n, seed, a )

  call r8ge_print ( n, n, a, '  The random R8GE matrix A:' )

  call r8ge_to_r8po ( n, a, b )

  call r8po_print ( n, b, '  The R8PO matrix' )

  return
end
subroutine r8po_det_test ( )

!*****************************************************************************80
!
!! R8PO_DET_TEST tests R8PO_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_DET_TEST'
  write ( *, '(a)' ) '  R8PO_DET computes the determinant of'
  write ( *, '(a)' ) '  a symmetric positive definite matrix'
  write ( *, '(a)' ) '  factored by R8PO_FA,.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  a(1:n,1:n) = 0.0D+00
  do i = 1, n
    a(i,i) = 2.0D+00
  end do
  do i = 1, n - 1
    a(i,i+1) = -1.0D+00
  end do
  do i = 2, n
    a(i-1,i) = -1.0D+00
  end do

  call r8po_print ( n, a, '  Matrix A:' )
!
!  Compute the Cholesky factor R.
!
  call r8po_fa ( n, a, r )
!
!  Compute the determinant of A.
!
  call r8po_det ( n, r, det )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant of A = ', det

  return
end
subroutine r8po_dif2_test ( )

!*****************************************************************************80
!
!! R8PO_DIF2_TEST tests R8PO_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_DIF2_TEST'
  write ( *, '(a)' ) '  R8PO_DIF2 sets the second difference matrix in R8PO format.'

  call r8po_dif2 ( n, a )

  call r8po_print ( n, a, '  The R8PO DIF2 matrix:' )
 
  return
end
subroutine r8po_fa_test ( )

!*****************************************************************************80
!
!! R8PO_FA_TEST tests R8PO_FA.
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

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) rtr(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_FA_TEST'
  write ( *, '(a)' ) '  R8PO_FA factors a positive definite symmetric'
  write ( *, '(a)' ) '  linear system,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do

  call r8po_print ( n, a, '  The matrix A:' )
!
!  Factor the matrix.
!
  call r8po_fa ( n, a, r )

  call r8ut_print ( n, n, r, '  The factor R (an R8UT matrix):' )
!
!  Compute the product R' * R.
!
  rtr(1:n,1:n) = matmul ( transpose ( r(1:n,1:n) ), r(1:n,1:n) )

  call r8ge_print ( n, n, rtr, '  The product R'' * R:' )

  return
end
subroutine r8po_indicator_test ( )

!*****************************************************************************80
!
!! R8PO_INDICATOR_TEST tests R8PO_INDICATOR.
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

  real ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8PO_INDICATOR sets up an R8PO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8po_indicator ( n, a )

  call r8po_print ( n, a, '  The R8PO indicator matrix:' )
 
  return
end
subroutine r8po_inverse_test ( )

!*****************************************************************************80
!
!! R8PO_INVERSE_TEST tests R8PO_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_INVERSE_TEST'
  write ( *, '(a)' ) '  R8PO_INVERSE computes the inverse of'
  write ( *, '(a)' ) '  a symmetric positive definite matrix'
  write ( *, '(a)' ) '  factored by R8PO_FA,'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do

  call r8po_print ( n, a, '  Matrix A:' )
!
!  Factor the matrix.
!
  call r8po_fa ( n, a, r )
!
!  Compute the inverse.
!
  call r8po_inverse ( n, r, b )

  call r8po_print ( n, b, '  Inverse matrix B:' )
!
!  Check.
!
  call r8po_mm ( n, a, b, c )

  call r8po_print ( n, c, '  Product A * B:' )

  return
end
subroutine r8po_ml_test ( )

!*****************************************************************************80
!
!! R8PO_ML_TEST tests R8PO_ML.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_ML_TEST'
  write ( *, '(a)' ) '  R8PO_ML can compute A*x for an R8PO matrix A'
  write ( *, '(a)' ) '  even after it has been factored by R8PO_FA.'

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8po_mv ( n, a, x, b )
!
!  Factor the matrix.
!
  call r8po_fa ( n, a, r )
!
!  Solve the linear system.
!
  x(1:n) = 0.0D+00
  call r8po_sl ( n, r, b, x )
 
  call r8vec_print ( n, x, '  Solution:' )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side, using the factored matrix.
!
  call r8po_ml ( n, r, x, b )
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  call r8po_sl ( n, r, b, x )
 
  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine r8po_mm_test ( )

!*****************************************************************************80
!
!! R8PO_MM_TEST tests R8PO_MM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8PO_MM_TEST'
  write ( *, '(a)' ) '  R8PO_MM computes the product of two'
  write ( *, '(a)' ) '  symmetric positive definite matrices.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
!
!  Set (the upper half of) matrix A.
!
  call r8po_zeros ( n, a )

  do i = 1, n
    a(i,i) = 2.0D+00
  end do

  do i = 1, n - 1
    a(i,i+1) = -1.0D+00
  end do

  call r8po_print ( n, a, '  Matrix A:' )
!
!  Set (the upper half of) matrix B.
!
  call r8po_zeros ( n, b )

  do i = 1, n
    b(i,i) = real ( i + i - 1, kind = 8 )
  end do

  do i = 1, n - 1
    b(i,i+1) = real ( i + i + 1 - 1, kind = 8 )
  end do

  call r8po_print ( n, b, '  Matrix B:' )
!
!  Compute the product.
!
  call r8po_mm ( n, a, b, c )

  call r8po_print ( n, c, '  Product matrix C = A * B:' )

  return
end
subroutine r8po_mv_test ( )

!*****************************************************************************80
!
!! R8PO_MV_TEST tests R8PO_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8PO_MV_TEST'
  write ( *, '(a)' ) '  R8PO_MV computes the product of a'
  write ( *, '(a)' ) '  symmetric positive definite matrix times a vector.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
!
!  Set (the upper half of) matrix A.
!
  call r8po_zeros ( n, a )

  do i = 1, n
    a(i,i) = 2.0D+00
  end do

  do i = 1, n - 1
    a(i,i+1) = -1.0D+00
  end do

  call r8po_print ( n, a, '  Matrix A:' )
!
!  Set the vector V.
!
  call r8vec_indicator1 ( n, v )

  call r8vec_print ( n, v, '  Vector V:' )
!
!  Compute the product.
!
  call r8po_mv ( n, a, v, w )

  call r8vec_print ( n, w, '  Product w = A * v:' )

  return
end
subroutine r8po_print_test ( )

!*****************************************************************************80
!
!! R8PO_PRINT_TEST tests R8PO_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_PRINT_TEST'
  write ( *, '(a)' ) '  R8PO_PRINT prints an R8PO matrix.'

  call r8po_indicator ( n, a )

  call r8po_print ( n, a, '  The R8PO matrix:' )
 
  return
end
subroutine r8po_print_some_test ( )

!*****************************************************************************80
!
!! R8PO_PRINT_SOME_TEST tests R8PO_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8PO_PRINT_SOME prints some of an R8PO matrix.'

  call r8po_indicator ( n, a )

  call r8po_print_some ( n, a, 2, 5, 6, 7, '  Rows 2-6, Columns 5-7:' )
 
  return
end
subroutine r8po_random_test ( )

!*****************************************************************************80
!
!! R8PO_RANDOM_TEST tests R8PO_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_RANDOM_TEST'
  write ( *, '(a)' ) '  R8PO_RANDOM computes a random positive definite'
  write ( *, '(a)' ) '  symmetric matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8po_random ( n, seed, a )

  call r8po_print ( n, a, '  The random R8PO matrix:' )
 
  return
end
subroutine r8po_sl_test ( )

!*****************************************************************************80
!
!! R8PO_SL_TEST tests R8PO_SL.
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

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_SL_TEST'
  write ( *, '(a)' ) '  R8PO_SL solves a linear system with R8PO matrix'
  write ( *, '(a)' ) '  after it has been factored buy R8PO_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set (the upper half of) matrix A.
!
  call r8po_zeros ( n, a )

  do i = 1, n
    a(i,i) = 2.0D+00
  end do

  do i = 1, n - 1
    a(i,i+1) = -1.0D+00
  end do

  call r8po_print ( n, a, '  Matrix A:' )
!
!  Factor the matrix.
!
  call r8po_fa ( n, a, r )
!
!  Set the right hand side.
!
  b(1:n-1) = 0.0D+00
  b(n) = real ( n + 1, kind = 8 )
  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Solve the linear system.
!
  call r8po_sl ( n, r, b, x )
 
  call r8vec_print ( n, x, '  Solution vector X:' )

  return
end
subroutine r8po_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8PO_TO_R8GE_TEST tests R8PO_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8PO_TO_R8GE converts an R8PO matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8po_random ( n, seed, a )

  call r8po_print ( n, a, '  The random R8PO matrix A:' )

  call r8ge_print ( n, n, a, '  The random R8PO matrix (printed by R8GE_PRINT):' )

  call r8po_to_r8ge ( n, a, b )

  call r8ge_print ( n, n, b, '  The random R8GE matrix (printed by R8GE_PRINT):' )

  return
end
subroutine r8po_zeros_test ( )

!*****************************************************************************80
!
!! R8PO_ZEROS_TEST tests R8PO_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PO_ZEROS_TEST'
  write ( *, '(a)' ) '  R8PO_ZEROS zeros out a matrix in R8PO format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8po_zeros ( n, a )

  call r8po_print ( n, a, '  The zero R8PO matrix:' )
 
  return
end
