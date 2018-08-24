program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_MAT_PRB.
!
!  Discussion:
!
!    TEST_MAT_PRB tests the TEST_MAT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MAT library.'
!
!  Utilities.
!
  call bvec_next_grlex_test ( ) 
  call legendre_zeros_test ( )
  call mertens_test ( )
  call moebius_test ( )
  call r8mat_is_eigen_left_test ( )
  call r8mat_is_eigen_right_test ( )
  call r8mat_is_llt_test ( )
  call r8mat_is_null_left_test ( )
  call r8mat_is_null_right_test ( )
  call r8mat_is_solution_test ( )
  call r8mat_norm_fro_test ( )
!
!  Important things.
!
  call test_analyze ( )
  call test_condition ( )
  call test_determinant ( )
  call test_eigen_left ( )
  call test_eigen_right ( )
  call test_inverse ( )
  call test_llt ( )
  call test_null_left ( )
  call test_null_right ( )
  call test_plu ( )
  call test_solution ( )
  call test_type ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bvec_next_grlex_test ( )

!*****************************************************************************80
!
!! BVEC_NEXT_GRLEX_TEST tests BVEC_NEXT_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
 
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BVEC_NEXT_GRLEX_TEST'
  write ( *, '(a)' ) '  BVEC_NEXT_GRLEX computes binary vectors in GRLEX order.'
  write ( *, '(a)' ) ''

  b(1:n) = 0

  do i = 0, 16
    write ( *, '(2x,i2,a)', advance = 'no' ) i, ':  '
    do j = 1, n
      write ( *, '(i1)', advance = 'no' ) b(j)
    end do
    write ( *, '(a)' ) ''
    call bvec_next_grlex ( n, b )
  end do

  return
end
subroutine legendre_zeros_test ( )

!*****************************************************************************80
!
!! LEGENDRE_ZEROS_TEST tests LEGENDRE_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: l(:)
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LEGENDRE_ZEROS_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_ZEROS computes the zeros of the N-th Legendre'
  write ( *, '(a)' ) '  polynomial.'

  do n = 1, 7
    allocate ( l(1:n) )
    call legendre_zeros ( n, l )
    call r8vec_print ( n, l, '  Legendre zeros' )
    deallocate ( l )
  end do

  return
end
subroutine mertens_test ( )

!*****************************************************************************80
!
!! MERTENS_TEST tests MERTENS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) mertens
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MERTENS_TEST'
  write ( *, '(a)' ) '  MERTENS computes the Mertens function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MERTENS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call mertens_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    c2 = mertens ( n )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine moebius_test ( )

!*****************************************************************************80
!
!! MOEBIUS_TEST tests MOEBIUS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOEBIUS_TEST'
  write ( *, '(a)' ) '  MOEBIUS computes the Moebius function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MOEBIUS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call moebius_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call moebius ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine r8mat_is_eigen_left_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_EIGEN_LEFT_TEST tests R8MAT_IS_EIGEN_LEFT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: k = 4
!
!  This is the CARRY ( 4.0, 4 ) matrix.
!
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
   0.13671875D+00,   0.05859375D+00,   0.01953125D+00,   0.00390625D+00, &
   0.60546875D+00,   0.52734375D+00,   0.39453125D+00,   0.25390625D+00, &
   0.25390625D+00,   0.39453125D+00,   0.52734375D+00,   0.60546875D+00, &
   0.00390625D+00,   0.01953125D+00,   0.05859375D+00,   0.13671875D+00 /), &
   (/ n, n /) )
  real ( kind = 8 ), dimension ( n ) :: lam = (/ &
     1.000000000000000D+00, &
     0.250000000000000D+00, &
     0.062500000000000D+00, &
     0.015625000000000D+00 /)
  real ( kind = 8 ) value
  real ( kind = 8 ), dimension ( n, k ) :: x = reshape ( (/ &
       1.0D+00, 11.0D+00, 11.0D+00,  1.0D+00, &
       1.0D+00,  3.0D+00, -3.0D+00, -1.0D+00, &
       1.0D+00, -1.0D+00, -1.0D+00,  1.0D+00, &
       1.0D+00, -3.0D+00,  3.0D+00, -1.0D+00 /), &
    (/ n, k /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_EIGEN_LEFT_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_EIGEN_LEFT tests the error in the left eigensystem'
  write ( *, '(a)' ) '    A'' * X - X * LAMBDA = 0'

  call r8mat_print ( n, n, a, '  Matrix A:' )
  call r8mat_print ( n, k, x, '  Eigenmatrix X:' )
  call r8vec_print ( n, lam, '  Eigenvalues LAM:' )

  call r8mat_is_eigen_left ( n, k, a, x, lam, value )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Frobenius norm of A''*X-X*LAMBDA is ', value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_EIGEN_LEFT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine r8mat_is_eigen_right_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_EIGEN_RIGHT_TEST tests R8MAT_IS_EIGEN_RIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: k = 4
!
!  This is the CARRY ( 4.0, 4 ) matrix.
!
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
   0.13671875D+00,   0.05859375D+00,   0.01953125D+00,   0.00390625D+00, &
   0.60546875D+00,   0.52734375D+00,   0.39453125D+00,   0.25390625D+00, &
   0.25390625D+00,   0.39453125D+00,   0.52734375D+00,   0.60546875D+00, &
   0.00390625D+00,   0.01953125D+00,   0.05859375D+00,   0.13671875D+00 /), &
   (/ n, n /) )
  real ( kind = 8 ), dimension ( n ) :: lam = (/ &
     1.000000000000000D+00, &
     0.250000000000000D+00, &
     0.062500000000000D+00, &
     0.015625000000000D+00 /)
  real ( kind = 8 ) value
  real ( kind = 8 ), dimension ( n, k ) :: x = reshape ( (/ &
       1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
       6.0D+00,  2.0D+00, -2.0D+00, -6.0D+00, &
      11.0D+00, -1.0D+00, -1.0D+00, 11.0D+00, &
       6.0D+00, -2.0D+00,  2.0D+00, -6.0D+00 /), &
    (/ n, k /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_EIGEN_RIGHT_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_EIGEN_RIGHT tests the error in the right eigensystem'
  write ( *, '(a)' ) '    A * X - X * LAMBDA = 0'

  call r8mat_print ( n, n, a, '  Matrix A:' )
  call r8mat_print ( n, k, x, '  Eigenmatrix X:' )
  call r8vec_print ( n, lam, '  Eigenvalues LAM:' )

  call r8mat_is_eigen_right ( n, k, a, x, lam, value )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Frobenius norm of A*X-X*LAMBDA is ', value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_EIGEN_RIGHT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine r8mat_is_llt_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_LLT_TEST tests R8MAT_IS_LLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), dimension ( 4, 4 ) :: a = reshape ( (/ &
    2.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 2.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /), (/ 4, 4 /) )
  real ( kind = 8 ) enorm
  real ( kind = 8 ), dimension ( 4, 4 ) :: l = reshape ( (/ &
    1.414213562373095D+00, 0.707106781186547D+00, &
    0.0D+00,               0.0D+00,               &
    0.0D+00,               1.224744871391589D+00, &
    0.816496580927726D+00, 0.0D+00,               &
    0.0D+00,               0.0D+00,               &
    1.154700538379251D+00, 0.866025403784439D+00, &
    0.0D+00,               0.0D+00,               &
    0.0D+00,               1.118033988749895D+00 /), (/ 4, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_LLT_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_LLT tests the error in a lower triangular'
  write ( *, '(a)' ) '  Cholesky factorization A = L * L'' by looking at'
  write ( *, '(a)' ) '    A - L * L'''

  call r8mat_print ( m, m, a, '  Matrix A:' )
  call r8mat_print ( m, n, l, '  Factor L:' )

  call r8mat_is_llt ( m, n, a, l, enorm )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Frobenius norm of A-L*L'' is ', enorm

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_LLT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine r8mat_is_null_left_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_NULL_LEFT_TEST tests R8MAT_IS_NULL_LEFT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
    1.0D+00, 4.0D+00, 7.0D+00, &
    2.0D+00, 5.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ) enorm
  real ( kind = 8 ), dimension ( m ) :: x = (/ 1.0D+00, -2.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_NULL_LEFT_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_NULL_LEFT tests whether the M vector X'
  write ( *, '(a)' ) '  is a left null vector of A, that is, x''*A=0.'

  call r8mat_print ( m, n, a, '  Matrix A:' )
  call r8vec_print ( m, x, '  Vector X:' )

  call r8mat_is_null_left ( m, n, a, x, enorm )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Frobenius norm of X''*A is ', enorm

  return
end
subroutine r8mat_is_null_right_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_NULL_RIGHT_TEST tests R8MAT_IS_NULL_RIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
    1.0D+00, 4.0D+00, 7.0D+00, &
    2.0D+00, 5.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ) enorm
  real ( kind = 8 ), dimension ( n ) :: x = (/ 1.0D+00, -2.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_NULL_RIGHT_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_NULL_RIGHT tests whether the N vector X'
  write ( *, '(a)' ) '  is a right null vector of A, that is, A*X=0.'

  call r8mat_print ( m, n, a, '  Matrix A:' )
  call r8vec_print ( n, x, '  Vector X:' )

  call r8mat_is_null_right ( m, n, a, x, enorm )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Frobenius norm of A*X is ', enorm

  return
end
subroutine r8mat_is_solution_test ( )

!*****************************************************************************80
!
!! R8MAT_IS_SOLUTION_TEST tests R8MAT_IS_SOLUTION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:,:)
  real ( kind = 8 ) enorm
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8MAT_IS_SOLUTION_TEST:'
  write ( *, '(a)' ) '  R8MAT_IS_SOLUTION tests whether X is the solution of'
  write ( *, '(a)' ) '  A*X=B by computing the Frobenius norm of the residual.'
!
!  Get random shapes.
!
  i4_lo = 1
  i4_hi = 10
  seed = 123456789
  m = i4_uniform_ab ( i4_lo, i4_hi, seed )
  n = i4_uniform_ab ( i4_lo, i4_hi, seed )
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
!
!  Get a random A.
!
  allocate ( a(1:m,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  call r8mat_uniform_ab ( m, n, r8_lo, r8_hi, seed, a )
!
!  Get a random X.
!
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  call r8mat_uniform_ab ( n, k, r8_lo, r8_hi, seed, x )
!
!  Compute B = A * X.
!
  allocate ( b(1:m,1:k) )
  b = matmul ( a, x )
!
!  Compute || A*X-B||
!
  call r8mat_is_solution ( m, n, k, a, x, b, enorm )
  
  write ( *, '(a)' ) ''
  write ( *, '(a,i2,a,i2)' ) '  A is ', m, ' by ', n
  write ( *, '(a,i2,a,i2)' ) '  X is ', n, ' by ', k
  write ( *, '(a,i2,a,i2)' ) '  B is ', m, ' by ', k
  write ( *, '(a,g14.6)' ) '  Frobenius error in A*X-B is ', enorm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )

  return
end
subroutine r8mat_norm_fro_test ( )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO_TEST tests R8MAT_NORM_FRO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8MAT_NORM_FRO_TEST'
  write ( *, '(a)' ) '  R8MAT_NORM_FRO computes a Frobenius norm of an R8MAT.'

  t1 = 0.0D+00
  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k, kind = 8 )
      t1 = t1 + real ( k * k, kind = 8 )
    end do
  end do

  t1 = sqrt ( t1 )

  call r8mat_print ( m, n, a, '  A:' )

  t2 = r8mat_norm_fro ( m, n, a )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Expected norm = ', t1
  write ( *, '(a,g14.6)' ) '  Computed norm = ', t2

  return
end
subroutine test_analyze ( )

!*****************************************************************************80
!
!! TEST_ANALYZE tests the matrix analysis functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 20 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ANALYZE'
  write ( *, '(a)' ) '  Analyze a matrix.'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  m = 3
  n = 3
  allocate ( a(1:m,1:n) )
  call a123 ( a )
  write ( *, '(a)' ) ''
  write ( *, '(2x,a)' ) trim ( title )
  write ( *, '(a)' ) ''
  call r8mat_analyze ( m, n, a )
  deallocate ( a )

  return
end
subroutine test_condition ( )

!*****************************************************************************80
!
!! TEST_CONDITION tests the condition number computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable :: b(:,:)
  real ( kind = 8 ) b_norm
  real ( kind = 8 ) beta
  real ( kind = 8 ) cond1
  real ( kind = 8 ) cond2
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_l1
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CONDITION'
  write ( *, '(a)' ) '  Compute the L1 condition number of an example of each'
  write ( *, '(a)' ) '  test matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    N      COND            COND'
  write ( *, '(a)' ) ' '
!
!  AEGERTER
!
  title = 'AEGERTER'
  n = 5
  call aegerter_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call aegerter ( n, a )
  call aegerter_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BAB
!
  title = 'BAB'
  n = 5
  seed = 123456789
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call bab_condition ( n, alpha, beta, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call bab ( n, alpha, beta, a )
  call bab_inverse ( n, alpha, beta, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BAUER
!
  title = 'BAUER'
  n = 6
  call bauer_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call bauer ( a )
  call bauer_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BIS
!
  title = 'BIS'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call bis_condition ( alpha, beta, n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call bis ( alpha, beta, n, n, a )
  call bis_inverse ( alpha, beta, n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BIW
!
  title = 'BIW'
  n = 5
  call biw_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call biw ( n, a )
  call biw_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BODEWIG
!
  title = 'BODEWIG'
  n = 4
  call bodewig_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call bodewig ( a )
  call bodewig_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  BOOTHROYD
!
  title = 'BOOTHROYD'
  n = 5
  call boothroyd_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call boothroyd ( n, a )
  call boothroyd_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  COMBIN
!
  title = 'COMBIN'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call combin_condition ( alpha, beta, n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  COMPANION
!
  title = 'COMPANION'
  n = 5
  seed = 123456789
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call companion_condition ( n, x, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call companion ( n, x, a )
  call companion_inverse ( n, x, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  CONEX1
!
  title = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex1_condition ( alpha, cond1 )
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  CONEX2
!
  title = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex2_condition ( alpha, cond1 )

  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  CONEX3
!
  title = 'CONEX3'
  n = 5
  call conex3_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  CONEX4
!
  title = 'CONEX4'
  n = 4
  call conex4_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call conex4 ( a )
  call conex4_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB2
!
  title = 'DAUB2'
  n = 4
  call daub2_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub2 ( n, a )
  call daub2_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB4
!
  title = 'DAUB4'
  n = 8
  call daub4_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub4 ( n, a )
  call daub4_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB6
!
  title = 'DAUB6'
  n = 12
  call daub6_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub6 ( n, a )
  call daub6_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB8
!
  title = 'DAUB8'
  n = 16
  call daub8_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub8 ( n, a )
  call daub8_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB10
!
  title = 'DAUB10'
  n = 20
  call daub10_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub10 ( n, a )
  call daub10_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DAUB12
!
  title = 'DAUB12'
  n = 24
  call daub12_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call daub12 ( n, a )
  call daub12_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DIAGONAL
!
  title = 'DIAGONAL'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  allocate ( x(1:n) )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call diagonal_condition ( n, x, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call diagonal ( n, n, x, a )
  call diagonal_inverse ( n, x, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  DIF2
!
  title = 'DIF2'
  n = 5
  call dif2_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call dif2 ( n, n, a )
  call dif2_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  DOWNSHIFT
!
  title = 'DOWNSHIFT'
  n = 5
  call downshift_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call downshift ( n, a )
  call downshift_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  EXCHANGE
!
  title = 'EXCHANGE'
  n = 5
  call exchange_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call exchange ( n, n, a )
  call exchange_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  FIBONACCI2
!
  title = 'FIBONACCI2'
  n = 5
  call fibonacci2_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call fibonacci2 ( n, a )
  call fibonacci2_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  GFPP
!
  title = 'GFPP'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call gfpp_condition ( n, alpha, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call gfpp ( n, alpha, a )
  call gfpp_inverse ( n, alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  GIVENS
!
  title = 'GIVENS'
  n = 5
  call givens_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call givens ( n, n, a )
  call givens_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  HANKEL_N
!
  title = 'HANKEL_N'
  n = 5
  call hankel_n_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call hankel_n ( n, a )
  call hankel_n_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  HARMAN
!
  title = 'HARMAN'
  n = 8
  call harman_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call harman ( a )
  call harman_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  HARTLEY
!
  title = 'HARTLEY'
  n = 5
  call hartley_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call hartley ( n, a )
  call hartley_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  IDENTITY
!
  title = 'IDENTITY'
  n = 5
  call identity_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call identity ( n, n, a )
  call identity_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  ILL3
!
  title = 'ILL3'
  n = 3
  call ill3_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call ill3 ( a )
  call ill3_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  JORDAN
!
  title = 'JORDAN'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call jordan_condition ( n, alpha, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call jordan ( n, n, alpha, a )
  call jordan_inverse ( n, alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  KERSHAW
!
  title = 'KERSHAW'
  n = 4
  call kershaw_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call kershaw ( a )
  call kershaw_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  LIETZKE
!
  title = 'LIETZKE'
  n = 5
  call lietzke_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call lietzke ( n, a )
  call lietzke_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  MAXIJ
!
  title = 'MAXIJ'
  n = 5
  call maxij_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call maxij ( n, n, a )
  call maxij_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  MINIJ
!
  title = 'MINIJ'
  n = 5
  call minij_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call minij ( n, n, a )
  call minij_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  ORTH_SYMM
!
  title = 'ORTH_SYMM'
  n = 5
  call orth_symm_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call orth_symm ( n, a )
  call orth_symm_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  OTO
!
  title = 'OTO'
  n = 5
  call oto_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call oto ( n, n, a )
  call oto_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  PASCAL1
!
  title = 'PASCAL1'
  n = 5
  call pascal1_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call pascal1 ( n, a )
  call pascal1_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  PASCAL3
!
  title = 'PASCAL3'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pascal3_condition ( n, alpha, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call pascal3 ( n, alpha, a )
  call pascal3_inverse ( n, alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  PEI
!
  title = 'PEI'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pei_condition ( alpha, n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call pei ( alpha, n, a )
  call pei_inverse ( alpha, n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  RODMAN
!
  title = 'RODMAN'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call rodman_condition ( n, alpha, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call rodman ( n, n, alpha, a )
  call rodman_inverse ( n, alpha, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  RUTIS1
!
  title = 'RUTIS1'
  n = 4
  call rutis1_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call rutis1 ( a )
  call rutis1_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  RUTIS2
!
  title = 'RUTIS2'
  n = 4
  call rutis2_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call rutis2 ( a )
  call rutis2_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  RUTIS3
!
  title = 'RUTIS3'
  n = 4
  call rutis3_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call rutis3 ( a )
  call rutis3_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  RUTIS5
!
  title = 'RUTIS5'
  n = 4
  call rutis5_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call rutis5 ( a )
  call rutis5_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  SUMMATION
!
  title = 'SUMMATION'
  n = 5
  call summation_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call summation ( n, n, a )
  call summation_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  SWEET1
!
  title = 'SWEET1'
  n = 6
  call sweet1_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call sweet1 ( a )
  call sweet1_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  SWEET2
!
  title = 'SWEET2'
  n = 6
  call sweet2_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call sweet2 ( a )
  call sweet2_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  SWEET3
!
  title = 'SWEET3'
  n = 6
  call sweet3_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call sweet3 ( a )
  call sweet3_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  SWEET4
!
  title = 'SWEET4'
  n = 13
  call sweet4_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call sweet4 ( a )
  call sweet4_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  TRI_UPPER
!
  title = 'TRI_UPPER'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call tri_upper_condition ( alpha, n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call tri_upper ( alpha, n, a )
  call tri_upper_inverse ( alpha, n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  UPSHIFT
!
  title = 'UPSHIFT'
  n = 5
  call upshift_condition ( n, cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call upshift ( n, a )
  call upshift_inverse ( n, b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  WILK03
!
  title = 'WILK03'
  n = 3
  call wilk03_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call wilk03 ( a )
  call wilk03_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  WILK04
!
  title = 'WILK04'
  n = 4
  call wilk04_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call wilk04 ( a )
  call wilk04_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  WILK05
!
  title = 'WILK05'
  n = 5
  call wilk05_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call wilk05 ( a )
  call wilk05_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )
!
!  WILSON
!
  title = 'WILSON'
  n = 4
  call wilson_condition ( cond1 )

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  call wilson ( a )
  call wilson_inverse ( b )
  a_norm = r8mat_norm_l1 ( n, n, a )
  b_norm = r8mat_norm_l1 ( n, n, b )
  cond2 = a_norm * b_norm

  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) title, n, cond1, cond2
  deallocate ( a )
  deallocate ( b )

  return
end
subroutine test_determinant ( )

!*****************************************************************************80
!
!! TEST_DETERMINANT tests the determinant computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) col_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) d4
  real ( kind = 8 ) d5
  real ( kind = 8 ) da
  real ( kind = 8 ) determ1
  real ( kind = 8 ) determ2
  real ( kind = 8 ) di
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  real ( kind = 8 ) prob
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  integer ( kind = 4 ) x_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  integer ( kind = 4 ) y_n
  real ( kind = 8 ) y_sum
  real ( kind = 8 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DETERMINANT'
  write ( *, '(a)' ) '  Compute the determinants of an example of each'
  write ( *, '(a)' ) '  test matrix.  Compare with the determinant routine,'
  write ( *, '(a)' ) '  if available.  Print the matrix Frobenius norm'
  write ( *, '(a)' ) '  for an estimate of magnitude.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    N      ' // &
                     'Determ          Determ         ||A||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  n = 3
  allocate ( a(1:n,1:n) )
  call a123 ( a )
  call a123_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  AEGERTER
!
  title = 'AEGERTER'
  n = 5
  allocate ( a(1:n,1:n) )
  call aegerter ( n, a )
  call aegerter_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ANTICIRCULANT
!
  title = 'ANTICIRCULANT'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTICIRCULANT
!
  title = 'ANTICIRCULANT'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTICIRCULANT
!
  title = 'ANTICIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTIHADAMARD
!
  title = 'ANTIHADAMARD'
  n = 5
  allocate ( a(1:n,1:n) )
  call antihadamard ( n, a )
  call antihadamard_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ANTISYMM_RANDOM
!
  title = 'ANTISYMM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call antisymm_random ( n, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  ANTISYMM_RANDOM
!
  title = 'ANTISYMM_RANDOM'
  n = 6
  allocate ( a(1:n,1:n) )
  key = 123456789
  call antisymm_random ( n, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  BAB
!
  title = 'BAB'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  allocate ( a(1:n,1:n) )
  call bab ( n, alpha, beta, a )
  call bab_determinant ( n, alpha, beta, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BAUER
!
  title = 'BAUER'
  n = 6
  allocate ( a(1:n,1:n) )
  call bauer ( a )
  call bauer_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BERNSTEIN
!
  title = 'BERNSTEIN'
  n = 5
  allocate ( a(1:n,1:n) )
  call bernstein ( n, a )
  call bernstein_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BIMARKOV_RANDOM
!
  title = 'BIMARKOV_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call bimarkov_random ( n, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  BIS
!
  title = 'BIS'
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  allocate ( a(1:n,1:n) )
  call bis ( alpha, beta, n, n, a )
  call bis_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BIW
!
  title = 'BIW'
  n = 5
  allocate ( a(1:n,1:n) )
  call biw ( n, a )
  call biw_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BODEWIG
!
  title = 'BODEWIG'
  n = 4
  allocate ( a(1:n,1:n) )
  call bodewig ( a )
  call bodewig_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BOOTHROYD
!
  title = 'BOOTHROYD'
  n = 5
  allocate ( a(1:n,1:n) )
  call boothroyd ( n, a )
  call boothroyd_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BORDERBAND
!
  title = 'BORDERBAND'
  n = 5
  allocate ( a(1:n,1:n) )
  call borderband ( n, a )
  call borderband_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CARRY
!
  title = 'CARRY'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  i4_lo = 2
  i4_hi = 20
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
  call carry ( n, k, a )
  call carry_determinant ( n, k, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CAUCHY
!
  title = 'CAUCHY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  call cauchy ( n, x, y, a )
  call cauchy_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  CHEBY_DIFF1
!
  title = 'CHEBY_DIFF1'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_diff1 ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_DIFF1
!
  title = 'CHEBY_DIFF1'
  n = 6
  allocate ( a(1:n,1:n) )
  call cheby_diff1 ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_T
!
  title = 'CHEBY_T'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_t ( n, a )
  call cheby_t_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_U
!
  title = 'CHEBY_U'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_u ( n, a )
  call cheby_u_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_VAN1
!
  title = 'CHEBY_VAN1'
  n = 5
  x_lo = -1.0D+00
  x_hi = +1.0D+00
  allocate ( x(1:n) )
  call r8vec_linspace ( n, x_lo, x_hi, x )
  allocate ( a(1:n,1:n) )
  call cheby_van1 ( n, x_lo, x_hi, n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  CHEBY_VAN2
!
  do n = 2, 10
    title = 'CHEBY_VAN2'
    allocate ( a(1:n,1:n) )
    call cheby_van2 ( n, a )
    call cheby_van2_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  CHEBY_VAN3
!
  title = 'CHEBY_VAN3'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_van3 ( n, a )
  call cheby_van3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHOW
!
  title = 'CHOW'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call chow ( alpha, beta, n, n, a )
  call chow_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT
!
  title = 'CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call circulant ( n, n, x, a )
  call circulant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  CIRCULANT2
!
  title = 'CIRCULANT2'
  n = 3
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT2
!
  title = 'CIRCULANT2'
  n = 4
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT2
!
  title = 'CIRCULANT2'
  n = 5
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT1
!
  title = 'CLEMENT1'
  n = 5
  allocate ( a(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT1
!
  title = 'CLEMENT1'
  n = 6
  allocate ( a(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT2
!
  title = 'CLEMENT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, y )
  call clement2 ( n, x, y, a )
  call clement2_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  CLEMENT2
!
  title = 'CLEMENT2'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, y )
  call clement2 ( n, x, y, a )
  call clement2_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  COMBIN
!
  title = 'COMBIN'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call combin ( alpha, beta, n, a )
  call combin_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  COMPANION
!
  title = 'COMPANION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call companion ( n, x, a )
  call companion_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  COMPLEX_I
!
  title = 'COMPLEX_I'
  n = 2
  allocate ( a(1:n,1:n) )
  call complex_i ( a )
  call complex_i_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX1
!
  title = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex1 ( alpha, a )
  call conex1_determinant ( alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX2
!
  title = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex2 ( alpha, a )
  call conex2_determinant ( alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX3
!
  title = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX4
!
  title = 'CONEX4'
  n = 4
  allocate ( a(1:n,1:n) )
  call conex4 ( a )
  call conex4_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONFERENCE
!
  title = 'CONFERENCE'
  n = 6
  allocate ( a(1:n,1:n) )
  call conference ( n, a )
  call conference_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CREATION
!
  title = 'CREATION'
  n = 5
  allocate ( a(1:n,1:n) )
  call creation ( n, n, a )
  call creation_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB2
!
  title = 'DAUB2'
  n = 4
  allocate ( a(1:n,1:n) )
  call daub2 ( n, a )
  call daub2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB4
!
  title = 'DAUB4'
  n = 8
  allocate ( a(1:n,1:n) )
  call daub4 ( n, a )
  call daub4_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB6
!
  title = 'DAUB6'
  n = 12
  allocate ( a(1:n,1:n) )
  call daub6 ( n, a )
  call daub6_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB8
!
  title = 'DAUB8'
  n = 16
  allocate ( a(1:n,1:n) )
  call daub8 ( n, a )
  call daub8_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB10
!
  title = 'DAUB10'
  n = 20
  allocate ( a(1:n,1:n) )
  call daub10 ( n, a )
  call daub10_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB12
!
  title = 'DAUB12'
  n = 24
  allocate ( a(1:n,1:n) )
  call daub12 ( n, a )
  call daub12_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIAGONAL
!
  title = 'DIAGONAL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call diagonal ( n, n, x, a )
  call diagonal_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  DIF1
!
  title = 'DIF1'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif1 ( n, n, a )
  call dif1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF1
!
  title = 'DIF1'
  n = 6
  allocate ( a(1:n,1:n) )
  call dif1 ( n, n, a )
  call dif1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF1CYCLIC
!
  title = 'DIF1CYCLIC'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif1cyclic ( n, a )
  call dif1cyclic_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF2
!
  title = 'DIF2'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif2 ( n, n, a )
  call dif2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF2CYCLIC
!
  title = 'DIF2CYCLIC'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif2cyclic ( n, a )
  call dif2cyclic_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DORR
!
  title = 'DORR'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call dorr ( alpha, n, a )
  call dorr_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DOWNSHIFT
!
  title = 'DOWNSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  call downshift ( n, a )
  call downshift_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EBERLEIN
!
  title = 'EBERLEIN'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call eberlein ( alpha, n, a )
  call eberlein_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EULERIAN
!
  title = 'EULERIAN'
  n = 5
  allocate ( a(1:n,1:n) )
  call eulerian ( n, n, a )
  call eulerian_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EXCHANGE
!
  title = 'EXCHANGE'
  n = 5
  allocate ( a(1:n,1:n) )
  call exchange ( n, n, a )
  call exchange_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI1
!
  title = 'FIBONACCI1'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call fibonacci1 ( n, alpha, beta, a )
  call fibonacci1_determinant ( n, alpha, beta, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI2
!
  title = 'FIBONACCI2'
  n = 5
  allocate ( a(1:n,1:n) )
  call fibonacci2 ( n, a )
  call fibonacci2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI3
!
  title = 'FIBONACCI3'
  n = 5
  allocate ( a(1:n,1:n) )
  call fibonacci3 ( n, a )
  call fibonacci3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIEDLER.
!
  title = 'FIEDLER'
  n = 7
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call fiedler ( n, n, x, a )
  call fiedler_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  FORSYTHE
!
  title = 'FORSYTHE'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call forsythe ( alpha, beta, n, a )
  call forsythe_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FORSYTHE
!
  title = 'FORSYTHE'
  n = 6
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call forsythe ( alpha, beta, n, a )
  call forsythe_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FOURIER_COSINE
!
  title = 'FOURIER_COSINE'
  n = 5
  allocate ( a(1:n,1:n) )
  call fourier_cosine ( n, a )
  call fourier_cosine_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FOURIER_SINE
!
  title = 'FOURIER_SINE'
  n = 5
  allocate ( a(1:n,1:n) )
  call fourier_sine ( n, a )
  call fourier_sine_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FRANK
!
  title = 'FRANK'
  n = 5
  allocate ( a(1:n,1:n) )
  call frank ( n, a )
  call frank_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GEAR
!
  do n = 4, 8
    title = 'GEAR'
    allocate ( a(1:n,1:n) )
    i4_lo = -n
    i4_hi = +n
    seed = 123456789
    ii = i4_uniform_ab ( i4_lo, i4_hi, seed )
    jj = i4_uniform_ab ( i4_lo, i4_hi, seed )
    call gear ( ii, jj, n, a )
    call gear_determinant ( ii, jj, n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  GFPP
!
  title = 'GFPP'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call gfpp ( n, alpha, a )
  call gfpp_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GIVENS.
!
  title = 'GIVENS'
  n = 5
  allocate ( a(1:n,1:n) )
  call givens ( n, n, a )
  call givens_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK316
!
  title = 'GK316'
  n = 5
  allocate ( a(1:n,1:n) )
  call gk316 ( n, a )
  call gk316_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK323
!
  title = 'GK323'
  n = 5
  allocate ( a(1:n,1:n) )
  call gk323 ( n, n, a )
  call gk323_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK324
!
  title = 'GK324'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call gk324 ( n, n, x, a )
  call gk324_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  GRCAR
!
  title = 'GRCAR'
  n = 5
  seed = 123456789
  allocate ( a(1:n,1:n) )
  i4_lo = 1
  i4_hi = n - 1
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
  call grcar ( n, n, k, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  HADAMARD
!
  title = 'HADAMARD'
  n = 5
  allocate ( a(1:n,1:n) )
  call hadamard ( n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  HANKEL
!
  title = 'HANKEL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:2*n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( 2 * n - 1, r8_lo, r8_hi, seed, x )
  call hankel ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  HANKEL_N
!
  title = 'HANKEL_N'
  n = 5
  allocate ( a(1:n,1:n) )
  call hankel_n ( n, a )
  call hankel_n_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HANOWA
!
  title = 'HANOWA'
  n = 6
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call hanowa ( alpha, n, a )
  call hanowa_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HARMAN
!
  title = 'HARMAN'
  n = 8
  allocate ( a(1:n,1:n) )
  call harman ( a )
  call harman_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HARTLEY
!
  title = 'HARTLEY'
  do n = 5, 8
    allocate ( a(1:n,1:n) )
    call hartley ( n, a )
    call hartley_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  HELMERT
!
  title = 'HELMERT'
  n = 5
  allocate ( a(1:n,1:n) )
  call helmert ( n, a )
  call helmert_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HELMERT2
!
  title = 'HELMERT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call helmert2 ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  HERMITE
!
  title = 'HERMITE'
  n = 5
  allocate ( a(1:n,1:n) )
  call hermite ( n, a )
  call hermite_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HERNDON
!
  title = 'HERNDON'
  n = 5
  allocate ( a(1:n,1:n) )
  call herndon ( n, a )
  call herndon_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HILBERT
!
  title = 'HILBERT'
  n = 5
  allocate ( a(1:n,1:n) )
  call hilbert ( n, n, a )
  call hilbert_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HOUSEHOLDER
!
  title = 'HOUSEHOLDER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call householder ( n, x, a )
  call householder_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  IDEM_RANDOM
!
  title = 'IDEM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  i4_lo = 0
  i4_hi = n
  rank = i4_uniform_ab ( i4_lo, i4_hi, seed )
  key = 123456789
  call idem_random ( n, rank, key, a )
  call idem_random_determinant ( n, rank, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IDENTITY
!
  title = 'IDENTITY'
  n = 5
  allocate ( a(1:n,1:n) )
  call identity ( n, n, a )
  call identity_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IJFACT1
!
  title = 'IJFACT1'
  n = 5
  allocate ( a(1:n,1:n) )
  call ijfact1 ( n, a )
  call ijfact1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IJFACT2
!
  title = 'IJFACT2'
  n = 5
  allocate ( a(1:n,1:n) )
  call ijfact2 ( n, a )
  call ijfact2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ILL3
!
  title = 'ILL3'
  n = 3
  allocate ( a(1:n,1:n) )
  call ill3 ( a )
  call ill3_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INTEGRATION
!
  title = 'INTEGRATION'
  n = 6
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call integration ( alpha, n, a )
  call integration_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INVOL
!
  title = 'INVOL'
  n = 5
  allocate ( a(1:n,1:n) )
  call invol ( n, a )
  call invol_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INVOL_RANDOM
!
  title = 'INVOL_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  i4_lo = 0
  i4_hi = n
  rank = i4_uniform_ab ( i4_lo, i4_hi, seed )
  key = 123456789
  call invol_random ( n, rank, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  JACOBI
!
  title = 'JACOBI'
  n = 5
  allocate ( a(1:n,1:n) )
  call jacobi ( n, n, a )
  call jacobi_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  JACOBI
!
  title = 'JACOBI'
  n = 6
  allocate ( a(1:n,1:n) )
  call jacobi ( n, n, a )
  call jacobi_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  JORDAN
!
  title = 'JORDAN'
  n = 6
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call jordan ( n, n, alpha, a )
  call jordan_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KAHAN
!
  title = 'KAHAN'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call kahan ( alpha, n, n, a )
  call kahan_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KERSHAW
!
  title = 'KERSHAW'
  n = 4
  allocate ( a(1:n,1:n) )
  call kershaw ( a )
  call kershaw_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KERSHAWTRI
!
  title = 'KERSHAWTRI'
  n = 5
  x_n = ( n + 1 ) / 2
  allocate ( a(1:n,1:n) )
  allocate ( x(1:x_n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( x_n, r8_lo, r8_hi, seed, x )
  call kershawtri ( n, x, a )
  call kershawtri_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  KMS
!
  title = 'KMS'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call kms ( alpha, n, n, a )
  call kms_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LAGUERRE
!
  title = 'LAGUERRE'
  n = 5
  allocate ( a(1:n,1:n) )
  call laguerre ( n, a )
  call laguerre_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LEGENDRE
!
  title = 'LEGENDRE'
  n = 5
  allocate ( a(1:n,1:n) )
  call legendre ( n, a )
  call legendre_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LEHMER
!
  title = 'LEHMER'
  n = 5
  allocate ( a(1:n,1:n) )
  call lehmer ( n, n, a )
  call lehmer_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LESLIE
!
  title = 'LESLIE'
  n = 4
  allocate ( a(1:n,1:n) )
  b =  0.025D+00
  di = 0.010D+00
  da = 0.100D+00
  call leslie ( b, di, da, a )
  call leslie_determinant ( b, di, da, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LESP
!
  title = 'LESP'
  n = 5
  allocate ( a(1:n,1:n) )
  call lesp ( n, n, a )
  call lesp_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LIETZKE
!
  title = 'LIETZKE'
  n = 5
  allocate ( a(1:n,1:n) )
  call lietzke ( n, a )
  call lietzke_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LIGHTS_OUT
!
  title = 'LIGHTS_OUT'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:row_num*col_num,1:row_num*col_num) )
  call lights_out ( row_num, col_num, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  LINE_ADJ
!
  title = 'LINE_ADJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call line_adj ( n, a )
  call line_adj_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LINE_ADJ
!
  title = 'LINE_ADJ'
  n = 6
  allocate ( a(1:n,1:n) )
  call line_adj ( n, a )
  call line_adj_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LINE_LOOP_ADJ
!
  title = 'LINE_LOOP_ADJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call line_loop_adj ( n, a )
  call line_loop_adj_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LOEWNER
!
  title = 'LOEWNER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( z(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, w )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, z )
  call loewner ( w, x, y, z, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( w )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  LOTKIN
!
  title = 'LOTKIN'
  n = 5
  allocate ( a(1:n,1:n) )
  call lotkin ( n, n, a )
  call lotkin_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MARKOV_RANDOM
!
  title = 'MARKOV_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call markov_random ( n, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  MAXIJ
!
  title = 'MAXIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call maxij ( n, n, a )
  call maxij_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MILNES
!
  title = 'MILNES'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call milnes ( n, n, x, a )
  call milnes_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  MINIJ
!
  title = 'MINIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call minij ( n, n, a )
  call minij_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER1
!
  title = 'MOLER1'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call moler1 ( alpha, n, n, a )
  call moler1_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER2
!
  title = 'MOLER2'
  n = 5
  allocate ( a(1:n,1:n) )
  call moler2 ( a )
  call moler2_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER3
!
  title = 'MOLER3'
  n = 5
  allocate ( a(1:n,1:n) )
  call moler3 ( n, n, a )
  call moler3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER4
!
  title = 'MOLER4'
  n = 4
  allocate ( a(1:n,1:n) )
  call moler4 ( a )
  call moler4_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  NEUMANN
!
  title = 'NEUMANN'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  call neumann ( row_num, col_num, a )
  call neumann_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ONE
!
  title = 'ONE'
  n = 5
  allocate ( a(1:n,1:n) )
  call one ( n, n, a )
  call one_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ORTEGA
!
  title = 'ORTEGA'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v1 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v2 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v3 )
  call ortega ( n, v1, v2, v3, a )
  call ortega_determinant ( n, v1, v2, v3, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
!
!  ORTH_RANDOM
!
  title = 'ORTH_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call orth_random ( n, key, a )
  call orth_random_determinant ( n, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ORTH_SYMM
!
  title = 'ORTH_SYMM'
  n = 5
  allocate ( a(1:n,1:n) )
  call orth_symm ( n, a )
  call orth_symm_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  OTO
!
  title = 'OTO'
  n = 5
  allocate ( a(1:n,1:n) )
  call oto ( n, n, a )
  call oto_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PARTER
!
  title = 'PARTER'
  n = 5
  allocate ( a(1:n,1:n) )
  call parter ( n, n, a )
  call parter_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL1
!
  title = 'PASCAL1'
  n = 5
  allocate ( a(1:n,1:n) )
  call pascal1 ( n, a )
  call pascal1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL2
!
  title = 'PASCAL2'
  n = 5
  allocate ( a(1:n,1:n) )
  call pascal2 ( n, a )
  call pascal2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL3
!
  title = 'PASCAL3'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pascal3 ( n, alpha, a )
  call pascal3_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PDS_RANDOM
!
  title = 'PDS_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call pds_random ( n, key, a )
  call pds_random_determinant ( n, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PEI
!
  title = 'PEI'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pei ( alpha, n, a )
  call pei_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PERMUTATION_RANDOM
!
  title = 'PERMUTATION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call permutation_random ( n, key, a )
  call permutation_random_determinant ( n, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PLU
!
  title = 'PLU'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( l(1:n,1:n) )
  allocate ( p(1:n,1:n) )
  allocate ( pivot(n) )
  allocate ( u(1:n,1:n) )
  seed = 123456789
  do i = 1, n
    i4_lo = i
    i4_hi = n
    pivot(i) = i4_uniform_ab ( i4_lo, i4_hi, seed )
  end do
  call plu ( n, pivot, a )
  call plu_determinant ( n, pivot, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( pivot )
  deallocate ( u )
!
!  POISSON
!
  title = 'POISSON'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  call poisson ( row_num, col_num, a )
  call poisson_determinant ( row_num, col_num, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PROLATE
!
  title = 'PROLATE'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call prolate ( alpha, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  RECTANGLE_ADJ
!
  title = 'RECTANGLE_ADJ'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  call rectangle_adj ( row_num, col_num, n, a )
  call rectangle_adj_determinant ( row_num, col_num, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REDHEFFER
!
  title = 'REDHEFFER'
  n = 5
  allocate ( a(1:n,1:n) )
  call redheffer ( n, a )
  call redheffer_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REF_RANDOM
!
  title = 'REF_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  prob = 0.65D+00
  key = 123456789
  call ref_random ( n, n, prob, key, a )
  call ref_random_determinant ( n, prob, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REF_RANDOM
!
  title = 'REF_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  prob = 0.85D+00
  key = 123456789
  call ref_random ( n, n, prob, key, a )
  call ref_random_determinant ( n, prob, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RIEMANN
!
  title = 'RIEMANN'
  n = 5
  allocate ( a(1:n,1:n) )
  call riemann ( n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  RING_ADJ
!
  do n = 1, 8
    title = 'RING_ADJ'
    allocate ( a(1:n,1:n) )
    call ring_adj ( n, a )
    call ring_adj_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  RIS
!
  title = 'RIS'
  n = 5
  allocate ( a(1:n,1:n) )
  call ris ( n, a )
  call ris_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RODMAN
!
  title = 'RODMAN'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call rodman ( n, n, alpha, a )
  call rodman_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ROSSER1
!
!  Note that while the correct determinant of this matrix is 0,
!  that value is very difficult to calculate correctly.  MATLAB
!  returns det ( A ) = -10611, for instance.
!
  title = 'ROSSER1'
  n = 8
  allocate ( a(1:n,1:n) )
  call rosser1 ( a )
  call rosser1_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ROUTH
!
  title = 'ROUTH'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n ) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call routh ( n, x, a )
  call routh_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  RUTIS1
!
  title = 'RUTIS1'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis1 ( a )
  call rutis1_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS2
!
  title = 'RUTIS2'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis2 ( a )
  call rutis2_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS3
!
  title = 'RUTIS3'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis3 ( a )
  call rutis3_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS4
!
  title = 'RUTIS4'
  n = 5
  allocate ( a(1:n,1:n) )
  call rutis4 ( n, a )
  call rutis4_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS5
!
  title = 'RUTIS5'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis5 ( a )
  call rutis5_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SCHUR_BLOCK
!
  title = 'SCHUR_BLOCK'
  n = 5
  x_n = ( n + 1 ) / 2
  y_n = n / 2
  allocate ( a(1:n,1:n) )
  allocate ( x(1:x_n) )
  allocate ( y(1:y_n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( x_n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( y_n, r8_lo, r8_hi, seed, y )
  call schur_block ( n, x, y, a )
  call schur_block_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  SKEW_CIRCULANT
!
  title = 'SKEW_CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call skew_circulant ( n, n, x, a )
  call skew_circulant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  SPLINE
!
  title = 'SPLINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call spline ( n, x, a )
  call spline_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  STIRLING
!
  title = 'STIRLING'
  n = 5
  allocate ( a(1:n,1:n) )
  call stirling ( n, n, a )
  call stirling_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  STRIPE
!
  title = 'STRIPE'
  n = 5
  allocate ( a(1:n,1:n) )
  call stripe ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SUMMATION
!
  title = 'SUMMATION'
  n = 5
  allocate ( a(1:n,1:n) )
  call summation ( n, n, a )
  call summation_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET1
!
  title = 'SWEET1'
  n = 6
  allocate ( a(1:n,1:n) )
  call sweet1 ( a )
  call sweet1_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET2
!
  title = 'SWEET2'
  n = 6
  allocate ( a(1:n,1:n) )
  call sweet2 ( a )
  call sweet2_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET3
!
  title = 'SWEET3'
  n = 6
  allocate ( a(1:n,1:n) )
  call sweet3 ( a )
  call sweet3_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET4
!
  title = 'SWEET4'
  n = 13
  allocate ( a(1:n,1:n) )
  call sweet4 ( a )
  call sweet4_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SYLVESTER
!
  title = 'SYLVESTER'
  n = 5
  x_n = ( n / 2 )
  y_n = n - ( n / 2 )
  allocate ( a(1:n,1:n) )
  allocate ( x(0:x_n) )
  allocate ( y(0:y_n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( x_n + 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( y_n + 1, r8_lo, r8_hi, seed, y )
  call sylvester ( n, x_n, x, y_n, y, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  SYLVESTER_KAC
!
  title = 'SYLVESTER_KAC'
  n = 5
  allocate ( a(1:n,1:n) )
  call sylvester_kac ( n, a )
  call sylvester_kac_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SYLVESTER_KAC
!
  title = 'SYLVESTER_KAC'
  n = 6
  allocate ( a(1:n,1:n) )
  call sylvester_kac ( n, a )
  call sylvester_kac_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SYMM_RANDOM
!
  title = 'SYMM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( d(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, d )
  key = 123456789
  call symm_random ( n, d, key, a )
  call symm_random_determinant ( n, d, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( d )
!
!  TOEPLITZ
!
  title = 'TOEPLITZ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:2*n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( 2 * n - 1, r8_lo, r8_hi, seed, x )
  call toeplitz ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  TOEPLITZ_5DIAG
!
  title = 'TOEPLITZ_5DIAG'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  d1 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d2 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d3 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d4 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d5 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call toeplitz_5diag ( n, d1, d2, d3, d4, d5, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TOEPLITZ_5S
!
  title = 'TOEPLITZ_5S'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TOEPLITZ_PDS
!
  title = 'TOEPLITZ_PDS'
  m = 3
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:m) )
  allocate ( y(1:m) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( m, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( m, r8_lo, r8_hi, seed, y )
  y_sum = sum ( y(1:m) )
  y(1:m) = y(1:m) / y_sum
  call toeplitz_pds ( m, n, x, y, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  TOURNAMENT_RANDOM
!
  title = 'TOURNAMENT_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  key = 123456789
  call tournament_random ( n, key, a )
  call tournament_random_determinant ( n, key, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRANSITION_RANDOM
!
  title = 'TRANSITION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = key
  call transition_random ( n, key, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TRENCH
!
  title = 'TRENCH'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call trench ( alpha, n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TRI_UPPER
!
  title = 'TRI_UPPER'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call tri_upper ( alpha, n, a )
  call tri_upper_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRIS
!
  title = 'TRIS'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call tris ( n, n, alpha, beta, gamma, a )
  call tris_determinant ( n, alpha, beta, gamma, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRIV
!
  title = 'TRIV'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n) )
  allocate ( z(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, z )
  call triv ( n, x, y, z, a )
  call triv_determinant ( n, x, y, z, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  TRIW
!
  title = 'TRIW'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  i4_lo = 0
  i4_hi = n - 1
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call triw ( alpha, k, n, a )
  call triw_determinant ( alpha, k, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  UPSHIFT
!
  title = 'UPSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  call upshift ( n, a )
  call upshift_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  VAND1
!
  title = 'VAND1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call vand1 ( n, x, a )
  call vand1_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  VAND2
!
  title = 'VAND2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call vand2 ( n, x, a )
  call vand2_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  WATHEN
!
  title = 'WATHEN'
  row_num = 5
  col_num = 5
  call wathen_order ( row_num, col_num, n )
  allocate ( a(1:n,1:n) )
  call wathen ( row_num, col_num, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  WILK03
!
  title = 'WILK03'
  n = 3
  allocate ( a(1:n,1:n) )
  call wilk03 ( a )
  call wilk03_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK04
!
  title = 'WILK04'
  n = 4
  allocate ( a(1:n,1:n) )
  call wilk04 ( a )
  call wilk04_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK05
!
  title = 'WILK05'
  n = 5
  allocate ( a(1:n,1:n) )
  call wilk05 ( a )
  call wilk05_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK12
!
  title = 'WILK12'
  n = 12
  allocate ( a(1:n,1:n) )
  call wilk12 ( a )
  call wilk12_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK20
!
  title = 'WILK20'
  n = 20
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call wilk20 ( alpha, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  WILK21
!
  title = 'WILK21'
  n = 21
  allocate ( a(1:n,1:n) )
  call wilk21 ( n, a )
  call wilk21_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILSON
!
  title = 'WILSON'
  n = 4
  allocate ( a(1:n,1:n) )
  call wilson ( a )
  call wilson_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ZERO
!
  title = 'ZERO'
  n = 5
  allocate ( a(1:n,1:n) )
  call zero ( n, n, a )
  call zero_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ZIELKE
!
  title = 'ZIELKE'
  n = 5
  allocate ( a(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  d1 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d2 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  d3 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call zielke ( n, d1, d2, d3, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )

  return
end
subroutine test_eigen_left ( )

!*****************************************************************************80
!
!! TEST_EIGEN_LEFT tests left eigensystems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable, dimension ( : ) :: d
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  real ( kind = 8 ), allocatable, dimension ( : ) :: lambda
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EIGEN_LEFT'
  write ( *, '(a)' ) '  Compute the Frobenius norm of the eigenvalue error:'
  write ( *, '(a)' ) '    X * A - LAMBDA * X'
  write ( *, '(a)' ) '  given K left eigenvectors X and eigenvalues LAMBDA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    N     K      ' // &
    '||A||          ||X*A-LAMBDA*X||'
  write ( *, '(a)' ) ''
!
!  A123
!
  title = 'A123'
  n = 3
  k = 3
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call a123 ( a )
  call a123_eigenvalues ( lambda )
  call a123_eigen_left ( x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CARRY
!
  title = 'CARRY'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  i4_lo = 2
  i4_hi = 20
  seed = 123456789
  i1 = i4_uniform_ab ( i4_lo, i4_hi, seed )
  call carry ( n, i1, a )
  call carry_eigenvalues ( n, i1, lambda )
  call carry_eigen_left ( n, i1, x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CHOW
!
  title = 'CHOW'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call chow ( alpha, beta, n, n, a )
  call chow_eigenvalues ( alpha, beta, n, lambda )
  call chow_eigen_left ( alpha, beta, n, x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  DIAGONAL
!
  title = 'DIAGONAL'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( d(1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, d )
  call diagonal ( n, n, d, a )
  call diagonal_eigenvalues ( n, d, lambda )
  call diagonal_eigen_left ( n, d, x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( d )
  deallocate ( lambda )
  deallocate ( x )
!
!  ROSSER1
!
  title = 'ROSSER1'
  n = 8
  k = 8
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rosser1 ( a )
  call rosser1_eigenvalues ( lambda )
  call rosser1_eigen_left ( x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  SYMM_RANDOM
!
  title = 'SYMM_RANDOM'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( d(1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, d )
  key = 123456789
  call symm_random ( n, d, key, a )
  call symm_random_eigenvalues ( n, d, key, lambda )
  call symm_random_eigen_left ( n, d, key, x )
  call r8mat_is_eigen_left ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( d )
  deallocate ( lambda )
  deallocate ( x )

  return
end
subroutine test_eigen_right ( )

!*****************************************************************************80
!
!! TEST_EIGEN_RIGHT tests right eigensystems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable, dimension ( : ) :: d
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  real ( kind = 8 ), allocatable, dimension ( : ) :: lambda
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EIGEN_RIGHT'
  write ( *, '(a)' ) '  Compute the Frobenius norm of the eigenvalue error:'
  write ( *, '(a)' ) '    A * X - X * LAMBDA'
  write ( *, '(a)' ) '  given K right eigenvectors X and eigenvalues LAMBDA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    N     K      ||A||          ||A*X-X*Lambda||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  n = 3
  k = 3
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  call a123 ( a )
  call a123_eigenvalues ( lambda )
  call a123_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  BAB
!
  title = 'BAB'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call bab ( n, alpha, beta, a )
  call bab_eigenvalues ( n, alpha, beta, lambda )
  call bab_eigen_right ( n, alpha, beta, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  BODEWIG
!
  title = 'BODEWIG'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  call bodewig ( a )
  call bodewig_eigenvalues ( lambda )
  call bodewig_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CARRY
!
  title = 'CARRY'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  i4_lo = 2
  i4_hi = 20
  seed = 123456789
  i1 = i4_uniform_ab ( i4_lo, i4_hi, seed )
  call carry ( n, i1, a )
  call carry_eigenvalues ( n, i1, lambda )
  call carry_eigen_right ( n, i1, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CHOW
!
  title = 'CHOW'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call chow ( alpha, beta, n, n, a )
  call chow_eigenvalues ( alpha, beta, n, lambda )
  call chow_eigen_right ( alpha, beta, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  COMBIN
!
  title = 'COMBIN'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call combin ( alpha, beta, n, a )
  call combin_eigenvalues ( alpha, beta, n, lambda )
  call combin_eigen_right ( alpha, beta, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  DIF2
!
  title = 'DIF2'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call dif2 ( n, n, a )
  call dif2_eigenvalues ( n, lambda )
  call dif2_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  EXCHANGE
!
  title = 'EXCHANGE'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call exchange ( n, n, a )
  call exchange_eigenvalues ( n, lambda )
  call exchange_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  IDEM_RANDOM
!
  title = 'IDEM_RANDOM'
  n = 5
  k = 5
  rank = 3
  key = 123456789
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call idem_random ( n, rank, key, a )
  call idem_random_eigenvalues ( n, rank, key, lambda )
  call idem_random_eigen_right ( n, rank, key, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  IDENTITY
!
  title = 'IDENTITY'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call identity ( n, n, a )
  call identity_eigenvalues ( n, lambda )
  call identity_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ILL3
!
  title = 'ILL3'
  n = 3
  k = 3
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call ill3 ( a )
  call ill3_eigenvalues ( lambda )
  call ill3_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  KERSHAW
!
  title = 'KERSHAW'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call kershaw ( a )
  call kershaw_eigenvalues ( lambda )
  call kershaw_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  KMS
!  Eigenvalue information requires 0 <= ALPHA <= 1.0.
!
  title = 'KMS'
  n = 5
  k = 5
  r8_lo = 0.0D+00
  r8_hi = +1.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call kms ( alpha, n, n, a )
  call kms_eigenvalues ( alpha, n, lambda )
  call kms_eigen_right ( alpha, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  LINE_ADJ
!
  title = 'LINE_ADJ'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  call line_adj ( n, a )
  call line_adj_eigenvalues ( n, lambda )
  call line_adj_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  LINE_LOOP_ADJ
!
  title = 'LINE_LOOP_ADJ'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  call line_loop_adj ( n, a )
  call line_loop_adj_eigenvalues ( n, lambda )
  call line_loop_adj_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ONE
!
  title = 'ONE'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call one ( n, n, a )
  call one_eigenvalues ( n, lambda )
  call one_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ORTEGA
!
  title = 'ORTEGA'
  n = 5
  k = n
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v1 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v2 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v3 )
  call ortega ( n, v1, v2, v3, a )
  call ortega_eigenvalues ( n, v1, v2, v3, lambda )
  call ortega_eigen_right ( n, v1, v2, v3, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
  deallocate ( x )
!
!  OTO
!
  title = 'OTO'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call oto ( n, n, a )
  call oto_eigenvalues ( n, lambda )
  call oto_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  PDS_RANDOM
!
  title = 'PDS_RANDOM'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  key = 123456789
  call pds_random ( n, key, a )
  call pds_random_eigenvalues ( n, key, lambda )
  call pds_random_eigen_right ( n, key, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  PEI
!
  title = 'PEI'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pei ( alpha, n, a )
  call pei_eigenvalues ( alpha, n, lambda )
  call pei_eigen_right ( alpha, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RODMAN
!
  title = 'RODMAN'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call rodman ( n, n, alpha, a )
  call rodman_eigenvalues ( n, alpha, lambda )
  call rodman_eigen_right ( n, alpha, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ROSSER1
!
  title = 'ROSSER1'
  n = 8
  k = 8
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rosser1 ( a )
  call rosser1_eigenvalues ( lambda )
  call rosser1_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS1
!
  title = 'RUTIS1'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis1 ( a )
  call rutis1_eigenvalues ( lambda )
  call rutis1_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS2
!
  title = 'RUTIS2'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis2 ( a )
  call rutis2_eigenvalues ( lambda )
  call rutis2_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS3
!  COMPLEX eigenvalues cannot be handled yet!
!
  if ( .false. ) then
    title = 'RUTIS3'
    n = 4
    k = 4
    allocate ( a(1:n,1:n) )
    allocate ( lambda(1:k) )
    allocate ( x(1:n,1:k) )
    call rutis3 ( a )
    call rutis3_eigenvalues ( lambda )
    call rutis3_eigen_right ( x )
!   call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
      title, n, k, norm_frobenius, error_frobenius
    deallocate ( a )
    deallocate ( lambda )
    deallocate ( x )

  end if
!
!  RUTIS5
!
  title = 'RUTIS5'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis5 ( a )
  call rutis5_eigenvalues ( lambda )
  call rutis5_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  SYLVESTER_KAC
!
  title = 'SYLVESTER_KAC'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call sylvester_kac ( n, a )
  call sylvester_kac_eigenvalues ( n, lambda )
  call sylvester_kac_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  SYMM_RANDOM
!
  title = 'SYMM_RANDOM'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  allocate ( d(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, d )
  key = 123456789
  call symm_random ( n, d, key, a )
  call symm_random_eigenvalues ( n, d, key, lambda )
  call symm_random_eigen_right ( n, d, key, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( d )
  deallocate ( lambda )
  deallocate ( x )
!
!  WILK12
!
  title = 'WILK12'
  n = 12
  k = 12
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call wilk12 ( a )
  call wilk12_eigenvalues ( lambda )
  call wilk12_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  WILSON
!
  title = 'WILSON'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call wilson ( a )
  call wilson_eigenvalues ( lambda )
  call wilson_eigen_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ZERO
!
  title = 'ZERO'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call zero ( n, n, a )
  call zero_eigenvalues ( n, lambda )
  call zero_eigen_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )

  return
end
subroutine test_inverse ( )

!*****************************************************************************80
!
!! TEST_INVERSE tests the inverse computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 8 ), allocatable, dimension ( : ) :: d
  real ( kind = 8 ) error_ab
  real ( kind = 8 ) error_ac
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) n
  real ( kind = 8 ) norma_frobenius
  real ( kind = 8 ) normc_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  integer ( kind = 4 ) x_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  integer ( kind = 4 ) y_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INVERSE'
  write ( *, '(a)' ) '  A = a test matrix of order N.'
  write ( *, '(a)' ) '  B = inverse as computed by a routine.'
  write ( *, '(a)' ) '  C = inverse as computed by R8MAT_INVERSE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A||    = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||C||    = Frobenius norm of C.'
  write ( *, '(a)' ) '  ||I-AC|| = Frobenius norm of I-A*C.'
  write ( *, '(a)' ) '  ||I-AB|| = Frobenius norm of I-A*B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    N    ' // &
                     '||A||      ||C||        ||I-AC||      ||I-AB||'
  write ( *, '(a)' ) ' '
!
!  AEGERTER
!
  title = 'AEGERTER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call aegerter ( n, a )
  call aegerter_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BAB
!
  title = 'BAB'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call bab ( n, alpha, beta, a )
  call bab_inverse ( n, alpha, beta, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BAUER
!
  title = 'BAUER'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call bauer ( a )
  call bauer_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BERNSTEIN
!
  title = 'BERNSTEIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call bernstein ( n, a )
  call bernstein_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BIS
!
  title = 'BIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call bis ( alpha, beta, n, n, a )
  call bis_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BIW
!
  title = 'BIW'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call biw ( n, a )
  call biw_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BODEWIG
!
  title = 'BODEWIG'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call bodewig ( a )
  call bodewig_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BOOTHROYD
!
  title = 'BOOTHROYD'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call boothroyd ( n, a )
  call boothroyd_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BORDERBAND
!
  title = 'BORDERBAND'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call borderband ( n, a )
  call borderband_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CARRY
!
  title = 'CARRY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  i4_lo = 2
  i4_hi = 20
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
  call carry ( n, k, a )
  call carry_inverse ( n, k, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CAUCHY
!
  title = 'CAUCHY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  call cauchy ( n, x, y, a )
  call cauchy_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  CHEBY_T
!
  title = 'CHEBY_T'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_t ( n, a )
  call cheby_t_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_U
!
  title = 'CHEBY_U'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_u ( n, a )
  call cheby_u_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_VAN2
!
  title = 'CHEBY_VAN2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_van2 ( n, a )
  call cheby_van2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_VAN3
!
  title = 'CHEBY_VAN3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_van3 ( n, a )
  call cheby_van3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHOW
!
  title = 'CHOW'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call chow ( alpha, beta, n, n, a )
  call chow_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CIRCULANT
!
  title = 'CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call circulant ( n, n, x, a )
  call circulant_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  CIRCULANT2
!
  title = 'CIRCULANT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CLEMENT1
!  N must be even.
!
  title = 'CLEMENT1'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CLEMENT2
!  N must be even.
!
  title = 'CLEMENT2'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, y )
  call clement2 ( n, x, y, a )
  call clement2_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  COMBIN
!
  title = 'COMBIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  COMPANION
!
  title = 'COMPANION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call companion ( n, x, a )
  call companion_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  COMPLEX_I
!
  title = 'COMPLEX_I'
  n = 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call complex_i ( a )
  call complex_i_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX1
!
  title = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX2
!
  title = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX3
!
  title = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX4
!
  title = 'CONEX4'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call conex4 ( a )
  call conex4_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONFERENCE
!  N-1 must be an odd prime or a power of an odd prime.
!
  title = 'CONFERENCE'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call conference ( n, a )
  call conference_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB2
!
  title = 'DAUB2'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub2 ( n, a )
  call daub2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB4
!
  title = 'DAUB4'
  n = 8
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub4 ( n, a )
  call daub4_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB6
!
  title = 'DAUB6'
  n = 12
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub6 ( n, a )
  call daub6_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB8
!
  title = 'DAUB8'
  n = 16
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub8 ( n, a )
  call daub8_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB10
!
  title = 'DAUB10'
  n = 20
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub10 ( n, a )
  call daub10_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB12
!
  title = 'DAUB12'
  n = 24
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub12 ( n, a )
  call daub12_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DIAGONAL
!
  title = 'DIAGONAL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call diagonal ( n, n, x, a )
  call diagonal_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  DIF1
!  N must be even.
!
  title = 'DIF1'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call dif1 ( n, n, a )
  call dif1_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DIF2
!
  title = 'DIF2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call dif2 ( n, n, a )
  call dif2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DORR
!
  title = 'DORR'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call dorr ( alpha, n, a )
  call dorr_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DOWNSHIFT
!
  title = 'DOWNSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call downshift ( n, a )
  call downshift_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  EULERIAN
!
  title = 'EULERIAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call eulerian ( n, n, a )
  call eulerian_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  EXCHANGE
!
  title = 'EXCHANGE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call exchange ( n, n, a )
  call exchange_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIBONACCI2
!
  title = 'FIBONACCI2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fibonacci2 ( n, a )
  call fibonacci2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIBONACCI3
!
  title = 'FIBONACCI3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fibonacci3 ( n, a )
  call fibonacci3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIEDLER.
!  The FIEDLER_INVERSE routine assumes the X vector is sorted.
!
  title = 'FIEDLER'
  n = 7
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call r8vec_sort_bubble_a ( n, x )
  call fiedler ( n, n, x, a )
  call fiedler_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  FORSYTHE
!
  title = 'FORSYTHE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call forsythe ( alpha, beta, n, a )
  call forsythe_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FOURIER_COSINE
!
  title = 'FOURIER_COSINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fourier_cosine ( n, a )
  call fourier_cosine_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FOURIER_SINE
!
  title = 'FOURIER_SINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fourier_sine ( n, a )
  call fourier_sine_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FRANK
!
  title = 'FRANK'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call frank ( n, a )
  call frank_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GFPP
!
  title = 'GFPP'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call gfpp ( n, alpha, a )
  call gfpp_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GIVENS
!
  title = 'GIVENS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call givens ( n, n, a )
  call givens_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK316
!
  title = 'GK316'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call gk316 ( n, a )
  call gk316_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK323
!
  title = 'GK323'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call gk323 ( n, n, a )
  call gk323_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK324
!
  title = 'GK324'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call gk324 ( n, n, x, a )
  call gk324_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  HANKEL_N
!
  title = 'HANKEL_N'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hankel_n ( n, a )
  call hankel_n_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HANOWA
!  N must be even.
!
  title = 'HANOWA'
  n = 6
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hanowa ( alpha, n, a )
  call hanowa_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HARMAN
!
  title = 'HARMAN'
  n = 8
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call harman ( a )
  call harman_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HARTLEY
!
  title = 'HARTLEY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hartley ( n, a )
  call hartley_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HELMERT
!
  title = 'HELMERT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call helmert ( n, a )
  call helmert_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HELMERT2
!
  title = 'HELMERT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call helmert2 ( n, x, a )
  call helmert2_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  HERMITE
!
  title = 'HERMITE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hermite ( n, a )
  call hermite_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HERNDON
!
  title = 'HERNDON'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call herndon ( n, a )
  call herndon_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HILBERT
!
  title = 'HILBERT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hilbert ( n, n, a )
  call hilbert_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HOUSEHOLDER
!
  title = 'HOUSEHOLDER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call householder ( n, x, a )
  call householder_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  IDENTITY
!
  title = 'IDENTITY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call identity ( n, n, a )
  call identity_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  ILL3
!
  title = 'ILL3'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call ill3 ( a )
  call ill3_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  INTEGRATION
!
  title = 'INTEGRATION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call integration ( alpha, n, a )
  call integration_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  INVOL
!
  title = 'INVOL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call invol ( n, a )
  call invol_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  JACOBI
!  N must be even.
!
  title = 'JACOBI'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call jacobi ( n, n, a )
  call jacobi_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  JORDAN
!
  title = 'JORDAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call jordan ( n, n, alpha, a )
  call jordan_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KAHAN
!
  title = 'KAHAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call kahan ( alpha, n, n, a )
  call kahan_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
   write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KERSHAW
!
  title = 'KERSHAW'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call kershaw ( a )
  call kershaw_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KERSHAWTRI
!
  title = 'KERSHAWTRI'
  n = 5
  x_n = ( n + 1 ) / 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:x_n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( x_n, r8_lo, r8_hi, seed, x )
  call kershawtri ( n, x, a )
  call kershawtri_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  KMS
!
  title = 'KMS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call kms ( alpha, n, n, a )
  call kms_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LAGUERRE
!
  title = 'LAGUERRE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call laguerre ( n, a )
  call laguerre_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LEGENDRE
!
  title = 'LEGENDRE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call legendre ( n, a )
  call legendre_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LEHMER
!
  title = 'LEHMER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lehmer ( n, n, a )
  call lehmer_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LESP
!
  title = 'LESP'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lesp ( n, n, a )
  call lesp_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LIETZKE
!
  title = 'LIETZKE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lietzke ( n, a )
  call lietzke_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LINE_ADJ
!  N must be even.
!
  title = 'LINE_ADJ'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call line_adj ( n, a )
  call line_adj_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LOTKIN
!
  title = 'LOTKIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lotkin ( n, n, a )
  call lotkin_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MAXIJ
!
  title = 'MAXIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call maxij ( n, n, a )
  call maxij_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MILNES
!
  title = 'MILNES'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call milnes ( n, n, x, a )
  call milnes_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  MINIJ
!
  title = 'MINIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call minij ( n, n, a )
  call minij_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MOLER1
!
  title = 'MOLER1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call moler1 ( alpha, n, n, a )
  call moler1_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MOLER3
!
  title = 'MOLER3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call moler3 ( n, n, a )
  call moler3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  ORTEGA
!
  title = 'ORTEGA'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v1 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v2 )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, v3 )
  call ortega ( n, v1, v2, v3, a )
  call ortega_inverse ( n, v1, v2, v3, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
!
!  ORTH_SYMM
!
  title = 'ORTH_SYMM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call orth_symm ( n, a )
  call orth_symm_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  OTO
!
  title = 'OTO'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call oto ( n, n, a )
  call oto_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PARTER
!
  title = 'PARTER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call parter ( n, n, a )
  call parter_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL1
!
  title = 'PASCAL1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call pascal1 ( n, a )
  call pascal1_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL2
!
  title = 'PASCAL2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call pascal2 ( n, a )
  call pascal2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL3
!
  title = 'PASCAL3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pascal3 ( n, alpha, a )
  call pascal3_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PDS_RANDOM
!
  title = 'PDS_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  key = 123456789
  call pds_random ( n, key, a )
  call pds_random_inverse ( n, key, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PEI
!
  title = 'PEI'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call pei ( alpha, n, a )
  call pei_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PERMUTATION_RANDOM
!
  title = 'PERMUTATION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  key = 123456789
  call permutation_random ( n, key, a )
  call permutation_random_inverse ( n, key, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PLU
!
  title = 'PLU'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( pivot(1:n) )
  seed = 123456789
  do i = 1, n
    i4_lo = i
    i4_hi = n
    pivot(i) = i4_uniform_ab ( i4_lo, i4_hi, seed )
  end do
  call plu ( n, pivot, a )
  call plu_inverse ( n, pivot, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( pivot )
!
!  RIS
!
  title = 'RIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call ris ( n, a )
  call ris_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RODMAN
!
  title = 'RODMAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call rodman ( n, n, alpha, a )
  call rodman_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS1
!
  title = 'RUTIS1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis1 ( a )
  call rutis1_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS2
!
  title = 'RUTIS2'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis2 ( a )
  call rutis2_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS3
!
  title = 'RUTIS3'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis3 ( a )
  call rutis3_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS4
!
  title = 'RUTIS4'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis4 ( n, a )
  call rutis4_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS5
!
  title = 'RUTIS5'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis5 ( a )
  call rutis5_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SCHUR_BLOCK
!
  title = 'SCHUR_BLOCK'
  n = 5
  x_n = ( n + 1 ) / 2
  y_n = n / 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:x_n) )
  allocate ( y(1:y_n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( x_n, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( y_n, r8_lo, r8_hi, seed, y )
  call schur_block ( n, x, y, a )
  call schur_block_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  SPLINE
!
  title = 'SPLINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call spline ( n, x, a )
  call spline_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  STIRLING
!
  title = 'STIRLING'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call stirling ( n, n, a )
  call stirling_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SUMMATION
!
  title = 'SUMMATION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call summation ( n, n, a )
  call summation_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SWEET1
!
  title = 'SWEET1'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call sweet1 ( a )
  call sweet1_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SWEET2
!
  title = 'SWEET2'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call sweet2 ( a )
  call sweet2_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SWEET3
!
  title = 'SWEET3'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call sweet3 ( a )
  call sweet3_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SWEET4
!
  title = 'SWEET4'
  n = 13
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call sweet4 ( a )
  call sweet4_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SYLVESTER_KAC
!  N must be even.
!
  title = 'SYLVESTER_KAC'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call sylvester_kac ( n, a )
  call sylvester_kac_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SYMM_RANDOM
!
  title = 'SYMM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( d(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, d )
  key = 123456789
  call symm_random ( n, d, key, a )
  call symm_random_inverse ( n, d, key, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( d )
!
!  TRI_UPPER
!
  title = 'TRI_UPPER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call tri_upper ( alpha, n, a )
  call tri_upper_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  TRIS
!
  title = 'TRIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  beta = r8_uniform_ab ( r8_lo, r8_hi, seed )
  gamma = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call tris ( n, n, alpha, beta, gamma, a )
  call tris_inverse ( n, alpha, beta, gamma, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  TRIV
!
  title = 'TRIV'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n) )
  allocate ( z(1:n-1) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, x )
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, y )
  call r8vec_uniform_ab ( n - 1, r8_lo, r8_hi, seed, z )
  call triv ( n, x, y, z, a )
  call triv_inverse ( n, x, y, z, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  TRIW
!
  title = 'TRIW'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  i4_lo = 0
  i4_hi = n - 1
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  k = i4_uniform_ab ( i4_lo, i4_hi, seed )
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call triw ( alpha, k, n, a )
  call triw_inverse ( alpha, k, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  UNITARY_RANDOM
!  Need to add C8_INVERSE(), or, more likely, push complex matrices out
!  to another package.
!
  title = 'UNITARY_RANDOM'
  n = 5
! allocate ( c8_a(1:n,1:n) )
! allocate ( c8_b(1:n,1:n) )
! allocate ( c8_c(1:n,1:n) )
! key = 123456789
! call unitary_random ( n, key, c8_a )
! call unitary_random_inverse ( n, key, c8_b )
! call c8mat_inverse ( n, c8_a, c8_c )
! call c8mat_is_inverse ( n, c8_a, c8_b, error_ab )
! call c8mat_is_inverse ( n, c8_a, c8_c, error_ac )
! norma_frobenius = c8mat_norm_fro ( n, n, a )
! normc_frobenius = c8mat_norm_fro ( n, n, c )
! write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
!   title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
! deallocate ( c8_a )
! deallocate ( c8_b )
! deallocate ( c8_c )
!
!  UPSHIFT
!
  title = 'UPSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call upshift ( n, a )
  call upshift_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  VAND1
!
  title = 'VAND1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call vand1 ( n, x, a )
  call vand1_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  VAND2
!
  title = 'VAND2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( n, r8_lo, r8_hi, seed, x )
  call vand2 ( n, x, a )
  call vand2_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  WILK03
!
  title = 'WILK03'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk03 ( a )
  call wilk03_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK04
!
  title = 'WILK04'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk04 ( a )
  call wilk04_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK05
!
  title = 'WILK05'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk05 ( a )
  call wilk05_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK21
!
  title = 'WILK21'
  n = 21
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk21 ( n, a )
  call wilk21_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILSON
!
  title = 'WILSON'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilson ( a )
  call wilson_inverse ( b )
  call 
r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g10.2,2x,g10.2,2x,g10.2,2x,g10.2)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine test_llt ( )

!*****************************************************************************80
!
!! TEST_LLT tests LLT factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ), allocatable :: l(:,:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST_LLT'
  write ( *, '(a)' ) '  A = a test matrix of order M by M'
  write ( *, '(a)' ) '  L is an M by N lower triangular Cholesky factor.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||A-LLT|| = Frobenius norm of A-L*L''.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Title                    M     N      ' // &
    '||A||            ||A-LLT||'
  write ( *, '(a)' ) ''
!
!  DIF2
!
  title = 'DIF2'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call dif2 ( m, n, a )
  call dif2_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  GIVENS
!
  title = 'GIVENS'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call givens ( m, n, a )
  call givens_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  KERSHAW
!
  title = 'KERSHAW'
  m = 4
  n = 4
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call kershaw ( a )
  call kershaw_llt ( l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  LEHMER
!
  title = 'LEHMER'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call lehmer ( n, n, a )
  call lehmer_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  MINIJ
!
  title = 'MINIJ'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call minij ( n, n, a )
  call minij_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  MOLER1
!
  title = 'MOLER1'
  m = 5
  n = 5
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call moler1 ( alpha, m, n, a )
  call moler1_llt ( alpha, n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  MOLER3
!
  title = 'MOLER3'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call moler3 ( m, n, a )
  call moler3_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  OTO
!
  title = 'OTO'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call oto ( m, n, a )
  call oto_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  PASCAL2
!
  title = 'PASCAL2'
  m = 5
  n = 5
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call pascal2 ( n, a )
  call pascal2_llt ( n, l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
!
!  WILSON
!
  title = 'WILSON'
  m = 4
  n = 4
  allocate ( a(1:m,1:m) )
  allocate ( l(1:m,1:n) )
  call wilson ( a )
  call wilson_llt ( l )
  call r8mat_is_llt ( m, n, a, l, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )

  return
end
subroutine test_null_left ( )

!*****************************************************************************80
!
!! TEST_NULL_LEFT tests left null vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) error_l2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ) norm_x_l2
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) r8vec_norm_l2
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NULL_LEFT'
  write ( *, '(a)' ) '  A = a test matrix of order M by N'
  write ( *, '(a)' ) '  x = an M vector, candidate for a left null vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||x|| = L2 norm of x.'
  write ( *, '(a)' ) '  ||x''*A||/||x|| = L2 norm of x''A over L2 norm of x.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    M     N      ' // &
               '||A||            ||x||        ||x''*A||/||x||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  m = 3
  n = 3
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  call a123 ( a )
  call a123_null_left ( x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CHEBY_DIFF1
!
  title = 'CHEBY_DIFF1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call cheby_diff1 ( n, a )
  call cheby_diff1_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CREATION
!
  title = 'CREATION'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call creation ( m, n, a )
  call creation_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1
!  Only has null vectors for M odd
!
  title = 'DIF1'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  call dif1 ( m, n, a )
  call dif1_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1CYCLIC
!
  title = 'DIF1CYCLIC'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif1cyclic ( n, a )
  call dif1cyclic_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF2CYCLIC
!
  title = 'DIF2CYCLIC'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  call dif2cyclic ( n, a )
  call dif2cyclic_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  EBERLEIN
!
  title = 'EBERLEIN'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call eberlein ( alpha, n, a )
  call eberlein_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  FIBONACCI1
!
  title = 'FIBONACCI1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  f1 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  f2 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call fibonacci1 ( n, f1, f2, a )
  call fibonacci1_null_left ( m, n, f1, f2, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  LAUCHLI
!
  title = 'LAUCHLI'
  m = 6
  n = m - 1
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call lauchli ( alpha, m, n, a )
  call lauchli_null_left ( alpha, m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  LINE_ADJ
!
  title = 'LINE_ADJ'
  m = 7
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call line_adj ( n, a )
  call line_adj_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  MOLER2
!
  title = 'MOLER2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call moler2 ( a )
  call moler2_null_left ( x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ONE
!
  title = 'ONE'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  call one ( m, n, a )
  call one_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  RING_ADJ
!  M must be a multiple of 4 for there to be a null vector.
!
  title = 'RING_ADJ'
  m = 12
  n = 12
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call ring_adj ( n, a )
  call ring_adj_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ROSSER1
!
  title = 'ROSSER1'
  m = 8
  n = 8
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call rosser1 ( a )
  call rosser1_null_left ( x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ZERO
!
  title = 'ZERO'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:m) )
  call zero ( m, n, a )
  call zero_null_left ( m, n, x )
  call r8mat_is_null_left ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( m, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )

  return
end
subroutine test_null_right ( )

!*****************************************************************************80
!
!! TEST_NULL_RIGHT tests right null vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) col_num
  real ( kind = 8 ) error_l2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ) norm_x_l2
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) r8vec_norm_l2
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NULL_RIGHT'
  write ( *, '(a)' ) '  A = a test matrix of order M by N'
  write ( *, '(a)' ) '  x = an N vector, candidate for a right null vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||x|| = L2 norm of x.'
  write ( *, '(a)' ) '  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    M     N      ' // &
               '||A||            ||x||        ||A*x||/||x||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  m = 3
  n = 3
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call a123 ( a )
  call a123_null_right ( x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ARCHIMEDES
!
  title = 'ARCHIMEDES'
  m = 7
  n = 8
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call archimedes ( a )
  call archimedes_null_right ( x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CHEBY_DIFF1
!
  title = 'CHEBY_DIFF1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call cheby_diff1 ( n, a )
  call cheby_diff1_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CREATION
!
  title = 'CREATION'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call creation ( m, n, a )
  call creation_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1
!  Only has null vectors for N odd.
!
  title = 'DIF1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif1 ( m, n, a )
  call dif1_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1CYCLIC
!
  title = 'DIF1CYCLIC'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif1cyclic ( n, a )
  call dif1cyclic_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF2CYCLIC
!
  title = 'DIF2CYCLIC'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif2cyclic ( n, a )
  call dif2cyclic_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  FIBONACCI1
!
  title = 'FIBONACCI1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  f1 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  f2 = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call fibonacci1 ( n, f1, f2, a )
  call fibonacci1_null_right ( m, n, f1, f2, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  HAMMING
!
  title = 'HAMMING'
  m = 5
  n = 2 ** m - 1
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call hamming ( m, n, a )
  call hamming_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  LINE_ADJ
!
  title = 'LINE_ADJ'
  m = 7
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call line_adj ( n, a )
  call line_adj_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  MOLER2
!
  title = 'MOLER2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call moler2 ( a )
  call moler2_null_right ( x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  NEUMANN
!
  title = 'NEUMANN'
  row_num = 5
  col_num = 5
  m = row_num * col_num
  n = row_num * col_num
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call neumann ( row_num, col_num, a )
  call neumann_null_right ( row_num, col_num, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ONE
!
  title = 'ONE'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call one ( m, n, a )
  call one_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  RING_ADJ
!  N must be a multiple of 4 for there to be a null vector.
!
  title = 'RING_ADJ'
  m = 12
  n = 12
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call ring_adj ( n, a )
  call ring_adj_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ROSSER1
!
  title = 'ROSSER1'
  m = 8
  n = 8
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call rosser1 ( a )
  call rosser1_null_right ( x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ZERO
!
  title = 'ZERO'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call zero ( m, n, a )
  call zero_null_right ( m, n, x )
  call r8mat_is_null_right ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )

  return
end
subroutine test_plu ( )

!*****************************************************************************80
!
!! TEST_PLU tests the PLU factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PLU'
  write ( *, '(a)' ) '  A = a test matrix of order M by N'
  write ( *, '(a)' ) '  P, L, U are the PLU factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||A-PLU|| = Frobenius norm of A-P*L*U.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    M     N      ' // &
               '||A||            ||A-PLU||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  m = 3
  n = 3
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call a123 ( a )
  call a123_plu ( p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  BODEWIG
!
  title = 'BODEWIG'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call bodewig ( a )
  call bodewig_plu ( p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  BORDERBAND
!
  title = 'BORDERBAND'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call borderband ( n, a )
  call borderband_plu ( n, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  DIF2
!
  title = 'DIF2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call dif2 ( m, n, a )
  call dif2_plu ( n, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  GFPP
!
  title = 'GFPP'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call gfpp ( n, alpha, a )
  call gfpp_plu ( n, alpha, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  GIVENS
!
  title = 'GIVENS'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call givens ( n, n, a )
  call givens_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  KMS
!
  title = 'KMS'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call kms ( alpha, m, n, a )
  call kms_plu ( alpha, n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  LEHMER
!
  title = 'LEHMER'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call lehmer ( n, n, a )
  call lehmer_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MAXIJ
!
  title = 'MAXIJ'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call maxij ( n, n, a )
  call maxij_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MINIJ
!
  title = 'MINIJ'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call minij ( n, n, a )
  call minij_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MOLER1
!
  title = 'MOLER1'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  alpha = r8_uniform_ab ( r8_lo, r8_hi, seed )
  call moler1 ( alpha, n, n, a )
  call moler1_plu ( alpha, n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MOLER3
!
  title = 'MOLER3'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call moler3 ( m, n, a )
  call moler3_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  OTO
!
  title = 'OTO'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call oto ( m, n, a )
  call oto_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  PLU
!
  title = 'PLU'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( l(1:n,1:n) )
  allocate ( p(1:n,1:n) )
  allocate ( pivot(n) )
  allocate ( u(1:n,1:n) )
  seed = 123456789
  do i = 1, n
    i4_lo = i
    i4_hi = n
    pivot(i) = i4_uniform_ab ( i4_lo, i4_hi, seed )
  end do
  call plu ( n, pivot, a )
  call plu_plu ( n, pivot, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  PASCAL2
!
  title = 'PASCAL2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call pascal2 ( n, a )
  call pascal2_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  VAND2
!
  title = 'VAND2'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  allocate ( x(1:m) )
  r8_lo = -5.0D+00
  r8_hi = +5.0D+00
  seed = 123456789
  call r8vec_uniform_ab ( m, r8_lo, r8_hi, seed, x )
  call vand2 ( m, x, a )
  call vand2_plu ( m, x, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
  deallocate ( x )
!
!  WILSON
!
  title = 'WILSON'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call wilson ( a )
  call wilson_plu ( p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )

  return
end
subroutine test_solution ( )

!*****************************************************************************80
!
!! TEST_SOLUTION tests the linear solution computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ) beta
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncol
  real ( kind = 8 ) norm_frobenius
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) r8_hi
  real ( kind = 8 ) r8_lo
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SOLUTION'
  write ( *, '(a)' ) '  Compute the Frobenius norm of the solution error:'
  write ( *, '(a)' ) '    A * X - B'
  write ( *, '(a)' ) '  given MxN matrix A, NxK solution X, MxK right hand side B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Title                    M     N     K      ||A||         ||A*X-B||'
  write ( *, '(a)' ) ' '
!
!  A123
!
  title = 'A123'
  m = 3
  n = 3
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call a123 ( a )
  call a123_rhs ( b )
  call a123_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  BODEWIG
!
  title = 'BODEWIG'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call bodewig ( a )
  call bodewig_rhs ( b )
  call bodewig_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  DIF2
!
  title = 'DIF2'
  m = 10
  n = 10
  k = 2
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call dif2 ( m, n, a )
  call dif2_rhs ( m, k, b )
  call dif2_solution ( n, k, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  FRANK
!
  title = 'FRANK'
  m = 10
  n = 10
  k = 2
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call frank ( n, a )
  call frank_rhs ( m, k, b )
  call frank_solution ( n, k, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  POISSON
!
  title = 'POISSON'
  nrow = 4
  ncol = 5
  m = nrow * ncol
  n = nrow * ncol
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call poisson ( nrow, ncol, a )
  call poisson_rhs ( nrow, ncol, b )
  call poisson_solution ( nrow, ncol, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILK03
!
  title = 'WILK03'
  m = 3
  n = 3
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilk03 ( a )
  call wilk03_rhs ( b )
  call wilk03_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILK04
!
  title = 'WILK04'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilk04 ( a )
  call wilk04_rhs ( b )
  call wilk04_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILSON
!
  title = 'WILSON'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilson ( a )
  call wilson_rhs ( b )
  call wilson_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )

  return
end
subroutine test_type ( )

!*****************************************************************************80
!
!! TEST_TYPE tests functions which test the type of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST_TYPE'
  write ( *, '(a)' ) '  Test functions that query the type of a matrix.'
!
!  R8MAT_IS_TRANSITION.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Title                    M     N     ||A||' // &
    '            ||Transition Error||'
  write ( *, '(a)' ) ''

  title = 'BODEWIG'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  call bodewig ( a )
  call r8mat_is_transition ( m, n, a, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_frobenius, error_frobenius
  deallocate ( a )

  title = 'SNAKES'
  m = 101
  n = 101
  allocate ( a(1:m,1:n) )
  call snakes ( a )
  call r8mat_is_transition ( m, n, a, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_frobenius, error_frobenius
  deallocate ( a )

  title = 'TRANSITION_RANDOM'
  m = 5
  n = 5
  key = 123456789
  allocate ( a(1:m,1:n) )
  call transition_random ( n, key, a )
  call r8mat_is_transition ( m, n, a, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_frobenius, error_frobenius
  deallocate ( a )

  return
end
