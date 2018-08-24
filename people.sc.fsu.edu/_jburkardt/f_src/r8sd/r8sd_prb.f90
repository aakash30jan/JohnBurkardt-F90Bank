program main

!*****************************************************************************80
!
!! MAIN is the main program for R8SD_PRB.
!
!  Discussion:
!
!    R8SD_PRB tests the R8SD library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8SD library.'

  call r8sd_cg_test ( )
  call r8sd_cg_test2 ( )
  call r8sd_dif2_test ( )
  call r8sd_indicator_test ( )
  call r8sd_mv_test ( )
  call r8sd_print_test ( )
  call r8sd_print_some_test ( )
  call r8sd_random_test ( )
  call r8sd_res_test ( )
  call r8sd_to_r8ge_test ( )
  call r8sd_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8sd_cg_test ( )

!*****************************************************************************80
!
!! R8SD_CG_TEST tests R8SD_CG.
!
!  Discussion:
!
!    NX and NY are the number of grid points in the X and Y directions.
!    N is the number of unknowns.
!    NDIAG is the number of nonzero diagonals we will store.  We only
!    store the main diagonal, and the superdiagonals.
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

  integer ( kind = 4 ), parameter :: ndiag = 3
  integer ( kind = 4 ), parameter :: nx = 10
  integer ( kind = 4 ), parameter :: ny = 10
  integer ( kind = 4 ), parameter :: n = nx * ny

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, nx /)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_CG_TEST'
  write ( *, '(a)' ) '  R8SD_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '  to a symmetric positive definite linear'
  write ( *, '(a)' ) '  system stored by diagonals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals is ', ndiag
  write ( *, '(a)' ) ' '
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0D+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0D+00
      else
        a(k,2) = -1.0D+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0D+00
      else
        a(k,3) = -1.0D+00
      end if

    end do
  end do
!
!  Print some of the matrix.
!
  call r8sd_print_some ( n, ndiag, offset, a, 1, 1, 10, 10, &
    '  First 10 rows and columns of matrix.' )
!
!  Set the desired solution.
!
  k = 0
  do j = 1, ny
    do i = 1, nx
      k = k + 1
      x(k) = real ( 10 * i + j, kind = 8 )
    end do
  end do
!
!  Compute the corresponding right hand side.
!
  call r8sd_mv ( n, n, ndiag, offset, a, x, b )

  call r8vec_print_some ( n, b, 10, '  Right hand side:' )
!
!  Set X to zero so no one accuses us of cheating.
!
  x(1:n) = 0.0D+00
!
!  Call the conjugate gradient method.
!
  call r8sd_cg ( n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call r8sd_mv ( n, n, ndiag, offset, a, x, b2 )
 
  err = maxval ( abs ( b2(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, 10, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual = ', err
!
!  Note that if we're not satisfied with the solution, we can
!  call again, using the computed X as our starting estimate.
!
!
!  Call the conjugate gradient method AGAIN.
!
  call r8sd_cg ( n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call r8sd_mv ( n, n, ndiag, offset, a, x, b2 )
 
  err = maxval ( abs ( b2(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, 10, '  Second attempt at solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual of second attempt = ', err

  return
end
subroutine r8sd_cg_test2 ( )

!*****************************************************************************80
!
!! R8SD_CG_TEST2 tests R8SD_CG.
!
!  Discussion:
!
!    This is a sample demonstration of how to compute some eigenvalues
!    and corresponding eigenvectors of a matrix.  The matrix is the
!    discretized Laplacian operator, which can be stored by diagonals,
!    and handled by the conjugate gradient method.
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

  integer ( kind = 4 ), parameter :: maxvec = 3
  integer ( kind = 4 ), parameter :: ndiag = 3
  integer ( kind = 4 ), parameter :: nx = 10
  integer ( kind = 4 ), parameter :: ny = 10
  integer ( kind = 4 ), parameter :: n = nx * ny
  real ( kind = 8 ), parameter :: pi = 3.141592653589D+00

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) del
  real ( kind = 8 ) dot
  real ( kind = 8 ) eval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ivec
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  real ( kind = 8 ) lambda_old
  real ( kind = 8 ) lamvec(maxvec)
  real ( kind = 8 ) norm
  integer ( kind = 4 ) nvec
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) vec(n,maxvec)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xnew(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_CG_TEST2'
  write ( *, '(a)' ) '  R8SD_CG is used for linear equation solving'
  write ( *, '(a)' ) '  in a demonstration of inverse iteration to'
  write ( *, '(a)' ) '  compute eigenvalues and eigenvectors of a'
  write ( *, '(a)' ) '  symmetric matrix stored by diagonals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are 25 of the smallest eigenvalues:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, J, eigenvalue(I,J):'
  write ( *, '(a)' ) ''

  do i = 1, min ( 5, nx )
    do j = 1, min ( 5, ny )
      eval = 4.0D+00 - 2.0D+00 * cos ( real ( i, kind = 8 ) * pi / real ( nx + 1, kind = 8 ) ) &
                     - 2.0D+00 * cos ( real ( j, kind = 8 ) * pi / real ( ny + 1, kind = 8 ) )
      write ( *, '(2i8,g14.6)' ) i, j, eval
    end do
  end do
!
!  OFFSET tells us where the nonzero diagonals are.  It does this
!  by recording how "high" or to the right the diagonals are from
!  the main diagonal.
!
  offset(1) =   0
  offset(2) =   1
  offset(3) =  nx
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0D+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0D+00
      else
        a(k,2) = -1.0D+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0D+00
      else
        a(k,3) = -1.0D+00
      end if

    end do
  end do

  nvec = 0
!
!  Set the starting eigenvector and eigenvalue estimates.
!
  do 

    lambda = 0.0D+00

    k = 0
    do j = 1, ny
      do i = 1, nx
        k = k + 1
        x(k) = 1.0D+00
      end do
    end do
!
!  Remove any components of previous eigenvectors.
!
    do ivec = 1, nvec
      dot = dot_product ( x(1:n), vec(1:n,ivec) )
      x(1:n) = x(1:n) - dot * vec(1:n,ivec)
    end do

    xnew(1:n) = x(1:n)
!
!  Iterate
!
    do iter = 1, 40

      norm = sqrt ( sum ( xnew(1:n)**2 ) )

      xnew(1:n) = xnew(1:n) / norm

      lambda_old = lambda
      lambda = 1.0D+00 / norm
!
!  Check for convergence.
!
      if ( 1 < iter ) then
        del = abs ( lambda - lambda_old )
        if ( del < 0.000001D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,g14.6)' ) 'Lambda estimate = ', lambda
          write ( *, '(a,i8)' ) 'Converged on step ', iter
          exit
        end if
      end if
!
!  Call the conjugate gradient method, solving
!    A * XNEW = X.
!
      x(1:n) = xnew(1:n)

      call r8sd_cg ( n, ndiag, offset, a, x, xnew  )

      do ivec = 1, nvec
        dot = dot_product ( xnew(1:n), vec(1:n,ivec) )
        xnew(1:n) = xnew(1:n) - dot * vec(1:n,ivec)
      end do

    end do

    nvec = nvec + 1
    lamvec(nvec) = lambda
    vec(1:n,nvec) = xnew(1:n)

    if ( maxvec <= nvec ) then
      exit
    end if

  end do

  return
end
subroutine r8sd_dif2_test ( )

!*****************************************************************************80
!
!! R8SD_DIF2_TEST tests R8SD_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ndiag = 2

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_DIF2_TEST'
  write ( *, '(a)' ) '  R8SD_DIF2 sets up an R8SD second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_dif2 ( n, n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD second difference matrix:' )

  return
end
subroutine r8sd_indicator_test ( )

!*****************************************************************************80
!
!! R8SD_INDICATOR_TEST tests R8SD_INDICATOR.
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
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8SD_INDICATOR sets up an R8SD indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_indicator ( n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD indicator matrix:' )

  return
end
subroutine r8sd_mv_test ( )

!*****************************************************************************80
!
!! R8SD_MV_TEST tests R8SD_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(N)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_MV_TEST'
  write ( *, '(a)' ) '  R8SD_MV computes b=A*x, where A is an R8SD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_indicator ( n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sd_mv ( n, n, ndiag, offset, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x' )

  return
end
subroutine r8sd_print_test ( )

!*****************************************************************************80
!
!! R8SD_PRINT_TEST tests R8SD_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_PRINT_TEST'
  write ( *, '(a)' ) '  R8SD_PRINT prints an R8SD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_indicator ( n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD matrix:' )

  return
end
subroutine r8sd_print_some_test ( )

!*****************************************************************************80
!
!! R8SD_PRINT_SOME_TEST tests R8SD_PRINT_SOME.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8SD_PRINT_SOME prints some of an R8SD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_indicator ( n, ndiag, offset, a )

  call r8sd_print_some ( n, ndiag, offset, a, 2, 3, 8, 7, '  Rows 2-8, Cols 3-7:' )

  return
end
subroutine r8sd_random_test ( )

!*****************************************************************************80
!
!! R8SD_RANDOM_TEST tests R8SD_RANDOM.
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

  integer ( kind = 4 ), parameter :: n = 15
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_RANDOM_TEST'
  write ( *, '(a)' ) '  R8SD_RANDOM randomizes an R8SD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  seed = 123456789
  call r8sd_random ( n, ndiag, offset, seed, a )

  call r8sd_print ( n, ndiag, offset, a, '  The random R8SD matrix:' )

  return
end
subroutine r8sd_res_test ( )

!*****************************************************************************80
!
!! R8SD_RES_TEST tests R8SD_RES.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 2

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1 /)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_RES_TEST'
  write ( *, '(a)' ) '  R8SD_RES computes a residual R=b-A*x'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals is ', ndiag
  write ( *, '(a)' ) ' '

  seed = 123456789
  call r8sd_random ( n, ndiag, offset, seed, a )
!
!  Print some of the matrix.
!
  call r8sd_print ( n, ndiag, offset, a, '  The R8SD matrix:' )
  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sd_mv ( n, n, ndiag, offset, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x' )
!
!  Make X2, a bad copy of X.
!
  do i = 1, n
    x2(i) = x(i) + 0.1D+00 * r8_uniform_01 ( seed )
  end do
  call r8vec_print ( n, x2, '  The defective vector X2:' )
!
!  Compute R = B-A*X2.
!
  call r8sd_res ( n, n, ndiag, offset, a, x2, b, r )
  
  call r8vec_print ( n, r, '  Residual b-A*x2:' )

  return
end
subroutine r8sd_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8SD_TO_R8GE_TEST tests R8SD_TO_R8GE.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8SD_TO_R8GER converts an R8SD matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_indicator ( n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD matrix:' )

  call r8sd_to_r8ge ( n, ndiag, offset, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8sd_zeros_test ( )

!*****************************************************************************80
!
!! R8SD_ZEROS_TEST tests R8SD_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SD_ZEROS_TEST'
  write ( *, '(a)' ) '  R8SD_ZEROS zeros an R8SD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call r8sd_zeros ( n, ndiag, offset, a )

  call r8sd_print ( n, ndiag, offset, a, '  The R8SD zero matrix:' )

  return
end
