program main

!*****************************************************************************80
!
!! MAIN is the main program for GEQP3_PRB.
!
!  Discussion:
!
!    GEQP3_PRB tests the GEQP3 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEQP3_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the GEQP3 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEQP3_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DGEQP3 for QR factorization of a real ( kind = 8 ) real matrix.
!
!  Discussion:
!
!    DGEQP3 requires that N <= M.
!
!    This example uses a 6x5 linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mmax = 8
  integer ( kind = 4 ), parameter :: nb = 64
  integer ( kind = 4 ), parameter :: nmax = 8
  integer ( kind = 4 ), parameter :: nrhsmx = 2

  integer ( kind = 4 ), parameter :: lda = mmax
  integer ( kind = 4 ), parameter :: ldb = mmax
  integer ( kind = 4 ), parameter :: lwork = 2 * nmax + ( nmax + 1 ) * nb

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nrhs = 2

  real ( kind = 8 ) a(lda,nmax)
  real ( kind = 8 ), save :: a_data(m,n) = reshape ( (/ &
    -0.09D+00, -1.56D+00, -1.48D+00, -1.09D+00,  0.08D+00, -1.59D+00, &
     0.14D+00,  0.20D+00, -0.43D+00,  0.84D+00,  0.55D+00, -0.72D+00, &
    -0.46D+00,  0.29D+00,  0.89D+00,  0.77D+00, -1.13D+00,  1.06D+00, &
     0.68D+00,  1.09D+00, -0.71D+00,  2.11D+00,  0.14D+00,  1.24D+00, &
     1.29D+00,  0.51D+00, -0.96D+00, -1.27D+00,  1.74D+00,  0.34D+00  &
    /), (/ m, n /) )
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(ldb,nrhsmx)
  real ( kind = 8 ), save :: b_data(m,nrhs) = reshape ( (/ &
    7.4D+00,  4.2D+00, -8.3D+00, 1.8D+00, 8.6D+00,  2.1D+00, &
    2.7D+00, -3.0D+00, -9.6D+00, 1.1D+00, 4.0D+00, -5.7D+00 /), &
    (/ m, nrhs /) )
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(nmax)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) rank
  real ( kind = 8 ) r(m,nrhs)
  real ( kind = 8 ) rnorm(nmax)
  real ( kind = 8 ) tau(nmax)
  real ( kind = 8 ), parameter :: tol = 0.01D+00
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) x(n,nrhs)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DGEQP3 computes the QR factorization,'
  write ( *, '(a)' ) '  with column pivoting,'
  write ( *, '(a)' ) '  of a M by N rectangular matrix A = Q * R,'
  write ( *, '(a)' ) '  with N <= M,'
  write ( *, '(a)' ) '  using real ( kind = 8 ) arithmetic.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Least squares solutions of A*X=B can then'
  write ( *, '(a)' ) '  be computed by calling:'
  write ( *, '(a)' ) '  DORMQR to compute QTB = Q'' * B, and'
  write ( *, '(a)' ) '  DTRSM to solve R * X = QTB.'

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  M = ', m
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a,i4)' ) '  NRHS = ', nrhs
  write ( *, '(a,g14.6)' ) '  TOL = ', tol
!
!  Transfer the packed arrays A_DATA and B_DATA into the loose arrays A and B.
!
  call r8mat_to_r8cmat ( lda, m, n,    a_data, a )
  call r8cmat_print ( lda, m, n, a, '  The R8CMAT A:' )

  call r8mat_to_r8cmat ( ldb, m, nrhs, b_data, b )
  call r8cmat_print ( ldb, m, nrhs, b, '  The R8CMAT B:' )
!
!  Set JPVT zero, which means that all columns are free.
!
  jpvt(1:n) = 0
!
!  Compute the QR factorization of A.
!
  call dgeqp3 ( m, n, a, lda, jpvt, tau, work, lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'GEQP3_PRB - Fatal error!'
    write ( *, '(a,i4)' ) '  DGEQP3 returned INFO = ', info
    stop 1
  end if
!
!  Determine and print the rank of the R factor,
!  using the relative tolerance TOL.
!
  rank = min ( m, n )

  do k = 1, min ( m, n )

    if ( abs ( a(k,k) ) <= tol * abs ( a(1,1) ) ) then
      rank = k - 1
      exit
    end if

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Estimated rank of A = ', rank
!
!  Compute C = (C1) = Q' * B, storing the result in B
!              (C2)
!
!  Here:
!    Q is of order MxM,
!    B is of order MxNRHS,
!    C will be of order MxNRHS.
!
  call dormqr ( 'Left', 'Transpose', m, nrhs, n, a, lda, tau, b, &
    ldb, work, lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'GEQP3_PRB - Fatal error!'
    write ( *, '(a,i4)' ) '  DORMQR returned INFO = ', info
    stop 1
  end if
!
!  Compute least-squares solutions by solving the system
!    R(1:K,1:K) * B = C1.
!  where
!    R is the leading KxK upper triangular submatrix in A,
!    C1 is a K by NRHS vector,
!    B is the K by NRHS solution.
!
  alpha = 1.0D+00
  k = rank

  call dtrsm ( 'Left', 'Upper', 'No transpose', 'Non-Unit', k, &
     nrhs, alpha, a, lda, b, ldb )
!
!  Estimate the square roots of the residual sums of squares 
!  (the 2-norm of each of the columns of C2)
!
  do j = 1, nrhs
    rnorm(j) = dnrm2 ( m - k, b(k+1:m,j), 1 )
  end do
!
!  Set the remaining elements of the solutions to zero (to give
!  the basic solutions.)
!
  b(k+1:n,1:nrhs) = 0.0D+00
!
!  Permute the least-squares solutions stored in B to give X = P*Y
!
  do j = 1, nrhs
    work(jpvt(1:n)) = b(1:n,j)
    b(1:n,j) = work(1:n)
  end do
!
!  Print the NRHS least-squares solution vectors of length N.
!
  call r8cmat_print ( ldb, n, nrhs, b, '  Least-squares solutions:' )
!
!  Print the square roots of the residual sums of squares
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Square root of residual sums of squares:'
  write ( *, '(a)' ) ' '
  write ( *, '(6x,6e14.6)' ) rnorm(1:nrhs)
!
!  Compute the NRHS residual vectors explicitly from the original data.
!
  r(1:m,1:nrhs) = matmul ( a_data(1:m,1:n), b(1:n,1:nrhs) )
  r(1:m,1:nrhs) = r(1:m,1:nrhs) - b_data(1:m,1:nrhs)
  call r8cmat_print ( m, m, nrhs, r, '  Residuals:' )

  return
end
subroutine r8cmat_print ( lda, m, n, a, title )

!*****************************************************************************80
!
!! R8CMAT_PRINT prints an R8CMAT.
!
!  Discussion:
!
!    An R8CMAT is an MxN array of R8's, stored with a leading dimension
!    of LD, and hence accessed either as a double indexed array:
!      (I,J) -> (I,J) 
!    or as a vector:
!      (I,J) -> (I+J*LD).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(LDA,N), the M by N matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  character ( len = * ) title

  call r8cmat_print_some ( lda, m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8cmat_print_some ( lda, m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CMAT_PRINT_SOME prints some of an R8CMAT.
!
!  Discussion:
!
!    An R8CMAT is an MxN array of R8's, stored with a leading dimension
!    of LD, and hence accessed either as a double indexed array:
!      (I,J) -> (I,J) 
!    or as a vector:
!      (I,J) -> (I+J*LD).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(LDA,N), the M by N matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_to_r8cmat ( lda, m, n, a1, a2 )

!*****************************************************************************80
!
!! R8MAT_TO_R8CMAT transfers data from an R8MAT to an R8CMAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, 
!    accessible as a vector:
!      (I,J) -> (I+J*M).
!    or as a doubly-dimensioned array, if declared A(M,N):
!      (I,J) -> A(I,J)      
!
!    An R8CMAT is an MxN array of R8's, stored with a leading dimension LD,
!    accessible as a vector:
!      (I,J) -> (I+J*LD).
!    or as a doubly-dimensioned array, if declared A(LD,N):
!      (I,J) -> A(I,J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A2.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!    M <= LDA.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) A1(M,N), the matrix to be copied.
!
!    Output, real ( kind = 8 ) A2(LDA,N), contains a copy of the
!    information in A1, in the MxN submatrix.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(m,n)
  real ( kind = 8 ) a2(lda,n)

  a2(1:m,1:n) = a1(1:m,1:n)
 
  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
