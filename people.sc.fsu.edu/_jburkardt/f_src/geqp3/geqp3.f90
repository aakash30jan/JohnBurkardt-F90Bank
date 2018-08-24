subroutine dcopy ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DCOPY copies a vector X to a vector Y.
!
!  Discussion:
!
!    This routine uses real ( kind = 8 ) real arithmetic.
!
!    The routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of DX.
!
!    Output, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of DY.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)

  if ( n <= 0 ) then
    return
  end if

  if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 7 )

    if ( m /= 0 ) then

      dy(1:m) = dx(1:m)

    end if

    do i = m + 1, n, 7
      dy(i) = dx(i)
      dy(i+1) = dx(i+1)
      dy(i+2) = dx(i+2)
      dy(i+3) = dx(i+3)
      dy(i+4) = dx(i+4)
      dy(i+5) = dx(i+5)
      dy(i+6) = dx(i+6)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( -n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( -n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, &
  ldc )

!*****************************************************************************80
!
!! DGEMM computes C = alpha * A * B and related operations.
!
!  Discussion:
!
!    DGEMM performs one of the matrix-matrix operations
!
!     C := alpha * op ( A ) * op ( B ) + beta * C,
!
!    where op ( X ) is one of
!
!      op ( X ) = X   or   op ( X ) = X',
!
!    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
!    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!    
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character TRANSA, specifies the form of op( A ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n', op ( A ) = A.
!    'T' or 't', op ( A ) = A'.
!    'C' or 'c', op ( A ) = A'.
!
!    Input, character TRANSB, specifies the form of op ( B ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n', op ( B ) = B.
!    'T' or 't', op ( B ) = B'.
!    'C' or 'c', op ( B ) = B'.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the  matrix op ( A ) 
!    and of the matrix C.  0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number  of columns of the matrix 
!    op ( B ) and the number of columns of the matrix C.  0 <= N.
!
!    Input, integer ( kind = 4 ) K, the number of columns of the matrix op ( A )
!    and the number of rows of the matrix op ( B ).  0 <= K.
!
!    Input, real ( kind = 8 ) ALPHA, the scalar multiplier 
!    for op ( A ) * op ( B ).
!
!    Input, real ( kind = 8 ) A(LDA,KA), where:
!    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
!    part of the array contains A;
!    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
!    K by M part of the array must contain the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of A as declared in 
!    the calling routine.  When TRANSA = 'N' or 'n' then LDA must be at least 
!    max ( 1, M ), otherwise LDA must be at least max ( 1, K ).
!
!    Input, real ( kind = 8 ) B(LDB,KB), where:
!    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
!    part of the array contains B;
!    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
!    n by k  part of the array must contain the matrix B.
!
!    Input, integer ( kind = 4 ) LDB, the first dimension of B as declared in
!    the calling routine.  When TRANSB = 'N' or 'n' then LDB must be at least 
!    max ( 1, K ), otherwise LDB must be at least max ( 1, N ).
!
!    Input, real ( kind = 8 ) BETA, the scalar multiplier for C.
!
!    Input, real ( kind = 8 ) C(LDC,N).
!    Before entry, the leading M by N part of this array must contain the 
!    matrix C, except when BETA is zero, in which case C need not be set 
!    on entry.
!    On exit, the array C is overwritten by the M by N matrix
!      alpha * op ( A ) * op ( B ) + beta * C.
!
!    Input, integer ( kind = 4 ) LDC, the first dimension of C as declared in 
!    the calling routine.  max ( 1, M ) <= LDC.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) ldc

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(ldb,*)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(ldc,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncola
  integer ( kind = 4 ) nrowa
  integer ( kind = 4 ) nrowb
  logical nota
  logical notb
  real ( kind = 8 ) temp
  character * 1 transa
  character * 1 transb
!
!  Set NOTA and NOTB as true if A and B respectively are not
!  transposed and set NROWA, NCOLA and NROWB as the number of rows
!  and columns of A and the number of rows of B respectively.
!
  nota = ( transa == 'N' .or. transa == 'n' )

  if ( nota ) then
    nrowa = m
    ncola = k
  else
    nrowa = k
    ncola = m
  end if

  notb = ( transb == 'N' .or. transb == 'n' )

  if ( notb ) then
    nrowb = k
  else
    nrowb = n
  end if
!
!  Test the input parameters.
!
  info = 0

  if ( transa /= 'N' .and. transa /= 'n' .and. &
       transa /= 'C' .and. transa /= 'c' .and. &
       transa /= 'T' .and. transa /= 't' ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input TRANSA had illegal value.'
    stop 1
  end if

  if ( transb /= 'N' .and. transb /= 'n' .and. &
       transb /= 'C' .and. transb /= 'c' .and. &
       transb /= 'T' .and. transb /= 't' ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input TRANSB had illegal value.'
    stop 1
  end if

  if ( m < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input M had illegal value.'
    stop 1
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input N had illegal value.'
    stop 1
  end if

  if ( k < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input K had illegal value.'
    stop 1
  end if

  if ( lda < max ( 1, nrowa ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDA had illegal value.'
    stop 1
  end if

  if ( ldb < max ( 1, nrowb ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDB had illegal value.'
    stop 1
  end if

  if ( ldc < max ( 1, m ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDC had illegal value.'
    stop 1
  end if
!
!  Quick return if possible.
!
  if ( m == 0 ) then
    return
  end if

  if ( n == 0 ) then
    return
  end if

  if ( ( alpha == 0.0D+00 .or. k == 0 ) .and. beta == 1.0D+00 ) then
    return
  end if
!
!  And if alpha is zero.
!
  if ( alpha == 0.0D+00 ) then
    if ( beta == 0.0D+00 ) then
      c(1:m,1:n) = 0.0D+00
    else
      c(1:m,1:n) = beta * c(1:m,1:n)
    end if
    return
  end if
!
!  Start the operations.
!
  if ( notb ) then
!
!  Form  C := alpha*A*B + beta*C.
!
    if ( nota ) then

      do j = 1, n

        if ( beta == 0.0D+00 ) then
          c(1:m,j) = 0.0D+00
        else if ( beta /= 1.0D+00 ) then
          c(1:m,j) = beta * c(1:m,j)
        end if

        do l = 1, k
          if ( b(l,j) /= 0.0D+00 ) then
            c(1:m,j) = c(1:m,j) + alpha * b(l,j) * a(1:m,l)
          end if
        end do

      end do
!
!  Form  C := alpha*A'*B + beta*C
!
    else

      do j = 1, n
        do i = 1, m

          temp = dot_product ( a(1:k,i), b(1:k,j) )

          if ( beta == 0.0D+00 ) then
            c(i,j) = alpha * temp
          else
            c(i,j) = alpha * temp + beta * c(i,j)
          end if

        end do
      end do

    end if
!
!  Form  C := alpha*A*B' + beta*C
!
  else

    if ( nota ) then

      do j = 1, n

        if ( beta == 0.0D+00 ) then
          c(1:m,j) = 0.0D+00
        else if ( beta /= 1.0D+00 ) then
          c(1:m,j) = beta * c(1:m,j)
        end if

        do l = 1, k
          if ( b(j,l) /= 0.0D+00 ) then
            c(1:m,j) = c(1:m,j) + alpha * b(j,l) * a(1:m,l)
          end if
        end do

      end do
!
!  Form  C := alpha*A'*B' + beta*C
!
    else

      do j = 1, n
        do i = 1, m

          temp = dot_product ( a(1:k,i), b(j,1:k) )
          if ( beta == 0.0D+00 ) then
            c(i,j) = alpha * temp
          else
            c(i,j) = alpha * temp + beta * c(i,j)
          end if
        end do
      end do

    end if

  end if

  return
end
subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DGEMV computes y := alpha * A * x + beta * y for general matrix A.
!
!  Discussion:
!
!    DGEMV performs one of the matrix-vector operations
!      y := alpha*A *x + beta*y
!    or
!      y := alpha*A'*x + beta*y,
!    where alpha and beta are scalars, x and y are vectors and A is an
!    m by n matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character TRANS, specifies the operation to be performed:
!    'N' or 'N'   y := alpha*A *x + beta*y.
!    'T' or 'T'   y := alpha*A'*x + beta*y.
!    'C' or 'C'   y := alpha*A'*x + beta*y.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!    0 <= N.
!
!    Input, real ( kind = 8 ) ALPHA, the scalar multiplier for A * x.
!
!    Input, real ( kind = 8 ) A(LDA,N).  The M x N subarray contains
!    the matrix A.
!
!    Input, integer ( kind = 4 ) LDA, the the first dimension of A as declared
!    in the calling routine.  max ( 1, M ) <= LDA.
!
!    Input, real ( kind = 8 ) X(*), an array containing the vector to be 
!    multiplied by the matrix A.  
!    If TRANS = 'N' or 'n', then X must contain N entries, stored in INCX 
!    increments in a space of at least ( 1 + ( N - 1 ) * abs ( INCX ) ) 
!    locations.
!    Otherwise, X must contain M entries, store in INCX increments
!    in a space of at least ( 1 + ( M - 1 ) * abs ( INCX ) ) locations.
!
!    Input, integer ( kind = 4 ) INCX, the increment for the elements of
!    X.  INCX must not be zero.
!
!    Input, real ( kind = 8 ) BETA, the scalar multiplier for Y.
!
!    Input/output, real ( kind = 8 ) Y(*), an array containing the vector to
!    be scaled and incremented by A*X.
!    If TRANS = 'N' or 'n', then Y must contain M entries, stored in INCY
!    increments in a space of at least ( 1 + ( M - 1 ) * abs ( INCY ) ) 
!    locations.
!    Otherwise, Y must contain N entries, store in INCY increments
!    in a space of at least ( 1 + ( N - 1 ) * abs ( INCY ) ) locations.
!
!    Input, integer ( kind = 4 ) INCY, the increment for the elements of
!    Y.  INCY must not be zero.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jx
  integer ( kind = 4 ) jy
  integer ( kind = 4 ) kx
  integer ( kind = 4 ) ky
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) leny
  logical lsame
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  character trans
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( trans, 'N' ) .and.  &
       .not. lsame ( trans, 'T' ) .and.  &
       .not. lsame ( trans, 'C' ) ) then
    info = 1
  else if ( m < 0 ) then
    info = 2
  else if ( n < 0 ) then
    info = 3
  else if ( lda < max ( 1, m ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  else if ( incy == 0 ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgemv', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. &
       ( n == 0 ) .or. &
       ( ( alpha == 0.0D+00 ) .and. ( beta == 1.0D+00 ) ) ) then
   return
  end if
!
!  Set LENX and LENY, the lengths of the vectors x and y, and set
!  up the start points in X and Y.
!
  if ( lsame ( trans, 'N' ) ) then
    lenx = n
    leny = m
  else
    lenx = m
    leny = n
  end if

  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( lenx - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( leny - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
!  First form  y := beta*y.
!
  if ( beta /= 1.0D+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0D+00 ) then
        y(1:leny) = 0.0D+00
      else
        y(1:leny) = beta * y(1:leny)
      end if
    else
      iy = ky
      if ( beta == 0.0D+00 ) then
        do i = 1, leny
          y(iy) = 0.0D+00
          iy = iy + incy
        end do
      else
        do i = 1, leny
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0D+00 ) then
    return
  end if
!
!  Form y := alpha*A*x + y.
!
  if ( lsame ( trans, 'N' ) ) then
    jx = kx
    if ( incy == 1 ) then
      do j = 1, n
        if ( x(jx) /= 0.0D+00 ) then
          temp = alpha * x(jx)
          do i = 1, m
            y(i) = y(i) + temp * a(i,j)
          end do
        end if
        jx = jx + incx
      end do
    else
      do j = 1, n
        if ( x(jx) /= 0.0D+00 ) then
          temp = alpha * x(jx)
          iy = ky
          do i = 1, m
            y(iy) = y(iy) + temp * a(i,j)
            iy = iy + incy
          end do
        end if
        jx = jx + incx
      end do
    end if
!
!  Form y := alpha*A'*x + y.
!
  else
    jy = ky
    if ( incx == 1 ) then
      do j = 1, n
        temp = 0.0D+00
        do i = 1, m
          temp = temp + a(i,j) * x(i)
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    else
      do j = 1, n
        temp = 0.0D+00
        ix = kx
        do i = 1, m
          temp = temp + a(i,j) * x(ix)
          ix = ix + incx
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    end if
  end if

  return
end
subroutine dgeqp3 ( m, n, a, lda, jpvt, tau, work, lwork, info )

!*****************************************************************************80
!
!! DGEQP3 computes a QR factorization with column pivoting using Level 3 BLAS.
!
!  Discussion:
!
!    The factorization has the form 
!      A*P = Q * R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On input, the M by N matrix A.
!    On exit, the upper triangle of the array contains the
!    min(M,N)-by-N upper trapezoidal matrix R; the elements below
!    the diagonal, together with the array TAU, represent the
!    orthogonal matrix Q as a product of min(M,N) elementary
!    reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, M ) <= LDA.
!
!    Input/output, integer ( kind = 4 ) JPVT(N), the column pivot vector.
!    On entry, if jpvt(J) /= 0, the J-th column of A is permuted
!    to the front of A*P (a leading column); if jpvt(J)=0,
!    the J-th column of A is a free column.
!    On exit, if jpvt(J)=K, then the J-th column of A*P was the
!    the K-th column of A.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Input/output, real ( kind = 8 ) WORK(max(1,LWORK)), workspace.
!    On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!    Input/output, integer ( kind = 4 ) LWORK, the dimension of the array WORK. 
!    3 * N + 1 <= LWORK.  For optimal performance, 2*N+(N+1)*NB <= LWORK, 
!    where NB is the optimal blocksize.  However, if LWORK is -1, a workspace
!    query is assumed.  The routine only calculates the optimal size of the 
!    WORK array, returns this value as the first entry of the WORK array, 
!    and no error message related to LWORK is issued by XERBLA.
!
!    Output, integer ( kind = 4 ) INFO, return flag.
!    0: successful exit
!    < 0: if (-INFO)-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) fjb
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 ), parameter :: inb = 1
  integer ( kind = 4 ), parameter :: inbmin = 2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iws
  integer ( kind = 4 ), parameter :: ixover = 3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jpvt(*)
  logical lquery
  integer ( kind = 4 ) lwkopt
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) m
  integer ( kind = 4 ) minmn
  integer ( kind = 4 ) minws
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbmin
  integer ( kind = 4 ) nfxd
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) sm
  integer ( kind = 4 ) sminmn
  integer ( kind = 4 ) sn
  real ( kind = 8 ) tau(*)
  integer ( kind = 4 ) topbmn
  real ( kind = 8 ) work(*)
!
!  Test input arguments.
!
  info = 0
  lquery = ( lwork == -1 )

  if ( m < 0 ) then
    info = -1
  else if ( n < 0 ) then
    info = -2
  else if ( lda < max ( 1, m ) ) then
    info = -4
  end if

  if ( info == 0 ) then

    minmn = min ( m, n )

    if ( minmn == 0 ) then
      iws = 1
      lwkopt = 1
    else
      iws = 3 * n + 1
      nb = ilaenv ( inb, 'DGEQRF', ' ', m, n, -1, -1 )
      lwkopt = 2 * n + ( n + 1 ) * nb
    end if

    work(1) = lwkopt

    if ( ( lwork < iws ) .and. .not. lquery ) then
      info = -8
    end if

  end if

  if ( info /= 0 ) then
    call xerbla ( 'DGEQP3', -info )
    return
  else if ( lquery ) then
    return
  end if
!
!  Quick return if possible.
!
  if ( minmn == 0 ) then
    return
  end if
!
!  Move initial columns up front.
!
  nfxd = 1
  do j = 1, n
    if ( jpvt(j) /= 0 ) then
      if ( j /= nfxd ) then
        call dswap ( m, a(1,j), 1, a(1,nfxd), 1 )
        jpvt(j) = jpvt( nfxd )
        jpvt(nfxd) = j
      else
        jpvt(j) = j
      end if
      nfxd = nfxd + 1
    else
      jpvt(j) = j
    end if
  end do

  nfxd = nfxd - 1
!
!  Factorize fixed columns.
!
!  Compute the QR factorization of fixed columns and update
!  remaining columns.
!
  if ( 0 < nfxd ) then
    na = min ( m, nfxd )
    call dgeqrf ( m, na, a, lda, tau, work, lwork, info )
    iws = max ( iws, int ( work(1) ) )
    if ( na < n ) then
      call dormqr ( 'Left', 'Transpose', m, n - na, na, a, lda, tau, &
        a(1,na+1), lda, work, lwork, info )
      iws = max ( iws, int ( work(1) ) )
    end if
  end if
!
!  Factorize free columns
!
  if ( nfxd < minmn ) then

    sm = m - nfxd
    sn = n - nfxd
    sminmn = minmn - nfxd
!
!  Determine the block size.
!
    nb = ilaenv ( inb, 'DGEQRF', ' ', sm, sn, -1, -1 )
    nbmin = 2
    nx = 0

    if ( ( 1 < nb ) .and. ( nb < sminmn ) ) then
!
!  Determine when to cross over from blocked to unblocked code.
!
      nx = max ( 0, ilaenv ( ixover, 'DGEQRF', ' ', sm, sn, -1, -1 ) )
!
!  Determine if workspace is large enough for blocked code.
!
      if ( nx < sminmn ) then

        minws = 2 * sn + ( sn + 1 ) * nb
        iws = max ( iws, minws )

        if ( lwork < minws ) then
!
!  Not enough workspace to use optimal nb: Reduce nb and
!  determine the minimum value of nb.
!
          nb = ( lwork - 2 * sn ) / ( sn + 1 )
          nbmin = max ( 2, ilaenv ( inbmin, 'DGEQRF', ' ', sm, sn, -1, -1 ) )
        end if
      end if
    end if
!
!  Initialize partial column norms. The first N elements of work
!  store the exact column norms.
!
    do j = nfxd + 1, n
      work(j) = dnrm2 ( sm, a( nfxd+1,j), 1 )
      work(n+j) = work(j)
    end do

    if ( nbmin <= nb .and. nb < sminmn .and. nx < sminmn ) then
!
!  Use blocked code initially.
!
      j = nfxd + 1
!
!  Compute factorization: while loop.
!
      topbmn = minmn - nx

      do while ( j <= topbmn )

        jb = min ( nb, topbmn - j + 1 )
!
!  Factorize JB columns among columns J:N.
!
        call dlaqps ( m, n-j+1, j-1, jb, fjb, a(1,j), lda, jpvt(j), tau(j), &
          work(j), work(n+j), work(2*n+1), work(2*n+jb+1), n-j+1 )

        j = j + fjb

      end do

    else

      j = nfxd + 1

    end if
!
!  Use unblocked code to factor the last or only block.
!
    if ( j <= minmn ) then
      call dlaqp2 ( m, n-j+1, j-1, a(1,j), lda, jpvt(j), &
        tau(j), work(j), work(n+j), work(2*n+1) )
    end if

  end if

  work(1) = iws

  return
end
subroutine dgeqrf ( m, n, a, lda, tau, work, lwork, info )

!*****************************************************************************80
!
!! DGEQRF computes a QR factorization of a real M by N matrix A = Q * R.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On input, the M by N matrix A.
!    On output, the elements on and above the diagonal of the array contain the
!    min(M,N) by N upper trapezoidal matrix R (R is upper triangular if N <= M);
!    the elements below the diagonal, with the array TAU, represent the 
!    orthogonal matrix Q as a product of elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, M ) <= LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Input/output, real ( kind = 8 ) WORK(max(1,LWORK)), workspace.
!    On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!    Input/output, integer ( kind = 4 ) LWORK, the dimension of the array WORK. 
!    max ( 1, N ) <= LWORK.  For optimum performance N * NB <= LWORK, 
!    where NB is the optimal blocksize.
!    If LWORK = -1 on input, then a workspace query is assumed; the routine
!    only calculates the optimal size of the WORK array, returns
!    this value as the first entry of the WORK array, and no error
!    message related to LWORK is issued by XERBLA.
!
!    Output, integer ( kind = 4 ) INFO, return flag.
!    0: successful exit
!    < 0: if (-INFO)-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) iinfo
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iws
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ldwork
  logical lquery
  integer ( kind = 4 ) lwkopt
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbmin
  integer ( kind = 4 ) nx
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
!
!  Test the input arguments
!
  info = 0
  nb = ilaenv ( 1, 'DGEQRF', ' ', m, n, -1, -1 )
  lwkopt = n * nb
  work(1) = lwkopt
  lquery = ( lwork == -1 )

  if ( m < 0 ) then
    info = -1
  else if ( n < 0 ) then
    info = -2
  else if ( lda < max ( 1, m ) ) then
    info = -4
  else if ( lwork < max ( 1, n ) .and. .not. lquery ) then
    info = -7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'DGEQRF', -info )
    return
  else if ( lquery ) then
    return
  end if
!
!  Quick return if possible
!
  k = min ( m, n )
  if ( k == 0 ) then
    work(1) = 1
    return
  end if

  nbmin = 2
  nx = 0
  iws = n
!
!  Determine when to cross over from blocked to unblocked code.
!
  if ( 1 < nb .and. nb < k ) then

    nx = max ( 0, ilaenv ( 3, 'DGEQRF', ' ', m, n, -1, -1 ) )
!
!  Determine if workspace is large enough for blocked code.
!
    if ( nx < k ) then

      ldwork = n
      iws = ldwork * nb
!
!  Not enough workspace to use optimal NB:  reduce NB and
!  determine the minimum value of NB.
!
      if ( lwork < iws ) then
        nb = lwork / ldwork
        nbmin = max ( 2, ilaenv ( 2, 'dgeqrf', ' ', m, n, -1, -1 ) )
      end if

    end if

  end if
!
!  Use blocked code initially.
!
  if ( nbmin <= nb .and. nb < k .and. nx < k ) then

    do i = 1, k - nx, nb

      ib = min ( k-i+1, nb )
!
!  Compute the QR factorization of the current block A(i:m,i:i+ib-1)
!
      call dgeqr2 ( m-i+1, ib, a(i,i), lda, tau(i), work, iinfo )
!
!  Form the triangular factor of the block reflector
!  H = H(i) H(i+1) . . . H(i+ib-1)
!
      if ( i + ib <= n ) then

        call dlarft ( 'Forward', 'Columnwise', m-i+1, ib, &
          a(i,i), lda, tau(i), work, ldwork )
!
!  Apply H' to A(i:m,i+ib:n) from the left
!
        call dlarfb ( 'Left', 'Transpose', 'Forward', 'Columnwise', m-i+1, &
          n-i-ib+1, ib, a(i,i), lda, work, ldwork, a(i,i+ib), &
          lda, work(ib+1), ldwork )

      end if

    end do

  else

    i = 1

  end if
!
!  Use unblocked code to factor the last or only block.
!
  if ( i <= k ) then
    call dgeqr2( m-i+1, n-i+1, a(i,i), lda, tau(i), work, iinfo )
  end if

  work(1) = iws

  return
end
subroutine dgeqr2 ( m, n, a, lda, tau, work, info )

!*****************************************************************************80
!
!! DGEQR2 computes a QR factorization of a real M by N matrix A = Q * R.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On input, the M by N matrix A.
!    On output, the elements on and above the diagonal of the array contain the
!    min(M,N) by N upper trapezoidal matrix R (R is upper triangular if N <= M);
!    the elements below the diagonal, with the array TAU, represent the 
!    orthogonal matrix Q as a product of elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, M ) <= LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Input/output, real ( kind = 8 ) WORK(N), workspace.
!
!    Output, integer ( kind = 4 ) INFO, return flag.
!    0: successful exit
!    < 0: if (-INFO)-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) aii
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
!
!  Test the input arguments
!
  info = 0
  if ( m < 0 ) then
    info = -1
  else if ( n < 0 ) then
    info = -2
  else if ( lda < max ( 1, m ) ) then
    info = -4
  end if

  if ( info /= 0 ) then
    call xerbla ( 'DGEQR2', -info )
    return
  end if

  k = min ( m, n )
!
!  Generate elementary reflector H(i) to annihilate A(i+1:m,i).
!
  do i = 1, k

    call dlarfg ( m-i+1, a(i,i), a(min(i+1,m),i), 1, tau(i) )
!
!  Apply H(i) to A(i:m,i+1:n) from the left.
!
    if ( i < n ) then
      aii = a(i,i)
      a(i,i) = 1.0D+00
      call dlarf ( 'Left', m-i+1, n-i, a(i,i), 1, tau(i), &
        a(i,i+1), lda, work )
      a(i,i) = aii
    end if

  end do

  return
end
subroutine dger ( m, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! DGER computes A := alpha*x*y' + A.
!
!  Discussion:
!
!    DGER performs the rank 1 operation
!      A := alpha*x*y' + A,
!    where alpha is a scalar, x is an m element vector, y is an n element
!    vector and A is an m by n matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!    0 <= N.
!
!    Input, real ( kind = 8 ) ALPHA, the scalar alpha.
!
!    Input, real ( kind = 8 ) X(1+(m-1)*abs(INCX)).
!    The incremented array X must contain the m element vector x.
!
!    Input, integer ( kind = 4 ) INCX, the increment for the elements of
!    X.  INCX must not be zero.
!
!    Input,  real ( kind = 8 ) Y(1+(n-1)*abs(INCY) ).
!    The incremented array Y must contain the n element vector y.
!
!    Input, integer ( kind = 4 ) INCY, the increment for the elements of
!    Y.  INCY must not be zero.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    Before entry, the leading m by n part of the array A must
!    contain the matrix of coefficients. On exit, A is
!    overwritten by the updated matrix.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of A as declared
!    in the calling program. max ( 1, M ) <= LDA.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) temp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jy
  integer ( kind = 4 ) kx
!
!  Test the input parameters.
!
  info = 0
  if ( m < 0 ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  else if ( incy == 0 ) then
    info = 7
  else if ( lda < max ( 1, m ) ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dger  ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. ( n == 0 ) .or. ( alpha == 0.0D+00 ) ) then
    return
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( 0 < incy ) then
    jy = 1
  else
    jy = 1 - ( n - 1 ) * incy
  end if

  if ( incx == 1 ) then
    do j = 1, n
      if ( y(jy) /= 0.0D+00 ) then
        temp = alpha * y(jy)
        do i = 1, m
          a(i,j) = a(i,j) + x(i) * temp
        end do
      end if
      jy = jy + incy
    end do
  else
    if ( 0 < incx ) then
      kx = 1
    else
      kx = 1 - ( m - 1 ) * incx
    end if
    do j = 1, n
      if ( y(jy) /= 0.0D+00 ) then
        temp = alpha * y(jy)
        ix = kx
        do i = 1, m
          a(i,j) = a(i,j) + x(ix) * temp
          ix = ix + incx
        end do
      end if
      jy = jy + incy
    end do
  end if

  return
end
function dlamc3 ( a, b )

!*****************************************************************************80
!
!! DLAMC3 forces A and B to be stored prior to doing the addition.
!
!  Discussion:
!
!    This function is used to try to avoid situations where the compiler
!    optimization makes certain software measurements inaccurate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, two values to be added.
!
!    Output, real ( kind = 8 ) DLAMC3, the sum of A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) dlamc3

  dlamc3 = a + b

  return
end
function dlamch ( cmach )

!*****************************************************************************80
!
!! DLAMCH determines real ( kind = 8 ) machine parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character CMACH, specifies the value to be returned:
!    'E' or 'e',   DLAMCH := eps
!    'S' or 's ,   DLAMCH := sfmin
!    'B' or 'b',   DLAMCH := base
!    'P' or 'p',   DLAMCH := eps*base
!    'N' or 'n',   DLAMCH := t
!    'R' or 'r',   DLAMCH := rnd
!    'M' or 'm',   DLAMCH := emin
!    'U' or 'u',   DLAMCH := rmin
!    'L' or 'l',   DLAMCH := emax
!    'O' or 'o',   DLAMCH := rmax
!    where
!    eps   = relative machine precision
!    sfmin = safe minimum, such that 1/sfmin does not overflow
!    base  = base of the machine
!    prec  = eps*base
!    t     = number of (base) digits in the mantissa
!    rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!    emin  = minimum exponent before (gradual) underflow
!    rmin  = underflow threshold - base^(emin-1)
!    emax  = largest exponent before overflow
!
!    Output, real ( kind = 8 ) DLAMCH, the requested value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  character cmach
  real ( kind = 8 ) eps
  real ( kind = 8 ) dlamch
  logical lsame
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) rnd
  real ( kind = 8 ) sfmin
  real ( kind = 8 ) small
  real ( kind = 8 ) rmach
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Assume rounding, not chopping. Always.
!
  rnd = one

  if ( rnd == one ) then
    eps = epsilon ( zero ) * 0.5D+00
  else
    eps = epsilon ( zero )
  end if

  if ( lsame ( cmach, 'E' ) ) then
    rmach = eps
  else if ( lsame ( cmach, 'S' ) ) then
    sfmin = tiny ( zero )
    small = one / huge ( zero )
!
!  Use SMALL plus a bit, to avoid the possibility of rounding
!  causing overflow when computing  1/sfmin.
!
    if ( sfmin <= small ) then
      sfmin = small * ( one + eps )
    end if
    rmach = sfmin
  else if ( lsame ( cmach, 'B' ) ) then
    rmach = radix ( zero )
  else if ( lsame ( cmach, 'P' ) ) then
    rmach = eps * radix ( zero )
  else if ( lsame ( cmach, 'N' ) ) then
    rmach = digits ( zero )
  else if ( lsame ( cmach, 'R' ) ) then
    rmach = rnd
  else if ( lsame ( cmach, 'M' ) ) then
    rmach = minexponent ( zero )
  else if ( lsame ( cmach, 'U' ) ) then
    rmach = tiny ( zero )
  else if ( lsame ( cmach, 'L' ) ) then
    rmach = maxexponent ( zero )
  else if ( lsame ( cmach, 'O' ) ) then
    rmach = huge ( zero )
  else
    rmach = zero
  end if

  dlamch = rmach

  return
end
function dlapy2 ( x, y )

!*****************************************************************************80
!
!! DLAPY2 returns sqrt ( X^2 + Y^2 ), while avoiding unnecessary overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, two values.
!
!    Output, real ( kind = 8 ) DLAPY2, the Euclidean norm of (X,Y).
!
  implicit none

  real ( kind = 8 ) dlapy2
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xabs
  real ( kind = 8 ) y
  real ( kind = 8 ) yabs
  real ( kind = 8 ) z

  xabs = abs ( x )
  yabs = abs ( y )
  w = max ( xabs, yabs )
  z = min ( xabs, yabs )

  if ( z == 0.0D+00 ) then
    dlapy2 = w
  else
    dlapy2 = w * sqrt ( 1.0D+00 + ( z / w ) * ( z / w ) )
  end if

  return
end
subroutine dlaqp2 ( m, n, offset, a, lda, jpvt, tau, vn1, vn2, work )

!*****************************************************************************80
!
!! DLAQP2 QR factors with column pivoting the block A(OFFSET+1:M,1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A. 
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A. 
!    0 <= N.
!
!    Input, integer ( kind = 4 ) OFFSET, the number of rows of the matrix A 
!    that must be pivoted but not factorized. 0 <= OFFSET.
!
!    Input/output, real ( kind = 8 ) A(LDA,M).
!    On entry, the M-by-N matrix A.
!    On exit, the upper triangle of block A(OFFSET+1:M,1:N) is 
!    the triangular factor obtained; the elements in block
!    A(OFFSET+1:M,1:N) below the diagonal, together with the
!    array TAU, represent the orthogonal matrix Q as a product of
!    elementary reflectors. Block A(1:OFFSET,1:N) has been
!    accordingly pivoted, but not factorized.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A. 
!    max ( 1, M ) <= LDA.
!
!    Input/output, integer ( kind = 4 ) JPVT(N).
!    On entry, if JPVT(i)  /=  0, the i-th column of A is permuted
!    to the front of A*P (a leading column); if JPVT(i) = 0,
!    the i-th column of A is a free column.
!    On exit, if JPVT(i) = k, then the i-th column of A*P
!    was the k-th column of A.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Input/output, real ( kind = 8 ) VN1(N), the partial column norms.
!
!    Input/output, real ( kind = 8 ) VN2(N), the exact column norms.
!
!    Input/output, real ( kind = 8 ) WORK(N), workspace.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) aii
  real ( kind = 8 ) dlamch
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) jpvt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) offpi
  integer ( kind = 4 ) offset
  real ( kind = 8 ), parameter :: one = 1.0D+00
  integer ( kind = 4 ) pvt
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) tol3z
  real ( kind = 8 ) vn1(*)
  real ( kind = 8 ) vn2(*)
  real ( kind = 8 ) work(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  mn = min ( m - offset, n )
  tol3z = sqrt ( dlamch ( 'Epsilon' ) )
!
!  Compute factorization.
!
  do i = 1, mn

    offpi = offset + i
!
!  Determine I-th pivot column and swap if necessary.
!
    pvt = ( i - 1 ) + idamax ( n-i+1, vn1(i), 1 )

    if ( pvt /= i ) then
      call dswap ( m, a(1,pvt), 1, a(1,i), 1 )
      itemp = jpvt( pvt )
      jpvt( pvt ) = jpvt(i)
      jpvt(i) = itemp
      vn1( pvt ) = vn1(i)
      vn2( pvt ) = vn2(i)
    end if
!
!  Generate elementary reflector H(i).
!
    if ( offpi < m ) then
      call dlarfg ( m-offpi+1, a(offpi,i), a( offpi+1,i), 1, tau(i) )
    else
      call dlarfg ( 1, a(m,i), a(m,i), 1, tau(i) )
    end if
!
!  Apply H(i)' to A(offset+i:m,i+1:n) from the left.
!
    if ( i < n ) then
      aii = a(offpi,i)
      a(offpi,i) = 1.0D+00
      call dlarf ( 'Left', m-offpi+1, n-i, a(offpi,i), 1, tau(i), &
        a(offpi,i+1), lda, work(1) )
      a(offpi,i) = aii
    end if
!
!  Update partial column norms.
!
    do j = i + 1, n
!
!  NOTE: The following 4 lines follow from the analysis in
!  Lapack Working Note 176.
!
      if ( vn1(j) /= 0.0D+00 ) then

        temp = 1.0D+00 - ( abs ( a(offpi,j) ) / vn1(j) )**2
        temp = max ( temp, 0.0D+00 )
        temp2 = temp * ( vn1(j) / vn2(j) )**2

        if ( temp2  <=  tol3z ) then
          if ( offpi < m ) then
            vn1(j) = dnrm2 ( m-offpi, a(offpi+1,j), 1 )
            vn2(j) = vn1(j)
          else
            vn1(j) = 0.0D+00
            vn2(j) = 0.0D+00
          end if
        else
          vn1(j) = vn1(j) * sqrt ( temp )
        end if

      end if

    end do

  end do

  return
end
subroutine dlaqps ( m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, &
  vn2, auxv, f, ldf )

!*****************************************************************************80
!
!! DLAQPS computes a step of QR factorization with column pivoting.
!
!  Discussion:
!
!    DLAQPS computes a step of QR factorization with column pivoting
!    of a real M-by-N matrix A by using Blas-3.  It tries to factorize
!    nb columns from A starting from the row offset+1, and updates all
!    of the matrix with Blas-3 xGEMM.
!
!    In some cases, due to catastrophic cancellations, it cannot
!    factorize nb columns.  Hence, the actual number of factorized
!    columns is returned in kb.
!
!    Block A(1:offset,1:N) is accordingly pivoted, but not factorized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A. 
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A. 
!    0 <= N.
!
!    Input, integer ( kind = 4 ) OFFSET, the number of rows of A that have 
!    been factorized in previous steps.
!
!    Input, integer ( kind = 4 ) NB, the number of columns to factorize.
!
!    Output, integer ( kind = 4 ) KB, the number of columns actually factorized.
!
!    Input/output, real ( kind = 8 ) A(LDA,N),
!    On entry, the M-by-N matrix A.
!    On exit, block A(offset+1:M,1:kb) is the triangular
!    factor obtained and block A(1:offset,1:N) has been
!    accordingly pivoted, but no factorized.
!    The rest of the matrix, block A(offset+1:M,kb+1:N) has been updated.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, M ) <= LDA.
!
!    Input/output, integer ( kind = 4 ) JPVT(N), the column pivot array.
!    jpvt(I) = K <==> Column K of the full matrix A has been
!    permuted into position I in AP.
!
!    Output, real ( kind = 8 ) TAU(KB), the scalar factors of the 
!    elementary reflectors.
!
!    Input/output, real ( kind = 8 ) VN1(N), the vector with the partial 
!    column norms.
!
!    Input/output, real ( kind = 8 ) VN2(N), the vector with the exact 
!    column norms.
!
!    Input/output, real ( kind = 8 ) AUXV(NB), auxiliary vector.
!
!    Input/output, real ( kind = 8 ) F(LDF,NB).
!    Matrix F' = L*Y'*A.
!
!    Input, integer ( kind = 4 ) LDF, the leading dimension of the array F.
!    max ( 1, N ) <= LDF.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldf

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) akk
  real ( kind = 8 ) auxv(*)
  real ( kind = 8 ) dlamch
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) f(ldf,*)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) lastrk
  integer ( kind = 4 ) lsticc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) offset
  real ( kind = 8 ), parameter :: one = 1.0D+00
  integer ( kind = 4 ) pvt
  integer ( kind = 4 ) rk
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) tol3z
  real ( kind = 8 ) vn1(*)
  real ( kind = 8 ) vn2(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  lastrk = min ( m, n + offset )
  lsticc = 0
  k = 0
  tol3z = sqrt ( dlamch ( 'Epsilon' ) )
!
!  Beginning of while loop.
!
  do while ( ( k < nb ) .and. ( lsticc == 0 ) )

    k = k + 1
    rk = offset + k
!
!  Determine ith pivot column and swap if necessary
!
    pvt = ( k - 1 ) + idamax ( n-k+1, vn1(k), 1 )

    if ( pvt /= k ) then
      call dswap ( m, a(1,pvt), 1, a(1,k), 1 )
      call dswap ( k-1, f(pvt,1), ldf, f(k,1), ldf )
      itemp = jpvt(pvt)
      jpvt(pvt) = jpvt(k)
      jpvt(k) = itemp
      vn1(pvt) = vn1(k)
      vn2(pvt) = vn2(k)
    end if
!
!  Apply previous Householder reflectors to column K:
!  A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
!
    if ( 1 < k ) then
      call dgemv ( 'No transpose', m-rk+1, k-1, -one, a(rk,1), &
        lda, f(k,1), ldf, one, a(rk,k), 1 )
    end if
!
!  Generate elementary reflector H(k).
!
    if ( rk < m ) then
      call dlarfg ( m-rk+1, a(rk,k), a(rk+1,k), 1, tau(k) )
    else
      call dlarfg ( 1, a(rk,k), a(rk,k), 1, tau(k) )
    end if

     akk = a(rk,k)
     a(rk,k) = 1.0D+00
!
!  Compute Kth column of F:
!  Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K).
!
    if ( k < n ) then
      call dgemv ( 'Transpose', m-rk+1, n-k, tau(k), a(rk,k+1 ), lda, &
        a(rk,k), 1, zero, f(k+1,k), 1 )
    end if
!
!  Padding F(1:K,K) with zeros.
!
    f(1:k,k) = 0.0D+00
!
!  Incremental updating of F:
!  F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)' * A(RK:M,K).
!
    if ( 1 < k ) then

      call dgemv ( 'Transpose', m-rk+1, k-1, -tau(k), a(rk,1), &
        lda, a(rk,k), 1, zero, auxv(1), 1 )

      call dgemv ( 'No transpose', n, k-1, one, f(1,1), ldf, &
        auxv(1), 1, one, f(1,k), 1 )

    end if
!
!  Update the current row of A:
!  A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.
!
    if ( k < n ) then
      call dgemv ( 'No transpose', n-k, k, -one, f(k+1,1), ldf, &
        a(rk,1), lda, one, a(rk,k+1), lda )
    end if
!
!  Update partial column norms.
!
    if ( rk < lastrk ) then
      do j = k + 1, n
        if ( vn1(j) /= 0.0D+00 ) then
!
!  NOTE: The following 4 lines follow from the analysis in
!  Lapack Working Note 176.
!
          temp = abs ( a(rk,j) ) / vn1(j)
          temp = max ( 0.0D+00, ( 1.0D+00 + temp ) * ( 1.0D+00 - temp ) )
          temp2 = temp * ( vn1(j) / vn2(j) )**2
          if ( temp2  <=  tol3z ) then
            vn2(j) = real ( lsticc, kind = 8 )
            lsticc = J
          else
            vn1(j) = vn1(j) * sqrt ( temp )
          end if
        end if
      end do
    end if

    a(rk,k) = akk

  end do

  kb = k
  rk = offset + kb
!
!  Apply the block reflector to the rest of the matrix:
!  A(offset+kb+1:M,kb+1:N) := A(offset+kb+1:M,kb+1:N) -
!  A(offset+kb+1:M,1:kb)*F(kb+1:N,1:kb)'.
!
  if ( kb < min ( n, m-offset ) ) then
    call dgemm ( 'No transpose', 'Transpose', m-rk, n-kb, kb, -one, &
      a(rk+1,1), lda, f(kb+1,1), ldf, one, a(rk+1,kb+1), lda )
  end if
!
!  Recomputation of difficult columns.
!
  do while ( 0 < lsticc )

    itemp = nint ( vn2( lsticc ) )
    vn1(lsticc) = dnrm2 ( m - rk, a(rk+1,lsticc), 1 )
!
!  NOTE: The computation of vn1( lsticc ) relies on the fact that 
!  SNRM2 does not fail on vectors with norm below the value of
!  SQRT(DLAMCH('S')) 
!
    vn2(lsticc) = vn1(lsticc)
    lsticc = itemp

  end do

  return
end
subroutine dlarf ( side, m, n, v, incv, tau, c, ldc, work )

!*****************************************************************************80
!
!! DLARF applies a real elementary reflector to an M by N matrix.
!
!  Discussion:
!
!    DLARF applies a real elementary reflector H to a real m by n matrix
!    C, from either the left or the right. H is represented in the form
!      H = I - tau * v * v'
!    where tau is a real scalar and v is a real vector.
!
!    If tau = 0, then H is taken to be the unit matrix.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE, specifies pre- or post-multiplication.
!    'L': form  H * C
!    'R': form  C * H
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix C.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix C.
!
!    Input, real ( kind = 8 ) V(1 + (M-1)*abs(INCV)) if SIDE = 'L', or
!    V(1 + (N-1)*abs(INCV)) if SIDE = 'R'.  This array contains the vector v 
!    in the representation of H.  V is not used if TAU = 0.
!
!    Input, integer ( kind = 4 ) INCV, the increment between elements of V. 
!    INCV must not be 0.
!
!    Input, real ( kind = 8 ) TAU, the value tau in the representation of H.
!
!    Input/output, real ( kind = 8 ) C(LDC,N).
!    On input, the m by n matrix C.
!    On output, C is overwritten by the matrix H * C if SIDE = 'L',
!    or C * H if SIDE = 'R'.
!
!    Input, integer ( kind = 4 ) LDC, the leading dimension of the array C. 
!    max ( 1, M ) <= LDC.
!
!    Input/output, real ( kind = 8 ) WORK(*), workspace.
!    Size must be at least:
!    * N if SIDE = 'L',
!    * M if SIDE = 'R'.
!
  implicit none

  integer ( kind = 4 ) ldc

  logical applyleft
  real ( kind = 8 ) c(ldc,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iladlc
  integer ( kind = 4 ) iladlr
  integer ( kind = 4 ) incv
  integer ( kind = 4 ) lastc
  integer ( kind = 4 ) lastv
  logical lsame
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  character side
  real ( kind = 8 ) tau
  real ( kind = 8 ) v(*)
  real ( kind = 8 ) work(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  applyleft = lsame ( side, 'L' )

  lastv = 0
  lastc = 0
!
!  Set up variables for scanning V.  LASTV begins pointing to the end of V.
!
  if ( tau /= 0.0D+00 ) then

    if ( applyleft ) then
      lastv = m
    else
      lastv = n
    end if

    if ( 0 < incv ) then
      i = 1 + ( lastv - 1 ) * incv
    else
      i = 1
    end if
!
!  Look for the last non-zero row in V.
!
    do while ( 0 < lastv .and. v(i) == 0.0D+00 )
      lastv = lastv - 1
      i = i - incv
    end do
!
!  Scan for the last non-zero column in C(1:lastv,:).
!
    if ( applyleft ) then

      lastc = iladlc ( lastv, n, c, ldc )
!
!  Scan for the last non-zero row in C(:,1:lastv).
!
    else

      lastc = iladlr ( m, lastv, c, ldc )

    end if

  end if
!
!  Note that LASTC == 0 renders the BLAS operations null; no special
!  case is needed at this level.
!
  if ( applyleft ) then
!
!  Form H * C
!
    if ( 0 < lastv ) then
!
!  w(1:lastc,1) := C(1:lastv,1:lastc)' * v(1:lastv,1)
!
      call dgemv ( 'Transpose', lastv, lastc, one, c, ldc, v, incv, &
        zero, work, 1 )
!
!  C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)'
!
      call dger ( lastv, lastc, -tau, v, incv, work, 1, c, ldc )

    end if

  else
!
!  Form  C * H
!
    if ( 0 < lastv ) then
!
!  w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
      call dgemv ( 'No transpose', lastc, lastv, one, c, ldc, &
        v, incv, zero, work, 1 )
!
!  C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)'
!
      call dger ( lastc, lastv, -tau, work, 1, v, incv, c, ldc )

    end if

  end if

  return
end
subroutine dlarfb ( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, &
  c, ldc, work, ldwork )

!*****************************************************************************80
!
!! DLARFB applies a block reflector to a general rectangular matrix.
!
!  Discussion:
!
!    DLARFB applies a real block reflector H or its transpose H' to a
!    real M by N matrix C, from either the left or the right.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE.
!    'L': apply H or H' from the Left
!    'R': apply H or H' from the Right
!
!    Input, character TRANS.
!    'N': apply H (No transpose)
!    'T': apply H' (Transpose)
!
!    Input, character DIRECT, indicates how H is formed from a product of 
!    elementary reflectors:
!    'F': H = H(1) H(2) . . . H(k) (Forward)
!    'B': H = H(k) . . . H(2) H(1) (Backward)
!
!    Input, character STOREV, indicates how the vectors which define the 
!    elementary reflectors are stored:
!    'C': Columnwise
!    'R': Rowwise
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix C.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix C.
!
!    Input, integer ( kind = 4 ) K, the order of the matrix T (= the number 
!    of elementary reflectors whose product defines the block reflector).
!
!    Input, real ( kind = 8 ) V(LDV,K) if STOREV = 'C',
!    (LDV,M) if STOREV = 'R' and SIDE = 'L'
!    V(LDV,N) if STOREV = 'R' and SIDE = 'R'.
!    The matrix V.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of the array V.
!    If STOREV = 'C' and SIDE = 'L', max(1,M) <= LDV;
!    if STOREV = 'C' and SIDE = 'R', max(1,N) <= LDV;
!    if STOREV = 'R', K <= LDV.
!
!    Input, real ( kind = 8 ) T(LDT,K), the triangular k by k matrix T in 
!    the representation of the block reflector.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of the array T. 
!    K <= LDT.
!
!    Input/output, real ( kind = 8 ) C(LDC,N).
!    On entry, the m by n matrix C.
!    On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!    Input, integer ( kind = 4 ) LDC, the leading dimension of the array C. 
!    max(1,M) <= LDC.
!
!    Input/output, real ( kind = 8 ) WORK(LDWORK,K), workspace.
!
!    Input, integer ( kind = 4 ) LDWORK, the leading dimension of WORK.
!    If SIDE = 'L', max(1,N) <= LDWORK;
!    if SIDE = 'R', max(1,M) <= LDWORK.
!
  implicit none

  integer ( kind = 4 ) ldc
  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) ldwork

  real ( kind = 8 ) c(ldc,*)
  character direct
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical lsame
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  character side
  character storev
  real ( kind = 8 ) t(ldt,*)
  character trans
  character transt
  real ( kind = 8 ) v(ldv,*)
  real ( kind = 8 ) work(ldwork,*)
!
!  Quick return if possible
!
  if ( m <= 0 .or. n <= 0 ) then
    return
  end if

  if ( lsame ( trans, 'N' ) ) then
    transt = 'T'
  else
    transt = 'N'
  end if

  if ( lsame ( storev, 'C' ) ) then

    if ( lsame ( direct, 'F' ) ) then
!
!  Let  V =  ( V1 )    (first K rows)
!            ( V2 )
!  where  V1  is unit lower triangular.
!
      if ( lsame ( side, 'L' ) ) then
!
!  Form  H * C  or  H' * C  where  C = ( C1 )
!                    ( C2 )
!
!  W := C' * V  =  (C1' * V1 + C2' * V2)  (stored in WORK)
!  W := C1'
!
        do j = 1, k
          call dcopy ( n, c(j,1), ldc, work(1,j), 1 )
        end do
!
!  W := W * V1
!
        call dtrmm ( 'Right', 'Lower', 'No transpose', 'Unit', n, &
          k, one, v, ldv, work, ldwork )
!
!  W := W + C2' * V2
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'No transpose', n, k, m-k, &
            one, c(k+1,1), ldc, v(k+1,1), ldv, one, work, ldwork )
        end if
!
!  W := W * T'  or  W * T
!
        call dtrmm ( 'Right', 'Upper', transt, 'Non-unit', n, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - V * W'
!  C2 := C2 - V2 * W'
!
        if ( k < m ) then
          call dgemm ( 'No transpose', 'Transpose', m-k, n, k, &
            -one, v(k+1,1), ldv, work, ldwork, one, c(k+1,1), ldc )
        end if
!
!  W := W * V1'
!
        call dtrmm ( 'Right', 'Lower', 'Transpose', 'Unit', n, k, &
          one, v, ldv, work, ldwork )
!
! C1 := C1 - W'
!
        do j = 1, k
          do i = 1, n
            c(j,i) = c(j,i) - work(i,j)
          end do
        end do

      else if ( lsame ( side, 'R' ) ) then
!
!  Form  C * H  or  C * H'  where  C = ( C1  C2 )
!  W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!  W := C1
!
        do j = 1, k
          call dcopy ( m, c(1,j), 1, work(1,j), 1 )
        end do
!
!  W := W * V1
!
        call dtrmm ( 'Right', 'Lower', 'No transpose', 'Unit', m, &
          k, one, v, ldv, work, ldwork )
!
!  W := W + C2 * V2
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'No transpose', m, k, n-k, &
            one, c(1,k+1), ldc, v(k+1,1), ldv, one, work, ldwork )
        end if
!
!  W := W * T  or  W * T'
!
        call dtrmm ( 'Right', 'Upper', trans, 'Non-unit', m, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - W * V'
!  C2 := C2 - W * V2'
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'Transpose', m, n-k, k, &
            -one, work, ldwork, v(k+1,1), ldv, one, c(1,k+1), ldc )
        end if
!
!  W := W * V1'
!
        call dtrmm ( 'Right', 'Lower', 'Transpose', 'Unit', m, k, &
          one, v, ldv, work, ldwork )
!
!  C1 := C1 - W
!
        do j = 1, k
          do i = 1, m
            c(i,j) = c(i,j) - work(i,j)
          end do
        end do

      end if

    else
!
!  Let  V =  ( V1 )
!            ( V2 )    (last K rows)
!  where  V2  is unit upper triangular.
!
      if ( lsame ( side, 'L' ) ) then
!
!  Form  H * C  or  H' * C  where  C = ( C1 )
!                                      ( C2 )
!
!  W := C' * V  =  (C1' * V1 + C2' * V2)  (stored in WORK)
!  W := C2'
!
        do j = 1, k
          call dcopy ( n, c(m-k+j,1), ldc, work(1,j), 1 )
        end do
!
!  W := W * V2
!
        call dtrmm ( 'Right', 'Upper', 'No transpose', 'Unit', n, &
          k, one, v(m-k+1,1), ldv, work, ldwork )
!
!  W := W + C1' * V1
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'No transpose', n, k, m-k, &
            one, c, ldc, v, ldv, one, work, ldwork )
        end if
!
!  W := W * T'  or  W * T
!
        call dtrmm ( 'Right', 'Lower', transt, 'Non-unit', n, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - V * W'
!  C1 := C1 - V1 * W'
!
        if ( k < m ) then
          call dgemm ( 'No transpose', 'Transpose', m-k, n, k, &
            -one, v, ldv, work, ldwork, one, c, ldc )
        end if
!
!  W := W * V2'
!
        call dtrmm ( 'Right', 'Upper', 'Transpose', 'Unit', n, k, &
          one, v(m-k+1,1), ldv, work, ldwork )
!
!  C2 := C2 - W'
!
        do j = 1, k
          do i = 1, n
            c(m-k+j,i) = c(m-k+j,i) - work(i,j)
          end do
        end do

      else if ( lsame ( side, 'R' ) ) then
!
!  Form  C * H  or  C * H'  where  C = ( C1  C2 )
!  W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!  W := C2
!
        do j = 1, k
          call dcopy ( m, c(1,n-k+j), 1, work(1,j), 1 )
        end do
!
!  W := W * V2
!
        call dtrmm ( 'Right', 'Upper', 'No transpose', 'Unit', m, &
          k, one, v(n-k+1,1), ldv, work, ldwork )
!
!  W := W + C1 * V1
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'No transpose', m, k, n-k, &
            one, c, ldc, v, ldv, one, work, ldwork )
        end if
!
!  W := W * T  or  W * T'
!
        call dtrmm ( 'Right', 'Lower', trans, 'Non-unit', m, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - W * V'
!  C1 := C1 - W * V1'
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'Transpose', m, n-k, k, &
            -one, work, ldwork, v, ldv, one, c, ldc )
        end if
!
!  W := W * V2'
!
        call dtrmm ( 'Right', 'Upper', 'Transpose', 'Unit', m, k, &
          one, v(n-k+1,1), ldv, work, ldwork )
!
!  C2 := C2 - W
!
        do j = 1, k
          do i = 1, m
            c(i,n-k+j) = c(i,n-k+j) - work(i,j)
          end do
        end do

      end if

    end if

  else if ( lsame ( storev, 'R' ) ) then

    if ( lsame ( direct, 'F' ) ) then
!
!  Let  V =  ( V1  V2 )    (V1: first K columns)
!  where  V1  is unit upper triangular.
!
      if ( lsame ( side, 'L' ) ) then
!
!  Form  H * C  or  H' * C  where  C = ( C1 )
!                                      ( C2 )
!  W := C' * V'  =  (C1' * V1' + C2' * V2') (stored in WORK)
!  W := C1'
!
        do j = 1, k
          call dcopy ( n, c(j,1), ldc, work(1,j), 1 )
        end do
!
!  W := W * V1'
!
        call dtrmm ( 'Right', 'Upper', 'Transpose', 'Unit', n, k, &
          one, v, ldv, work, ldwork )
!
!  W := W + C2' * V2'
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'Transpose', n, k, m-k, one, &
            c(k+1,1), ldc, v(1,k+1), ldv, one, work, ldwork )
        end if
!
!  W := W * T'  or  W * T
!
        call dtrmm ( 'Right', 'Upper', transt, 'Non-unit', n, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - V' * W'
!  C2 := C2 - V2' * W'
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'Transpose', m-k, n, k, -one, &
            v(1,k+1), ldv, work, ldwork, one, c(k+1,1), ldc )
        end if
!
!  W := W * V1
!
        call dtrmm ( 'Right', 'Upper', 'No transpose', 'Unit', n, &
          k, one, v, ldv, work, ldwork )
!
!  C1 := C1 - W'
!
        do j = 1, k
          do i = 1, n
            c(j,i) = c(j,i) - work(i,j)
          end do
        end do

      else if ( lsame ( side, 'R' ) ) then
!
!  Form  C * H  or  C * H'  where  C = ( C1  C2 )
!  W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!  W := C1
!
        do j = 1, k
          call dcopy ( m, c(1,j), 1, work(1,j), 1 )
        end do
!
!  W := W * V1'
!
        call dtrmm ( 'Right', 'Upper', 'Transpose', 'Unit', m, k, &
          one, v, ldv, work, ldwork )
!
!  W := W + C2 * V2'
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'Transpose', m, k, n-k, &
            one, c(1,k+1), ldc, v(1,k+1), ldv, one, work, ldwork )
        end if
!
!  W := W * T  or  W * T'
!
        call dtrmm ( 'Right', 'Upper', trans, 'Non-unit', m, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - W * V
!  C2 := C2 - W * V2
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'No transpose', m, n-k, k, &
            -one, work, ldwork, v(1,k+1), ldv, one, c(1,k+1), ldc )
        end if
!
!  W := W * V1
!
        call dtrmm ( 'Right', 'Upper', 'No transpose', 'Unit', m, &
          k, one, v, ldv, work, ldwork )
!
!  C1 := C1 - W
!
        c(1:m,1:k) = c(1:m,1:k) - work(1:m,1:k)

      end if

    else
!
!  Let  V =  ( V1  V2 )    (V2: last K columns)
!  where  V2  is unit lower triangular.
!
      if ( lsame ( side, 'L' ) ) then
!
!  Form  H * C  or  H' * C  where  C = ( C1 )
!                                      ( C2 )
!  W := C' * V'  =  (C1' * V1' + C2' * V2') (stored in WORK)
!  W := C2'
!
        do j = 1, k
          call dcopy ( n, c(m-k+j,1), ldc, work(1,j), 1 )
        end do
!
!  W := W * V2'
!
        call dtrmm ( 'Right', 'Lower', 'Transpose', 'Unit', n, k, &
          one, v(1,m-k+1), ldv, work, ldwork )
!
!  W := W + C1' * V1'
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'Transpose', n, k, m-k, one, &
            c, ldc, v, ldv, one, work, ldwork )
        end if
!
!  W := W * T'  or  W * T
!
        call dtrmm ( 'Right', 'Lower', transt, 'Non-unit', n, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - V' * W'
!  C1 := C1 - V1' * W'
!
        if ( k < m ) then
          call dgemm ( 'Transpose', 'Transpose', m-k, n, k, -one, &
            v, ldv, work, ldwork, one, c, ldc )
        end if
!
!  W := W * V2
!
        call dtrmm ( 'Right', 'Lower', 'No transpose', 'Unit', n, &
          k, one, v(1,m-k+1), ldv, work, ldwork )
!
! C2 := C2 - W'
!
        do j = 1, k
          do i = 1, n
            c(m-k+j,i) = c(m-k+j,i) - work(i,j)
          end do
        end do

      else if ( lsame ( side, 'R' ) ) then
!
!  Form  C * H  or  C * H'  where  C = ( C1  C2 )
!  W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!  W := C2
!
        do j = 1, k
          call dcopy ( m, c(1,n-k+j), 1, work(1,j), 1 )
        end do
!
!  W := W * V2'
!
        call dtrmm ( 'Right', 'Lower', 'Transpose', 'Unit', m, k, &
          one, v(1,n-k+1), ldv, work, ldwork )
!
!  W := W + C1 * V1'
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'Transpose', m, k, n-k, &
            one, c, ldc, v, ldv, one, work, ldwork )
        end if
!
!  W := W * T  or  W * T'
!
        call dtrmm ( 'Right', 'Lower', trans, 'Non-unit', m, k, &
          one, t, ldt, work, ldwork )
!
!  C := C - W * V
!  C1 := C1 - W * V1
!
        if ( k < n ) then
          call dgemm ( 'No transpose', 'No transpose', m, n-k, k, &
            -one, work, ldwork, v, ldv, one, c, ldc )
        end if
!
!  W := W * V2
!
        call dtrmm ( 'Right', 'Lower', 'No transpose', 'Unit', m, &
          k, one, v(1,n-k+1), ldv, work, ldwork )
!
!  C1 := C1 - W
!
        c(1:m,n-k+1:n) = c(1:m,n-k+1:n) - work(1:m,1:k)

      end if

    end if

  end if

  return
end
subroutine dlarfg ( n, alpha, x, incx, tau )

!*****************************************************************************80
!
!! DLARFG generates a real elementary reflector H of order N.
!
!  Discussion:
!
!    DLARFG generates a real elementary reflector H of order n, such that
!      H * ( alpha ) = ( beta ),   H' * H = I.
!          (   x   )   (   0  )
!    where alpha and beta are scalars, and x is an (n-1)-element real
!    vector.  H is represented in the form
!      H = I - tau * ( 1 ) * ( 1 v' ) ,
!                    ( v )
!    where tau is a real scalar and v is a real (n-1)-element vector.
!
!    If the elements of x are all zero, then tau = 0 and H is taken to be
!    the unit matrix.
!
!    Otherwise  1 <= tau <= 2.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the elementary reflector.
!
!    Input/output, real ( kind = 8 ) ALPHA.
!    On entry, the value alpha.
!    On exit, the value beta.
!
!    Input/output, real ( kind = 8 ) X(1+(N-2)*abs(INCX))
!    On entry, the vector x.
!    On exit, it is overwritten with the vector v.
!
!    Input, integer ( kind = 4 ) INCX, the increment between elements of X.
!    0 < INCX.
!
!    Output, real ( kind = 8 ) TAU, the value tau.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) dlamch
  real ( kind = 8 ) dlapy2
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) knt
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) rsafmn
  real ( kind = 8 ) safmin
  real ( kind = 8 ) tau
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xnorm
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  if ( n <= 1 ) then
    tau = 0.0D+00
    return
  end if

  xnorm = dnrm2 ( n - 1, x, incx )

  if ( xnorm == 0.0D+00 ) then
!
!  H = I.
!
    tau = 0.0D+00

  else
!
!  General case.
!
    beta = - sign ( dlapy2 ( alpha, xnorm ), alpha )
    safmin = dlamch ( 'S' ) / dlamch ( 'E' )
    knt = 0
!
!  XNORM, BETA may be inaccurate; scale X and recompute them.
!
    if ( abs ( beta ) < safmin ) then

      rsafmn = 1.0D+00 / safmin

      do
        knt = knt + 1
        call dscal ( n - 1, rsafmn, x, incx )
        beta = beta * rsafmn
        alpha = alpha * rsafmn
        if ( safmin <= abs ( beta ) ) then
          exit
        end if
      end do
!
!  New BETA is at most 1, at least SAFMIN.
!
      xnorm = dnrm2 ( n-1, x, incx )
      beta = - sign ( dlapy2 ( alpha, xnorm ), alpha )

    end if

    tau = ( beta - alpha ) / beta
    call dscal ( n - 1, one / ( alpha - beta ), x, incx )
!
!  If ALPHA is subnormal, it may lose relative accuracy.
!
    do j = 1, knt
      beta = beta * safmin
    end do
    alpha = beta
  end if

  return
end
subroutine dlarft ( direct, storev, n, k, v, ldv, tau, t, ldt )

!*****************************************************************************80
!
!! DLARFT forms the triangular factor T of a real block reflector H of order N.
!
!  Discussion:
!
!    DLARFT forms the triangular factor T of a real block reflector H
!    of order N, which is defined as a product of K elementary reflectors.
!
!    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!    If STOREV = 'C', the vector which defines the elementary reflector
!    H(i) is stored in the I-th column of the array V, and
!      H  =  I - V * T * V'
!
!    If STOREV = 'R', the vector which defines the elementary reflector
!    H(i) is stored in the I-th row of the array V, and
!      H  =  I - V' * T * V
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 September 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character DIRECT, specifies the order in which the elementary 
!    reflectors are multiplied to form the block reflector:
!    'F': H = H(1) H(2) . . . H(k) (Forward)
!    'B': H = H(k) . . . H(2) H(1) (Backward)
!
!    Input, character STOREV, specifies how the vectors which define the 
!    elementary reflectors are stored:
!    'C': columnwise
!    'R': rowwise
!
!    Input, integer ( kind = 4 ) N, the order of the block reflector H. 
!    0 <= N.
!
!    Input, integer ( kind = 4 ) K, the order of the triangular factor 
!    T (= the number of elementary reflectors). 1 <= K.
!
!    Input, real ( kind = 8 ) V(*,*).
!    * V(LDV,K) if STOREV = 'C', or 
!    * V(LDV,N) if STOREV = 'R'.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of the array V.
!    If STOREV = 'C', max(1,N) <= LDV; if STOREV = 'R', K <= LDV.
!
!    Input, real ( kind = 8 ) TAU(K).  TAU(i) must contain the scalar factor 
!    of the elementary reflector H(i).
!
!    Output, real ( kind = 8 ) T(LDT,K), the k by k triangular factor T of 
!    the block reflector.  If DIRECT = 'F', T is upper triangular; if 
!    DIRECT = 'B', T is lower triangular. 
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of the array T. 
!    K <= LDT.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) ldv

  character direct
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lastv
  logical lsame
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  integer ( kind = 4 ) prevlastv
  character storev
  real ( kind = 8 ) t(ldt,*)
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) v(ldv,*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  if ( lsame ( direct, 'F' ) ) then

    prevlastv = n

    do i = 1, k

      prevlastv = max ( i, prevlastv )
!
!  H(i)  =  I
!
      if ( tau(i) == 0.0D+00 ) then

        t(1:i,i) = 0.0D+00
!
!  General case
!
      else
!
!  Skip any trailing zeros.
!
        if ( lsame ( storev, 'C' ) ) then

          do lastv = n, i + 1, -1
            if ( v(lastv,i) /= 0.0D+00 ) then
              exit
            end if
          end do

          do j = 1, i - 1
            t(j,i) = - tau(i) * v(i,j)
          end do

          j = min ( lastv, prevlastv )
!
!  T(1:i-1,i) := - tau(i) * V(i+1:j,1:i-1)' * V(i+1:j,i)
!
          call dgemv ( 'Transpose', j-i, i-1, -tau(i), v(i+1,1), ldv, &
            v(i+1,i), 1, one, t(1,i), 1 )
!
!  Skip any trailing zeros.
!
        else

          do lastv = n, i + 1, -1
            if ( v(i,lastv) /= 0.0D+00 ) then
              exit
            end if
          end do

          do j = 1, i - 1
            t(j,i) = -tau(i) * v(j,i)
          end do   

          j = min ( lastv, prevlastv )
!
!  T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)'
!
          call dgemv ( 'No transpose', i-1, j-i, -tau(i), v(1,i+1), ldv, &
            v(i,i+1 ), ldv, one, t(1,i), 1 )

        end if
!
!  T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
        call dtrmv ( 'Upper', 'No transpose', 'Non-unit', i-1, t, ldt, &
          t(1,i), 1 )

        t(i,i) = tau(i)

        if ( 1 < i ) then
          prevlastv = max ( prevlastv, lastv )
        else
          prevlastv = lastv
        end if

      end if

    end do

  else

    prevlastv = 1
!
!  H(i) = I.
!
    do i = k, 1, -1

      if ( tau(i) == 0.0D+00 ) then

        t(i:k,i) = 0.0D+00
!
!  General case.
!
      else

        if ( i < k ) then
!
!  Skip any leading zeros.
!
          if ( lsame ( storev, 'C' ) ) then

            do lastv = 1, i - 1
              if ( v(lastv,i) /= 0.0D+00 ) then
                exit
              end if
            end do

            do j = i + 1, k
              t(j,i) = -tau(i) * v(n-k+i,j)
            end do   

            j = max ( lastv, prevlastv )
!
!  T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)' * V(j:n-k+i,i)
!
            call dgemv ( 'Transpose', n-k+i-j, k-i, -tau(i), v(j,i+1), &
              ldv, v(j,i), 1, one, t(i+1,i), 1 )
!
!  Skip any leading zeros.
!
          else

            do lastv = 1, i - 1
              if ( v(i,lastv) /= 0.0D+00 ) then
                exit
              end if
            end do

            do j = i + 1, k
              t(j,i) = -tau(i) * v(j,n-k+i)
            end do
   
            j = max ( lastv, prevlastv )
!
!  T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)'
!
            call dgemv ( 'No transpose', k-i, n-k+i-j, -tau(i), v(i+1,j), &
              ldv, v(i,j), ldv, one, t(i+1,i), 1 )

          end if
!
!  T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
          call dtrmv ( 'Lower', 'No transpose', 'Non-unit', k-i, &
            t(i+1,i+1), ldt, t(i+1,i), 1 )

          if ( 1 < i ) then
            prevlastv = min ( prevlastv, lastv )
          else
            prevlastv = lastv
          end if

        end if

        t(i,i) = tau(i)

      end if

    end do

  end if

  return
end
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!    This routine uses real ( kind = 8 ) real arithmetic.
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
  implicit none

  real ( kind = 8 ) absxi
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm
  real ( kind = 8 ) scale
  real ( kind = 8 ) ssq
  real ( kind = 8 ) x(*)

  if ( n < 1 .or. incx < 1 ) then

    norm  = 0.0D+00

  else if ( n == 1 ) then

    norm  = abs ( x(1) )

  else

    scale = 0.0D+00
    ssq = 1.0D+00

    do ix = 1, 1 + ( n - 1 )*incx, incx
      if ( x(ix) /= 0.0D+00 ) then
        absxi = abs ( x(ix) )
        if ( scale < absxi ) then
          ssq = 1.0D+00 + ssq * ( scale / absxi )**2
          scale = absxi
        else
          ssq = ssq + ( absxi / scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt ( ssq )
  end if

  dnrm2 = norm

  return
end
subroutine dorm2r ( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )

!*****************************************************************************80
!
!! DORM2R overwrites M by N matrix C with Q*C, Q'*C, C*Q or C*Q'.
!
!  Discussion:
!
!    DORM2R overwrites the general real m by n matrix C with
!      Q * C  if SIDE = 'L' and TRANS = 'N', or
!      Q'* C  if SIDE = 'L' and TRANS = 'T', or
!      C * Q  if SIDE = 'R' and TRANS = 'N', or
!      C * Q' if SIDE = 'R' and TRANS = 'T',
!    where Q is a real orthogonal matrix defined as the product of k
!    elementary reflectors
!      Q = H(1) H(2) . . . H(k)
!    as returned by DGEQRF.  Q is of order m if SIDE = 'L' and of order n
!    if SIDE = 'R'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE:
!    'L': apply Q or Q' from the Left
!    'R': apply Q or Q' from the Right
!
!    Input, character TRANS:
!    'N': apply Q  (No transpose)
!    'T': apply Q' (Transpose)
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix C. 
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix C. 
!    0 <= N.
!
!    Input, integer ( kind = 4 ) K, the number of elementary reflectors 
!    whose product defines the matrix Q.
!    If SIDE = 'L', 0 <= K <= M;
!    if SIDE = 'R', 0 <= K <= N;
!
!    Input, real ( kind = 8 ) A(LDA,K), the matrix.
!    The i-th column must contain the vector which defines the
!    elementary reflector H(i), for i = 1,2,...,k, as returned by
!    DGEQRF in the first k columns of its array argument A.
!    A is modified by the routine but restored on exit.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    If SIDE = 'L', max(1,M) <= LDA;
!    if SIDE = 'R', max(1,N) <= LDA.
!
!    Input, real ( kind = 8 ) TAU(K).
!    TAU(i) must contain the scalar factor of the elementary
!    reflector H(i), as returned by DGEQRF.
!
!    Input/output, real ( kind = 8 ) C(LDC,N).
!    On entry, the m by n matrix C.
!    On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!    Input, integer ( kind = 4 ) LDC, the leading dimension of the array C.
!    max(1,M) <= LDC.
!
!    Input/output, real (kind = 8 ) WORK(*), workspace.
!    Size must be at least:
!    * N if SIDE = 'L', 
!    * M if SIDE = 'R'.
!
!    Output, integer ( kind = 4 ) INFO, return flag:
!    = 0: successful exit
!    < 0: the (-INFO)-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldc

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) aii
  real ( kind = 8 ) c(ldc,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) info
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  logical left
  logical lsame
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mi
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nq
  logical notran
  character side
  real ( kind = 8 ) tau(*)
  character trans
  real ( kind = 8 ) work(*)
!
!  Test the input arguments.
!
  info = 0
  left = lsame ( side, 'L' )
  notran = lsame ( trans, 'N' )
!
!  NQ is the order of Q.
!
  if ( left ) then
    nq = m
  else
    nq = n
  end if

  if ( .not. left .and. .not. lsame ( side, 'R' ) ) then
    info = -1
  else if ( .not. notran .and. .not. lsame ( trans, 'T' ) ) then
    info = -2
  else if ( m < 0 ) then
    info = -3
  else if ( n < 0 ) then
    info = -4
  else if ( k < 0 .or. nq < k ) then
    info = -5
  else if ( lda < max ( 1, nq ) ) then
    info = -7
  else if ( ldc < max ( 1, m ) ) then
    info = -10
  end if

  if ( info /= 0 ) then
    call xerbla ( 'DORM2R', -info )
    return
  end if
!
!  Quick return if possible
!
  if ( m == 0 .or. n == 0 .or. k == 0 ) then
    return
  end if

  if ( (       left .and. .not. notran ) .or. &
       ( .not. left .and.       notran ) ) then
    i1 = 1
    i2 = k
    i3 = 1
  else
    i1 = k
    i2 = 1
    i3 = -1
  end if

  if ( left ) then
    ni = n
    jc = 1
  else
    mi = m
    ic = 1
  end if

  do i = i1, i2, i3
!
!  H(i) is applied to C(i:m,1:n).
!
    if ( left ) then

      mi = m - i + 1
      ic = i
!
!  H(i) is applied to C(1:m,i:n).
!
    else

      ni = n - i + 1
      jc = i

    end if
!
!  Apply H(i).
!
     aii = a(i,i)
     a(i,i) = 1.0D+00
     call dlarf ( side, mi, ni, a(i,i), 1, tau(i), c(ic,jc), ldc, work )
     a(i,i) = aii

  end do

  return
end
subroutine dormqr ( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, &
  info )

!*****************************************************************************80
!
!! DORMQR replaces the rectangular matrix C by Q*C, Q'*C, C*Q or C*Q'.
!
!  Discussion:
!
!    DORMQR overwrites the general real M-by-N matrix C with
!
!                   SIDE = 'L'   SIDE = 'R'
!      TRANS = 'N':   Q  * C        C * Q
!      TRANS = 'T':   Q' * C        C * Q'
!
!    where Q is a real orthogonal matrix defined as the product of k
!    elementary reflectors
!      Q = H(1) H(2) . . . H(k)
!    as returned by DGEQRF.  Q is of order M if SIDE = 'L' and of order N
!    if SIDE = 'R'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE:
!    'L': apply Q or Q' from the Left;
!    'R': apply Q or Q' from the Right.
!
!    Input, character TRANS:
!    'N':  No transpose, apply Q;
!    'T':  Transpose, apply Q'.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix C. 
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix C.
!    0 <= N.
!
!    Input, integer ( kind = 4 ) K, the number of elementary reflectors 
!    whose product defines the matrix Q.
!      If SIDE = 'L', M >= K >= 0;
!      if SIDE = 'R', N >= K >= 0.
!
!    Input, real ( kind = 8 ) A(LDA,K).
!    The i-th column must contain the vector which defines the
!    elementary reflector H(i), for I = 1,2,...,K, as returned by
!    DGEQRF in the first k columns of its array argument A.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    If SIDE = 'L', max(1,M) <= LDA;
!    if SIDE = 'R', max(1,N) <= LDA.
!
!    Input, real ( kind = 8 ) TAU(K), TAU(I) contains the scalar factor 
!    of the elementary reflector H(i), as returned by DGEQRF.
!
!    Input/output, real ( kind = 8 ) C(LDC,N).
!    On entry, the M-by-N matrix C.
!    On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!    Input, integer ( kind = 4 ) LDC, the leading dimension of the array C. 
!    max ( 1, M ) <= LDC.
!
!    Input/output, real ( kind = 8 ) WORK(MAX(1,LWORK)), workspace.
!    On exit, if INFO = 0, WORK(1) returns the optimal value of LWORK.
!
!    Input, integer ( kind = 4 ) LWORK, the dimension of the array WORK.
!    If SIDE = 'L', max(1,N) <= LWORK;
!    if SIDE = 'R', max(1,M) <= LWORK.
!    For optimum performance N*NB <= LWORK if SIDE = 'L', and
!    M*NB <= LWORK if SIDE = 'R', where NB is the optimal blocksize.
!    If LWORK = -1, then a workspace query is assumed; the routine
!    only calculates the optimal size of the work array, returns
!    this value as the first entry of the work array, and no error
!    message related to Lwork is issued by XERBLA.
!
!    Output, integer ( kind = 4 ) INFO, return flag.
!    0:  successful exit.
!    < 0:  if (-INFO)-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbmax = 64

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldc
  integer ( kind = 4 ), parameter :: ldt = nbmax + 1

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) c(ldc,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) iinfo
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iws
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ldwork
  logical left
  logical lquery
  logical lsame
  integer ( kind = 4 ) lwkopt
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mi
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbmin
  integer ( kind = 4 ) ni
  logical notran
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nw
  character side
  real ( kind = 8 ) t(ldt,nbmax)
  real ( kind = 8 ) tau(*)
  character trans
  real ( kind = 8 ) work(*)
!
!  Test the input arguments.
!
  info = 0
  left = lsame ( side, 'L' )
  notran = lsame ( trans, 'N' )
  lquery = ( lwork == -1 )
!
!  NQ is the order of Q and NW is the minimum dimension of WORK.
!
  if ( left ) then
    nq = m
    nw = n
  else
    nq = n
    nw = m
  end if

  if ( .not. left .and. .not. lsame ( side, 'R' ) ) then
    info = -1
  else if ( .not. notran .and. .not. lsame ( trans, 'T' ) ) then
    info = -2
  else if ( m < 0 ) then
    info = -3
  else if ( n < 0 ) then
    info = -4
  else if ( k < 0 .or. nq < k ) then
    info = -5
  else if ( lda < max ( 1, nq ) ) then
    info = -7
  else if ( ldc < max ( 1, m ) ) then
    info = -10
  else if ( lwork < max ( 1, nw ) .and. .not. lquery ) then
    info = -12
  end if
!
!  Determine the block size.  NB may be at most NBMAX, where NBMAX
!  is used to define the local array T.
!
  if ( info == 0 ) then
    nb = min ( nbmax, ilaenv ( 1, 'DORMQR', side // trans, m, n, k, -1 ) )
    lwkopt = max ( 1, nw ) * nb
    work(1) = lwkopt
  end if

  if ( info /= 0 ) then
    call xerbla ( 'DORMQR', -info )
    return
  else if ( lquery ) then
    return
  end if
!
!  Quick return if possible.
!
  if ( m == 0 .or. n == 0 .or. k == 0 ) then
     work(1) = 1
     return
  end if

  nbmin = 2
  ldwork = nw

  if ( 1 < nb .and. nb < k ) then
    iws = nw * nb
    if ( lwork < iws ) then
      nb = lwork / ldwork
      nbmin = max ( 2, ilaenv ( 2, 'DORMQR', side // trans, m, n, k, - 1 ) )
    end if
  else
    iws = nw
  end if
!
!  Use unblocked code.
!
  if ( nb < nbmin .or. k <= nb ) then

    call dorm2r ( side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo )
!
!  Use blocked code.
!
  else

    if ( ( left .and. .not. notran ) .or. &
         ( .not. left .and. notran ) ) then
      i1 = 1
      i2 = k
      i3 = nb
    else
      i1 = ( ( k - 1 ) / nb ) * nb + 1
      i2 = 1
      i3 = -nb
    end if

    if ( left ) then
      ni = n
      jc = 1
    else
      mi = m
      ic = 1
    end if

    do i = i1, i2, i3

      ib = min ( nb, k-i+1 )
!
!  Form the triangular factor of the block reflector
!  H = H(i) H(i+1) . . . H(i+ib-1)
!
      call dlarft ( 'Forward', 'Columnwise', nq-i+1, ib, a(i,i), lda, tau(i), &
        t, ldt )
!
!  H or H' is applied to C(i:m,1:n)
!
      if ( left ) then

        mi = m - i + 1
        ic = i
!
!  H or H' is applied to C(1:m,i:n)
!
      else

        ni = n - i + 1
        jc = i

      end if
!
!  Apply H or H'
!
      call dlarfb ( side, trans, 'Forward', 'Columnwise', mi, ni, ib, &
        a(i,i), lda, t, ldt, c(ic,jc), ldc, work, ldwork )

    end do

  end if

  work(1) = lwkopt

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!
!    This routine uses real ( kind = 8 ) real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m + 1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! DSWAP interchanges two vectors.
!
!  Discussion:
!
!    This routine uses real ( kind = 8 ) real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m + 1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )

!*****************************************************************************80
!
!! DTRMM performs B:=A*B or B:=B*A, A triangular, B rectangular.
!
!  Discussion:
!
!    This routine performs one of the matrix-matrix operations
!      B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!      op( A ) = A   or   op( A ) = A'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE, specifies whether op(A) multiplies B from
!    the left or right as follows:
!    'L' or 'l': B := alpha*op( A )*B.
!    'R' or 'r': B := alpha*B*op( A ).
!
!    Input, character UPLO, specifies whether the matrix A is an upper or
!    lower triangular matrix as follows:
!    'U' or 'u': A is an upper triangular matrix.
!    'L' or 'l': A is a lower triangular matrix.
!
!    Input, character TRANS, specifies the form of op( A ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n': op( A ) = A.
!    'T' or 't': op( A ) = A'.
!    'C' or 'c': op( A ) = A'.
!
!    Input, character DIAG, specifies whether or not A is unit triangular
!    as follows:
!    'U' or 'u': A is assumed to be unit triangular.
!    'N' or 'n': A is not assumed to be unit triangular.
!
!    Input, integer ( kind = 4 ) M, the number of rows of B.  0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of B.  
!    0 <= N.
!
!    Input, real ( kind = 8 ) ALPHA, the scalar  alpha.  When alpha is
!    zero, A is not referenced and B need not be set before entry.
!
!    Input, real ( kind = 8 ) A(LDA,K), where k is m when  SIDE = 'L' or 'l'  
!    and is  n  when  SIDE = 'R' or 'r'.
!    Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!    upper triangular part of the array  A must contain the upper
!    triangular matrix  and the strictly lower triangular part of
!    A is not referenced.
!    Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!    lower triangular part of the array  A must contain the lower
!    triangular matrix  and the strictly upper triangular part of
!    A is not referenced.
!    Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!    A  are not referenced either,  but are assumed to be  unity.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of A as declared
!    in the calling program.  When SIDE = 'L' or 'l' then LDA must be at 
!    least max ( 1, M ); when SIDE = 'R' or 'r', LDA must be at least 
!    max ( 1, N ).
!
!    Input/output, real ( kind = 8 ) B(LDB,N).
!    Before entry, the leading m by n part of the array  B must contain 
!    the matrix  B, and on exit is overwritten by the transformed matrix.
!
!    Input, integer ( kind = 4 ) LDB, the first dimension of B as declared
!    in  the  calling program.   max ( 1, M ) <= LDB.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(ldb,*)
  character diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical lsame
  logical lside
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nounit
  integer ( kind = 4 ) nrowa
  real ( kind = 8 ) one
  parameter ( one = 1.0D+00 )
  character side
  real ( kind = 8 ) temp
  character transa
  character uplo
  logical upper
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Test the input parameters.
!
  lside  = lsame ( side  , 'l' )

  if ( lside ) then
    nrowa = m
  else
    nrowa = n
  end if

  nounit = lsame ( diag  , 'n' )
  upper  = lsame ( uplo  , 'u' )

  info = 0
  if (      ( .not. lside                ) .and. &
            ( .not. lsame ( side  , 'r' ) )      ) then
    info = 1
  else if ( ( .not. upper                ) .and. &
              ( .not. lsame ( uplo  , 'l' ) )      ) then
    info = 2
  else if ( ( .not. lsame ( transa, 'n' ) ) .and. &
              ( .not. lsame ( transa, 't' ) ) .and. &
             ( .not. lsame ( transa, 'c' ) )      ) then
    info = 3
  else if ( ( .not. lsame ( diag  , 'u' ) ) .and. &
              ( .not. lsame ( diag  , 'n' ) )      ) then
    info = 4
  else if ( m < 0 ) then
    info = 5
  else if ( n < 0 ) then
    info = 6
  else if ( lda < max ( 1, nrowa ) ) then
    info = 9
  else if ( ldb < max ( 1, m ) ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtrmm ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if
!
!  And when alpha is zero.
!
  if ( alpha == 0.0D+00 ) then
    b(1:m,1:n) = 0.0D+00
    return
  end if
!
!  Start the operations.
!
  if ( lside ) then
!
!  Form  B := alpha*A*B.
!
    if ( lsame ( transa, 'n' ) ) then

      if ( upper ) then
        do j = 1, n
          do k = 1, m
            if ( b(k,j) /= 0.0D+00 ) then
              temp = alpha * b(k,j)
              do i = 1, k - 1
                b(i,j) = b(i,j) + temp * a(i,k)
              end do
              if ( nounit ) then
                temp = temp * a(k,k)
              end if
              b(k,j) = temp
            end if
          end do
        end do
      else
        do j = 1, n
          do k = m, 1, -1
            if ( b(k,j) /= 0.0D+00 ) then
              temp = alpha * b(k,j)
              b(k,j) = temp
              if ( nounit ) then
                b(k,j) = b(k,j) * a(k,k)
              end if
              do i = k + 1, m
                b(i,j) = b(i,j) + temp * a(i,k)
              end do
            end if
          end do
        end do
      end if
!
!  Form  B := alpha*A'*B.
!
    else

      if ( upper ) then
        do j = 1, n
          do i = m, 1, -1
            temp = b(i,j)
            if ( nounit ) then
              temp = temp * a(i,i)
            end if
            do k = 1, i - 1
              temp = temp + a(k,i) * b(k,j)
            end do
            b(i,j) = alpha * temp
          end do
        end do
      else
        do j = 1, n
          do i = 1, m
            temp = b(i,j)
            if ( nounit ) then
              temp = temp * a(i,i)
            end if
            do k = i + 1, m
              temp = temp + a(k,i) * b(k,j)
            end do
            b(i,j) = alpha * temp
          end do
        end do
      end if
    end if

  else
!
!  Form  B := alpha*B*A.
!
    if ( lsame ( transa, 'n' ) ) then

      if ( upper ) then

        do j = n, 1, -1
          temp = alpha
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = 1, m
            b(i,j) = temp * b(i,j)
          end do
          do k = 1, j - 1
            if ( a(k,j) /= 0.0D+00 ) then
              temp = alpha * a(k,j)
              do i = 1, m
                b(i,j) = b(i,j) + temp * b(i,k)
              end do
            end if
          end do
        end do

      else

        do j = 1, n
          temp = alpha
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = 1, m
            b(i,j) = temp * b(i,j)
          end do
          do k = j + 1, n
            if ( a(k,j) /= 0.0D+00 ) then
              temp = alpha * a(k,j)
              do i = 1, m
                b(i,j) = b(i,j) + temp * b(i,k)
              end do
            end if
          end do
        end do

      end if
!
!  Form  B := alpha*B*A'.
!
    else

      if ( upper ) then
 
        do k = 1, n
          do j = 1, k - 1
            if ( a(j,k) /= 0.0D+00 ) then
              temp = alpha * a(j,k)
              do i = 1, m
                b(i,j) = b(i,j) + temp * b(i,k)
              end do
            end if
          end do
          temp = alpha
          if ( nounit ) then
            temp = temp * a(k,k)
          end if
          if ( temp /= one ) then
            do i = 1, m
              b(i,k) = temp * b(i,k)
            end do
          end if
        end do

      else

        do k = n, 1, -1
          do j = k + 1, n
            if ( a(j,k) /= 0.0D+00 ) then
              temp = alpha * a(j,k)
              do i = 1, m
                b(i,j) = b(i,j) + temp * b(i,k)
              end do
            end if
          end do
          temp = alpha
          if ( nounit ) then
            temp = temp * a(k,k)
          end if
          if ( temp /= one ) then
            do i = 1, m
              b(i,k) = temp * b(i,k)
            end do
          end if
        end do

      end if

    end if

  end if

  return
end
subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )

!*****************************************************************************80
!
!! DTRMV computes x: = A*x or x = A'*x for a triangular matrix A.
!
!  Discussion:
!
!    DTRMV performs one of the matrix-vector operations
!
!      x := A*x,   or   x := A'*x,
!
!    where x is an n element vector and  A is an n by n unit, or non-unit,
!    upper or lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character UPLO, specifies whether the matrix is an upper or
!    lower triangular matrix as follows:
!    'U' or 'U': A is an upper triangular matrix.
!    'L' or 'L': A is a lower triangular matrix.
!
!    Input, character TRANS, specifies the operation to be performed as
!    follows:
!    'N' or 'n': x := A*x.
!    'T' or 't': x := A'*x.
!    'C' or 'c': x := A'*x.
!
!    Input, character DIAG, specifies whether or not A is unit
!    triangular as follows:
!    'U' or 'u': A is assumed to be unit triangular.
!    'N' or 'n': A is not assumed to be unit triangular.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!    0 <= N.
!
!    Input, real ( kind = 8 ) A(LDA,N).
!    Before entry with  UPLO = 'U' or 'U', the leading n by n
!    upper triangular part of the array A must contain the upper
!    triangular matrix and the strictly lower triangular part of
!    A is not referenced.
!    Before entry with UPLO = 'L' or 'L', the leading n by n
!    lower triangular part of the array A must contain the lower
!    triangular matrix and the strictly upper triangular part of
!    A is not referenced.
!    Note that when  DIAG = 'U' or 'U', the diagonal elements of
!    A are not referenced either, but are assumed to be unity.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of A as declared
!    in the calling program. max ( 1, N ) <= LDA.
!
!    Input/output, real ( kind = 8 ) X(1+(N-1)*abs( INCX)).
!    Before entry, the incremented array X must contain the n
!    element vector x. On exit, X is overwritten with the
!    tranformed vector x.
!
!    Input, integer ( kind = 4 ) INCX, the increment for the elements of
!    X.  INCX must not be zero.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  character diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jx
  integer ( kind = 4 ) kx
  integer ( kind = 4 ) n
  logical nounit
  logical lsame
  real ( kind = 8 ) temp
  character trans
  character uplo
  real ( kind = 8 ) x(*)
!
!  Test the input parameters.
!
  info = 0
  if  ( .not. lsame ( uplo , 'U' ) .and.  &
        .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
            .not. lsame ( trans, 'T' ) .and.  &
            .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
            .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( lda < max ( 1, n ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtrmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX  too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form x := A*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0D+00 ) then
            temp = x(j)
            do i = 1, j - 1
              x(i) = x(i) + temp * a(i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(j,j)
            end if
          end if
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0D+00 ) then
            temp = x(jx)
            ix = kx
            do i = 1, j - 1
              x(ix) = x(ix) + temp * a(i,j)
              ix = ix + incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(j,j)
            end if
          end if
          jx = jx + incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0D+00 ) then
            temp = x(j)
            do i = n, j + 1, -1
              x(i) = x(i) + temp * a(i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(j,j)
            end if
          end if
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          if ( x(jx) /= 0.0D+00 ) then
            temp = x(jx)
            ix = kx
            do i = n, j + 1, -1
              x(ix) = x(ix) + temp * a(i,j)
              ix = ix - incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(j,j)
            end if
          end if
          jx = jx - incx
        end do
      end if
    end if
  else
!
!  Form x := A'*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j - 1, 1, -1
            temp = temp + a(i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        jx = kx + ( n - 1 ) * incx
        do j = n, 1, -1
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j - 1, 1, -1
            ix = ix   - incx
            temp = temp + a(i,j) * x(ix)
          end do
          x(jx) = temp
          jx = jx - incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j + 1, n
            temp = temp + a(i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j + 1, n
            ix = ix + incx
            temp = temp + a(i,j) * x(ix)
          end do
          x(jx) = temp
          jx = jx + incx
        end do
      end if
    end if
  end if

  return
end
subroutine dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )

!*****************************************************************************80
!
!! DTRSM performs B:=inv(A)*C or B:=C*inv(A), B and C rectangular, A triangular.
!
!  Discussion:
!
!    DTRSM  solves one of the matrix equations
!      op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!      op( A ) = A   or   op( A ) = A'.
!    The matrix X is overwritten on B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character SIDE, specifies whether op( A ) appears on the left
!    or right of X as follows:
!    'L' or 'l': op( A )*X = alpha*B.
!    'R' or 'r': X*op( A ) = alpha*B.
!
!    Input, character UPLO, specifies whether the matrix A is an upper or
!    lower triangular matrix as follows:
!    'U' or 'u': A is an upper triangular matrix.
!    'L' or 'l': A is a lower triangular matrix.
!
!    Input, character TRANSA, specifies the form of op( A ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n': op( A ) = A.
!    'T' or 't': op( A ) = A'.
!    'C' or 'c': op( A ) = A'.
!
!    Input, character DIAG, specifies whether or not A is unit triangular
!    as follows:
!    'U' or 'u': A is assumed to be unit triangular.
!    'N' or 'n': A is not assumed to be unit triangular.
!
!    Input, integer ( kind = 4 ) M, the number of rows of B.  0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of B.  0 <= N.
!
!    Input, real ( kind = 8 ) ALPHA, the scalar  alpha.  When alpha is
!    zero then A is not referenced and  B need not be set before entry.
!
!    Input, real ( kind = 8 ) A(LDA,K) where K is M when  SIDE = 'L' or 'l'  
!    and K is N when  SIDE = 'R' or 'r'.
!    Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!    upper triangular part of the array  A must contain the upper
!    triangular matrix  and the strictly lower triangular part of
!    A is not referenced.
!    Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!    lower triangular part of the array  A must contain the lower
!    triangular matrix  and the strictly upper triangular part of
!    A is not referenced.
!    Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!    A  are not referenced either,  but are assumed to be  unity.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of A as declared
!    in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!    LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!    then LDA must be at least max( 1, n ).
!
!    Input/output, real ( kind = 8 ) B(LDB,N).
!    Before entry, the leading m by n part of the array B must
!    contain the right-hand side matrix B, and on exit is
!    overwritten by the solution matrix X.
!
!    Input, integer ( kind = 4 ) LDB, the first dimension of B as declared
!    in the calling program.  LDB must be at least max ( 1, M ).
!
  implicit none

  integer lda
  integer ldb

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(ldb,*)
  character diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical lsame
  logical lside
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nounit
  integer ( kind = 4 ) nrowa
  real ( kind = 8 ), parameter :: one = 1.0D+00
  character side
  real ( kind = 8 ) temp
  character transa
  character uplo
  logical upper
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Test the input parameters.
!
  lside  = lsame ( side, 'L' )

  if ( lside ) then
    nrowa = m
  else
    nrowa = n
  end if

  nounit = lsame ( diag, 'N' )
  upper = lsame ( uplo, 'U' )

  info = 0

  if (      ( .not. lside                 ) .and. &
            ( .not. lsame ( side  , 'R' ) )      ) then
    info = 1
  else if ( ( .not. upper                 ) .and. &
            ( .not. lsame ( uplo  , 'L' ) )      ) then
    info = 2
  else if ( ( .not. lsame ( transa, 'n' ) ) .and. &
            ( .not. lsame ( transa, 't' ) ) .and. &
            ( .not. lsame ( transa, 'C' ) )      ) then
    info = 3
  else if ( ( .not. lsame ( diag  , 'u' ) ) .and. &
            ( .not. lsame ( diag  , 'n' ) )      ) then
    info = 4
  else if ( m < 0 ) then
    info = 5
  else if ( n < 0 ) then
    info = 6
  else if ( lda < max ( 1, nrowa ) ) then
    info = 9
  else if ( ldb < max ( 1, m ) ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'DTRSM ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if
!
!  and when alpha is zero.
!
  if ( alpha == 0.0D+00 ) then
    b(1:m,1:n) = 0.0D+00
    return
  end if
!
!  Start the operations.
!
  if ( lside ) then
!
!  Form  B := alpha*inv( a )*B.
!
    if ( lsame ( transa, 'n' ) ) then

      if ( upper ) then
        do j = 1, n
          if ( alpha /= one ) then
            do i = 1, m
              b(i,j) = alpha * b(i,j)
            end do
          end if
          do k = m, 1, -1
            if ( b(k,j) /= 0.0D+00 ) then
              if ( nounit ) then
                b(k,j) = b(k,j) / a(k,k)
              end if
              do i = 1, k - 1
                b(i,j) = b(i,j) - b(k,j) * a(i,k)
              end do
            end if
          end do
        end do
      else
        do j = 1, n
          if ( alpha /= one ) then
            do i = 1, m
              b(i,j) = alpha * b(i,j)
            end do
          end if
          do k = 1, m
            if ( b(k,j) /= 0.0D+00 ) then
              if ( nounit ) then
                b(k,j) = b(k,j) / a(k,k)
              end if
              do i = k + 1, m
                b(i,j) = b(i,j) - b(k,j) * a(i,k)
              end do
            end if
          end do
        end do
      end if
!
!  Form  B := alpha*inv( A' )*B.
!
    else

      if ( upper ) then
        do j = 1, n
          do i = 1, m
            temp = alpha * b(i,j)
            do k = 1, i - 1
              temp = temp - a(k,i) * b(k,j)
            end do
            if ( nounit ) then
              temp = temp / a(i,i)
            end if
            b(i,j) = temp
          end do
        end do
      else
        do j = 1, n
          do i = m, 1, -1
            temp = alpha * b(i,j)
            do k = i + 1, m
              temp = temp - a(k,i) * b(k,j)
            end do
            if ( nounit ) then
              temp = temp / a(i,i)
            end if
            b(i,j) = temp
          end do
        end do
      end if
    end if
!
!  Form  B := alpha*B*inv( A ).
!
  else

    if ( lsame ( transa, 'n' ) ) then

      if ( upper ) then

        do j = 1, n
          if ( alpha /= one ) then
            do i = 1, m
              b(i,j) = alpha * b(i,j)
            end do
          end if
          do k = 1, j - 1
            if ( a(k,j) /= 0.0D+00 ) then
              do i = 1, m
                b(i,j) = b(i,j) - a(k,j) * b(i,k)
              end do
            end if
          end do
          if ( nounit ) then
            temp = one / a(j,j)
            do i = 1, m
              b(i,j) = temp * b(i,j)
            end do
          end if
        end do

      else

        do j = n, 1, -1
          if ( alpha /= one ) then
            do i = 1, m
              b(i,j) = alpha * b(i,j)
            end do
          end if
          do k = j + 1, n
            if ( a(k,j) /= 0.0D+00 ) then
              do i = 1, m
                b(i,j) = b(i,j) - a(k,j) * b(i,k)
              end do
            end if
          end do
          if ( nounit ) then
            temp = one / a(j,j)
            do i = 1, m
              b(i,j) = temp * b(i,j)
            end do
          end if
        end do
      end if
!
!  Form  B := alpha*B*inv( A' ).
!
    else

      if ( upper ) then
        do k = n, 1, -1
          if ( nounit ) then
            temp = one / a(k,k)
            do i = 1, m
              b(i,k) = temp * b(i,k)
            end do
          end if
          do j = 1, k - 1
            if ( a(j,k) /= 0.0D+00 ) then
              temp = a(j,k)
              do i = 1, m
                b(i,j) = b(i,j) - temp * b(i,k)
              end do
            end if
          end do
          if ( alpha /= one ) then
            do i = 1, m
              b(i,k) = alpha * b(i,k)
            end do
          end if
        end do
      else
        do k = 1, n
          if ( nounit ) then
            temp = one / a(k,k)
            do i = 1, m
              b(i,k) = temp * b(i,k)
            end do
          end if
          do j = k + 1, n
            if ( a(j,k) /= 0.0D+00 ) then
              temp = a(j,k)
              do i = 1, m
                b(i,j) = b(i,j) - temp * b(i,k)
              end do
            end if
          end do
          if ( alpha /= one ) then
            do i = 1, m
              b(i,k) = alpha * b(i,k)
            end do
          end if
        end do
      end if
    end if
  end if

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX indexes the array element of maximum absolute value.
!
!  Discussion:
!
!    This routine uses real ( kind = 8 ) real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of SX.
!
!    Output, integer ( kind = 4 ) IDAMAX, the index of the element of SX of 
!    maximum absolute value.
!
  implicit none

  real ( kind = 8 ) dmax
  real ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n

  idamax = 0

  if ( n < 1 .or. incx <= 0 ) then
    return
  end if

  idamax = 1

  if ( n == 1 ) then
    return
  end if

  if ( incx == 1 ) then

    dmax = abs ( dx(1) )

    do i = 2, n
      if ( dmax < abs ( dx(i) ) ) then
        idamax = i
        dmax = abs ( dx(i) )
      end if
    end do

  else

    ix = 1
    dmax = abs ( dx(1) )
    ix = ix + incx

    do i = 2, n
      if ( dmax < abs ( dx(ix) ) ) then
        idamax = i
        dmax = abs ( dx(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function ieeeck ( ispec, zero, one )

!*****************************************************************************80
!
!! IEEECK verifies that infinity and NAN arithmetic are safe.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISPEC, specifies whether to test just for 
!    infinity arithmetic or also NAN arithmetic.
!    0: Verify infinity arithmetic only.
!    1: Verify infinity and NAN arithmetic.
!
!    Input, real ( kind = 8 ) ZERO, must contain the value 0.0.
!    This is passed to prevent the compiler from optimizing away this code.
!
!    Input, real ( kind = 8 ) ONE, must contain the value 1.0.
!    This is passed to prevent the compiler from optimizing away this code.
!
!    Output, integer ( kind = 4 ) IEEECK, the result.
!    0:  Arithmetic failed to produce the correct answers
!    1:  Arithmetic produced the correct answers
!
  implicit none

  integer ( kind = 4 ) ieeeck
  integer ( kind = 4 ) ispec
  real ( kind = 8 ) nan1
  real ( kind = 8 ) nan2
  real ( kind = 8 ) nan3
  real ( kind = 8 ) nan4
  real ( kind = 8 ) nan5
  real ( kind = 8 ) nan6
  real ( kind = 8 ) neginf
  real ( kind = 8 ) negzro
  real ( kind = 8 ) newzro
  real ( kind = 8 ) one
  real ( kind = 8 ) posinf
  real ( kind = 8 ) zero

  ieeeck = 1

  posinf = one / zero

  if ( posinf <= one ) then
    ieeeck = 0
    return
  end if

  neginf = - one / zero

  if ( zero <= neginf ) then
    ieeeck = 0
    return
  end if

  negzro = one / ( neginf + one )

  if ( negzro /= zero ) then
    ieeeck = 0
    return
  end if

  neginf = one / negzro

  if ( zero <= neginf ) then
    ieeeck = 0
    return
  end if

  newzro = negzro + zero

  if ( newzro /= zero ) then
    ieeeck = 0
    return
  end if

  posinf = one / newzro

  if ( posinf <= one ) then
    ieeeck = 0
    return
  end if

  neginf = neginf * posinf

  if ( zero <= neginf ) then
    ieeeck = 0
    return
  end if

  posinf = posinf * posinf

  if ( posinf <= one ) then
    ieeeck = 0
    return
  end if
!
!  Return if we were only asked to check infinity arithmetic.
!
  if ( ispec == 0 ) then
    return
  end if

  nan1 = posinf + neginf

  nan2 = posinf / neginf

  nan3 = posinf / posinf

  nan4 = posinf * zero

  nan5 = neginf * negzro

  nan6 = nan5 * zero

  if ( nan1 == nan1 ) then
    ieeeck = 0
    return
  end if

  if ( nan2 == nan2 ) then
    ieeeck = 0
    return
  end if

  if ( nan3 == nan3 ) then
    ieeeck = 0
    return
  end if

  if ( nan4 == nan4 ) then
    ieeeck = 0
    return
  end if

  if ( nan5 == nan5 ) then
    ieeeck = 0
    return
  end if

  if ( nan6 == nan6 ) then
    ieeeck = 0
    return
  end if

  return
end
function iladlc ( m, n, a, lda )

!*****************************************************************************80
!
!! ILADLC scans matrix A for its last nonzero column.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!
!    Input  real ( kind = 8 ) A(LDA,N), the M by N array.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, m ) <= LDA.
!
!    Output, integer ( kind = 4 ) ILADLC, the index of the last nonzero 
!    column of A.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iladlc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  iladlc = -1
!
!  Quick test for the common case where one corner is non-zero.
!
  if ( n == 0 ) then
    iladlc = n
  else if ( a(1,n) /= 0.0D+00 .or. a(m,n) /= 0.0D+00 ) then
    iladlc = n
!
!  Now scan each column from the end, returning with the first non-zero.
!
  else
    do iladlc = n, 1, -1
      do i = 1, m
        if ( a(i,iladlc) /= 0.0D+00 ) then
          return
        end if
      end do
    end do
  end if

  return
end
function iladlr ( m, n, a, lda )

!*****************************************************************************80
!
!! ILADLR returns the index of the last nonzero row of a matrix.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!
!    Input  real ( kind = 8 ) A(LDA,N), the M by N array.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, m ) <= LDA.
!
!    Output, integer ( kind = 4 ) ILADLR, the index of the last nonzero 
!    row of A.
!
  implicit none

  integer ( kind = 4 ) lda

  real ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iladlr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  iladlr = -1
!
!  Quick test for the common case where one corner is non-zero.
!
  if ( m == 0 ) then
    iladlr = m
  else if ( a(m,1) /= 0.0D+00 .or. a(m,n) /= 0.0D+00 ) then
    iladlr = m
!
!  Scan up each column tracking the last zero row seen.
!
  else

    iladlr = 0
    do j = 1, n
      i = m
      do while ( a(max(i,1),j) == 0.0D+00 .and. 1 <= i )
        i = i - 1
      end do
      iladlr = max ( iladlr, i )
    end do

  end if

  return
end
function ilaenv ( ispec, name, opts, n1, n2, n3, n4 ) 

!*****************************************************************************80
!
!! ILAENV returns problem-dependent parameters.
!
!  Discussion:
!
!    ILAENV is called from the LAPACK routines to choose problem-dependent
!    parameters for the local environment.  See ISPEC for a description of
!    the parameters.
!
!    ILAENV returns an integer.  If the value returned is greater than or
!    equal to zero, it has returned the value of the parameter specified
!    by ISPEC.  Otherwise, the (-ILAENV)-th input argument was illegal.
!
!    This version provides a set of parameters which should give good,
!    but not optimal, performance on many of the currently available
!    computers.  Users are encouraged to modify this subroutine to set
!    the tuning parameters for their particular machine using the option
!    and problem size information in the arguments.
!
!    This routine will not function correctly if it is converted to all
!    lower case.  Converting it to all upper case is allowed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2014
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISPEC, specifies the parameter to be 
!    returned.
!    1: the optimal blocksize; if this value is 1, an unblocked
!       algorithm will give the best performance.
!    2: the minimum block size for which the block routine
!       should be used; if the usable block size is less than
!       this value, an unblocked routine should be used.
!    3: the crossover point (in a block routine, for N less
!       than this value, an unblocked routine should be used)
!    4: the number of shifts, used in the nonsymmetric
!       eigenvalue routines (DEPRECATED)
!    5: the minimum column dimension for blocking to be used;
!       rectangular blocks must have dimension at least k by m,
!       where k is given by ilaenv(2,...) and m by ilaenv(5,...)
!    6: the crossover point for the SVD (when reducing an m by n
!       matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!       this value, a QR factorization is used first to reduce
!       the matrix to a triangular form.)
!    7: the number of processors
!    8: the crossover point for the multishift QR method
!       for nonsymmetric eigenvalue problems (DEPRECATED)
!    9: maximum size of the subproblems at the bottom of the
!       computation tree in the divide-and-conquer algorithm
!       (used by xGELSD and xGESDD)
!    10: ieee NaN arithmetic can be trusted not to trap
!    11: infinity arithmetic can be trusted not to trap
!    12 <= ISPEC <= 16: xHSEQR or one of its subroutines,
!        see IPARMQ for detailed explanation
!
!    Input, character ( len = * ) NAME, the name of the calling subroutine, 
!    in either upper case or lower case.
!
!    Input, character ( len = * ) OPTS, the options to the 
!    subroutine NAME, concatenated into a single character string.  For 
!    example, UPLO = 'U', TRANS = 'T', and DIAG = 'N' for a triangular routine 
!    would be specified as OPTS = 'UTN'.
!
!    Input, integer ( kind = 4 ) N1, N2, N3, N4, problem dimensions for
!    subroutine NAME.  Not all of these may be required.
!
!    Output, integer ( kind = 4 ) ILAENV, the value of the requested
!    parameter.
!
  character c1
  character ( len = 2 ) c2
  character ( len = 3 ) c3
  character ( len = 2 ) c4
  logical cname
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ieeeck
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 ) iparmq
  integer ( kind = 4 ) ispec
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  character ( len = * ) name
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbmin
  integer ( kind = 4 ) nx
  real ( kind = 8 ), parameter :: one = 1.0D+00
  character ( len = * ) opts
  logical sname
  character  ( len = 6 ) subnam
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  if ( 1 <= ispec .and. ispec <= 3 ) then

    ilaenv = 1
    subnam = name
    ic = ichar ( subnam(1:1) )
    iz = ichar ( 'Z' )
!
!  ASCII character set
!
    if ( iz == 90 .or. iz == 122 ) then

      if ( 97 <= ic .and. ic <= 122 ) then
        subnam(1:1) = char ( ic - 32 )
        do i = 2, 6
          ic = ichar ( subnam(i:i) )
          if ( 97 <= ic .and. ic <= 122 ) then
            subnam(i:i) = char ( ic - 32 )
          end if
        end do
      end if
!
!  EBCDIC character set
!
    else if ( iz == 233 .or. iz == 169 ) then

      if ( ( 129 <= ic .and. ic <= 137 ) .or. &
           ( 145 <= ic .and. ic <= 153 ) .or. &
           ( 162 <= ic .and. ic <= 169 ) ) then
        subnam(1:1) = char ( ic + 64 )
        do i = 2, 6
          ic = ichar ( subnam(i:i) )
          if ( ( 129 <= ic .and. ic <= 137 ) .or. &
               ( 145 <= ic .and. ic <= 153 ) .or. &
               ( 162 <= ic .and. ic <= 169 ) ) then
            subnam(i:i) = char ( ic + 64 )
          end if
        end do
      end if
!
!  Prime machines:  ASCII+128
!
    else if ( iz == 218 .or. iz == 250 ) then

      if ( 225 <= ic .and. ic <= 250 ) then
        subnam(1:1) = char ( ic-32 )
        do i = 2, 6
          ic = ichar ( subnam(i:i) )
          if ( 225 <= ic .and. ic <= 250 ) then
            subnam(i:i) = char ( ic - 32 )
          end if
        end do
      end if

    end if

    c1 = subnam(1:1)
    sname = c1 == 'S' .or. c1 == 'D'
    cname = c1 == 'C' .or. c1 == 'Z'

    if ( .not. ( cname .or. sname ) ) then
      return
    end if

    c2 = subnam(2:3)
    c3 = subnam(4:6)
    c4 = c3(2:3)

  end if
!
!  Invalid value for ISPEC
!
  if ( ispec <= 0 ) then

    ilaenv = -1
!
!  ISPEC = 1:  block size
!
!  In these examples, separate code is provided for setting nb for
!  real and complex.  We assume that nb will take the same value in
!  single or double precision.
!
  else if ( ispec == 1 ) then

    nb = 1

    if ( c2 == 'GE' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      else if ( c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. &
                c3 == 'QLF' ) then
        if ( sname ) then
          nb = 32
        else
          nb = 32
        end if
      else if ( c3 == 'HRD' ) then
        if ( sname ) then
          nb = 32
        else
          nb = 32
        end if
      else if ( c3 == 'BRD' ) then
        if ( sname ) then
          nb = 32
        else
          nb = 32
        end if
      else if ( c3 == 'TRI' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      end if

    else if ( c2 == 'PO' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      end if

    else if ( c2 == 'SY' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      else if ( sname .and. c3 == 'TRD' ) then
        nb = 32
      else if ( sname .and. c3 == 'GST' ) then
        nb = 64
      end if

    else if ( cname .and. c2 == 'HE' ) then

      if ( c3 == 'TRF' ) then
        nb = 64
      else if ( c3 == 'TRD' ) then
        nb = 32
      else if ( c3 == 'GST' ) then
        nb = 64
      end if

    else if ( sname .and. c2 == 'OR' ) then

      if ( c3(1:1) == 'G' ) then

        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nb = 32
        end if

      else if ( c3(1:1) == 'M' ) then

        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nb = 32
        end if

      end if

    else if ( cname .and. c2 == 'UN' ) then

      if ( c3(1:1) == 'G' ) then

        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. & 
             c4 == 'BR' ) then
          nb = 32
        end if

      else if ( c3(1:1) == 'M' ) then

        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nb = 32
        end if

      end if

    else if ( c2 == 'GB' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          if ( n4 <= 64 ) then
            nb = 1
          else
            nb = 32
          end if
        else
          if ( n4 <= 64 ) then
            nb = 1
          else
            nb = 32
          end if
        end if
      end if

    else if ( c2 == 'PB' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          if ( n2 <= 64 ) then
            nb = 1
          else
            nb = 32
          end if
        else
          if ( n2 <= 64 ) then
            nb = 1
          else
            nb = 32
          end if
        end if
      end if

    else if ( c2 == 'TR' ) then

      if ( c3 == 'TRI' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      end if

    else if ( c2 == 'LA' ) then

      if ( c3 == 'UUM' ) then
        if ( sname ) then
          nb = 64
        else
          nb = 64
        end if
      end if

    else if ( sname .and. c2 == 'ST' ) then

      if ( c3 == 'EBZ' ) then
        nb = 1
      end if

    end if

    ilaenv = nb
!
!  ISPEC = 2:  minimum block size
!
  else if ( ispec == 2 ) then

    nbmin = 2

    if ( c2 == 'GE' ) then

      if ( c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. &
           c3 == 'QLF' ) then
        if ( sname ) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if ( c3 == 'HRD' ) then
        if ( sname ) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if ( c3 == 'BRD' ) then
        if ( sname ) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if ( c3 == 'TRI' ) then
        if ( sname ) then
          nbmin = 2
        else
          nbmin = 2
        end if
      end if

    else if ( c2 == 'SY' ) then

      if ( c3 == 'TRF' ) then
        if ( sname ) then
          nbmin = 8
        else
          nbmin = 8
        end if
      else if ( sname .and. c3 == 'TRD' ) then
        nbmin = 2
      end if

    else if ( cname .and. c2 == 'HE' ) then

      if ( c3 == 'TRD' ) then
        nbmin = 2
      end if

    else if ( sname .and. c2 == 'OR' ) then

      if ( c3(1:1) == 'G' ) then
        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nbmin = 2
        end if
      else if ( c3(1:1) == 'M' ) then
        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nbmin = 2
        end if
      end if

    else if ( cname .and. c2 == 'UN' ) then

      if ( c3(1:1) == 'G' ) then
        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nbmin = 2
        end if
      else if ( c3(1:1) == 'M' ) then
        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nbmin = 2
        end if
      end if

    end if

    ilaenv = nbmin
!
!  ISPEC = 3:  crossover point
!
  else if ( ispec == 3 ) then

    nx = 0

    if ( c2 == 'GE' ) then

      if ( c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. c3 == 'QLF' ) then
        if ( sname ) then
          nx = 128
        else
          nx = 128
        end if
      else if ( c3 == 'HRD' ) then
        if ( sname ) then
          nx = 128
        else
          nx = 128
        end if
      else if ( c3 == 'BRD' ) then
        if ( sname ) then
          nx = 128
        else
          nx = 128
        end if
      end if

    else if ( c2 == 'SY' ) then

      if ( sname .and. c3 == 'TRD' ) then
        nx = 32
      end if

    else if ( cname .and. c2 == 'HE' ) then

      if ( c3 == 'TRD' ) then
        nx = 32
      end if

    else if ( sname .and. c2 == 'OR' ) then

      if ( c3(1:1) == 'G' ) then

        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nx = 128
        end if

      end if

    else if ( cname .and. c2 == 'UN' ) then

      if ( c3(1:1) == 'G' ) then
        if ( c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. &
             c4 == 'QL' .or. c4 == 'HR' .or. c4 == 'TR' .or. &
             c4 == 'BR' ) then
          nx = 128
        end if
      end if

    end if

    ilaenv = nx
!
!  ISPEC = 4:  number of shifts (used by xHSEQR)
!
  else if ( ispec == 4 ) then

    ilaenv = 6
!
!  ISPEC = 5:  minimum column dimension (not used).
!
  else if ( ispec == 5 ) then

    ilaenv = 2
!
!  ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD).
!
  else if ( ispec == 6 ) then

    ilaenv = int ( real ( min ( n1, n2 ), kind = 4 ) * 1.6E+00 )
!
!  ISPEC = 7:  number of processors (not used).
!
  else if ( ispec == 7 ) then

    ilaenv = 1
!
!  ISPEC = 8:  crossover point for multishift (used by xHSEQR).
!
  else if ( ispec == 8 ) then

    ilaenv = 50
!
!  ISPEC = 9:  maximum size of the subproblems at the bottom of the
!  computation tree in the divide-and-conquer algorithm
!  (used by xGELSD and xGESDD).
!
  else if ( ispec == 9 ) then

    ilaenv = 25
!
!  ISPEC = 10: ieee NaN arithmetic can be trusted not to trap.
!
  else if ( ispec == 10 ) then

    ilaenv = 1
    if ( ilaenv == 1 ) then
      ilaenv = ieeeck ( 1, zero, one )
    end if
!
!  ISPEC = 11: infinity arithmetic can be trusted not to trap.
!
  else if ( ispec == 11 ) then

    ilaenv = 1
    if ( ilaenv == 1 ) then
      ilaenv = ieeeck ( 0, zero, one )
    end if
!
!  12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
!
  else if ( 12 <= ispec .and. ispec <= 16 ) then

    ilaenv = iparmq ( ispec, name, opts, n1, n2, n3, n4 )
!
!  Illegal value.
!
  else

    ilaenv = -1

  end if

  return
end
function iparmq ( ispec, name, opts, n, ilo, ihi, lwork )

!*****************************************************************************80
!
!! IPARMQ sets problem and machine dependent parameters.
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
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISPEC, specifies which tunable parameter 
!    IPARMQ should return.
!    ISPEC=12: (INMIN)  Matrices of order nmin or less
!    are sent directly to xLAHQR, the implicit
!    double shift QR algorithm.  NMIN must be
!    at least 11.
!    ISPEC=13: (INWIN)  Size of the deflation window.
!    This is best set greater than or equal to
!    the number of simultaneous shifts NS.
!    Larger matrices benefit from larger deflation
!    windows.
!    ISPEC=14: (INIBL) Determines when to stop nibbling and
!    invest in an (expensive) multi-shift QR sweep.
!    If the aggressive early deflation subroutine
!    finds LD converged eigenvalues from an order
!    NW deflation window and LD.GT.(NW*NIBBLE)/100,
!    then the next QR sweep is skipped and early
!    deflation is applied immediately to the
!    remaining active diagonal block.  Setting
!    IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!    multi-shift QR sweep whenever early deflation
!    finds a converged eigenvalue.  Setting
!    IPARMQ(ISPEC=14) greater than or equal to 100
!    prevents TTQRE from skipping a multi-shift
!    QR sweep.
!    ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!    a multi-shift QR iteration.
!    ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!    following meanings.
!    0:  During the multi-shift QR sweep,
!    xLAQR5 does not accumulate reflections and
!    does not use matrix-matrix multiply to
!    update the far-from-diagonal matrix
!    entries.
!    1:  During the multi-shift QR sweep,
!    xLAQR5 and/or xLAQRaccumulates reflections and uses
!    matrix-matrix multiply to update the
!    far-from-diagonal matrix entries.
!    2:  During the multi-shift QR sweep.
!    xLAQR5 accumulates reflections and takes
!    advantage of 2-by-2 block structure during
!    matrix-matrix multiplies.
!    (If xTRMM is slower than xGEMM, then
!    IPARMQ(ISPEC=16)=1 may be more efficient than
!    IPARMQ(ISPEC=16)=2 despite the greater level of
!    arithmetic work implied by the latter choice.)
!    
!    Input, character ( len = * ) NAME, the name of the calling subroutine.
!    
!    Input, character ( len = * ) OPTS, a concatenation of the string 
!    arguments to TTQRE.
!    
!    Input, integer ( kind = 4 ) N, the order of the Hessenberg matrix H.
!    
!    Input, integer ( kind = 4 ) ILO, IHI.   It is assumed that H is already
!    upper triangular in rows and columns 1:ILO-1 and IHI+1:N.
!    
!    Input, integer ( kind = 4 ) LWORK, the amount of workspace available.
!
!    Output, integer ( kind = 4 ) IPARMQ, the value of the requested parameter.
!
  implicit none

  integer ( kind = 4 ), parameter :: iacc22 = 16
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), parameter :: inibl = 14
  integer ( kind = 4 ), parameter :: inmin = 12
  integer ( kind = 4 ), parameter :: inwin = 13
  integer ( kind = 4 ) iparmq
  integer ( kind = 4 ), parameter :: ishfts = 15
  integer ( kind = 4 ) ispec
  integer ( kind = 4 ), parameter :: k22min = 14
  integer ( kind = 4 ), parameter :: kacmin = 14
  integer ( kind = 4 ), parameter :: knwswp = 500
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) n
  character ( len = * ) name
  integer ( kind = 4 ) nh
  integer ( kind = 4 ), parameter :: nibble = 14
  integer ( kind = 4 ), parameter :: nmin = 75
  integer ( kind = 4 ) ns
  character ( len = * ) opts
  real ( kind = 4 ), parameter :: two = 2.0E+00
!
!  Set the number of simultaneous shifts.
!
  if ( ( ispec == ishfts ) .or. &
       ( ispec == inwin ) .or. &
       ( ispec == iacc22 ) ) then

    nh = ihi - ilo + 1
    ns = 2

    if ( 30 <= nh ) then
      ns = 4
    end if

    if ( 60 <= nh ) then
      ns = 10
    end if

    if ( 150 <= nh ) then
      ns = max ( 10, nh / nint ( log ( real ( nh ) ) / log ( two ) ) )
    end if

    if ( 590 <= nh ) then
      ns = 64
    end if

    if ( 3000 <= nh ) then
      ns = 128
    end if

    if ( 6000 <= nh ) then
      ns = 256
    end if

    ns = max ( 2, ns - mod ( ns, 2 ) )

  end if
!
!  Matrices of order smaller than NMIN get sent
!  to xLAHQR, the classic double shift algorithm.
!  This must be at least 11.
!
  if ( ispec == inmin ) then

    iparmq = nmin
!
!  INIBL: skip a multi-shift qr iteration and
!  whenever aggressive early deflation finds
!  at least (NIBBLE*(window size)/100) deflations.
!
  else if ( ispec == inibl ) then

    iparmq = nibble
!
!  NSHFTS: The number of simultaneous shifts.
!
  else if ( ispec == ishfts ) then

    iparmq = ns
!
!  NW: deflation window size.
!
  else if ( ispec == inwin ) then

    if ( nh <= knwswp ) then
      iparmq = ns
    else
      iparmq = ( 3 * ns ) / 2
    end if
!
!  IACC22: Whether to accumulate reflections
!  before updating the far-from-diagonal elements
!  and whether to use 2-by-2 block structure while
!  doing it.  A small amount of work could be saved
!  by making this choice dependent also upon the
!  NH=IHI-ILO+1.
!
  else if ( ispec == iacc22 ) then

    iparmq = 0
    if ( kacmin <= ns ) then
      iparmq = 1
    end if

    if ( k22min <= ns ) then
      iparmq = 2
    end if
!
!  Invalid value of ISPEC.
!
  else
    iparmq = -1
  end if

  return
end
function lsame ( ca, cb )

!*****************************************************************************80
!
!! LSAME returns TRUE if CA is the same letter as CB regardless of case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, character CA, CB, the character to compare.
!
!    Output, logical LSAME, is TRUE if the characters are equal,
!    disregarding case.
!
  implicit none

  character ca
  character cb
  integer ( kind = 4 ) inta
  integer ( kind = 4 ) intb
  logical lsame
  integer ( kind = 4 ) zcode
!
!  Test if the characters are equal
!
  lsame = ( ca == cb )

  if ( lsame ) then
    return
  end if
!
!  Now test for equivalence if both characters are alphabetic.
!
  zcode = ichar ( 'Z' )
!
!  Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!  machines, on which ICHAR returns a value with bit 8 set.
!  ICHAR('A') on Prime machines returns 193 which is the same as
!  ICHAR('A') on an EBCDIC machine.
!
  inta = ichar ( ca )
  intb = ichar ( cb )

  if ( zcode == 90 .or. zcode == 122 ) then
!
!  ASCII is assumed - zcode is the ASCII code of either lower or
!  upper case 'Z'.
!
    if ( 97 <= inta .and. inta <= 122 ) then
      inta = inta - 32
    end if

    if ( 97 <= intb .and. intb <= 122 ) then
      intb = intb - 32
    end if

  else if ( zcode == 233 .or. zcode == 169 ) then
!
!  EBCDIC is assumed - zcode is the EBCDIC code of either lower or
!  upper case 'Z'.
!
    if ( 129 <= inta .and. inta <= 137 .or. &
         145 <= inta .and. inta <= 153 .or. &
         162 <= inta .and. inta <= 169 ) then
      inta = inta + 64
    end if

    if ( 129 <= intb .and. intb <= 137 .or. &
         145 <= intb .and. intb <= 153 .or. &
         162 <= intb .and. intb <= 169 ) then
      intb = intb + 64
    end if

  else if ( zcode == 218 .or. zcode == 250 ) then
!
!  ASCII is assumed, on Prime machines - zcode is the ASCII code
!  plus 128 of either lower or upper case 'Z'.
!
    if ( 225 <= inta .and. inta <= 250 ) then
      inta = inta - 32
    end if

    if ( 225 <= intb .and. intb <= 250 ) then
      intb = intb - 32
    end if

  end if

  lsame = ( inta == intb )

  return
end
subroutine xerbla ( srname, info )

!*****************************************************************************80
!
!! XERBLA is an error handler for the LAPACK routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage, 
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, character ( len = * ) SRNAME, the name of the routine
!    which called XERBLA.
!
!    Input, integer ( kind = 4 ) INFO, the position of the invalid parameter in
!    the parameter list of the calling routine.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a,a,a,i2,a)' ) ' ** On entry to ', srname, &
    ' parameter number ', info, ' had an illegal value.'

  stop
end
