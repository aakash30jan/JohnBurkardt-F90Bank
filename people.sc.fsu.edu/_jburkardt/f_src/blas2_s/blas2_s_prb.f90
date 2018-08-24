program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS2_S_PRB.
!
!  Discussion:
!
!    BLAS2_S_PRB tests the BLAS library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS2_S_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS2_S_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SGEMV.
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
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), allocatable :: a(:,:)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character trans
  real ( kind = 4 ), allocatable :: x(:)
  real ( kind = 4 ), allocatable :: y(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a general matrix A,'
  write ( *, '(a)' ) '  SGEMV computes y := alpha * A * x + beta * y'
  write ( *, '(a)' ) '  or             y := alpha * A'' * x + beta * y.'
!
!  y = alpha * A * x + beta * y
!
  trans = 'N'
  m = 5
  n = 4
  alpha = 2.0E+00
  lda = m
  allocate ( a(1:m,1:n) )
  call r4mat_test ( trans, lda, m, n, a )
  allocate ( x(1:n) )
  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do
  incx = 1
  beta = 3.0E+00
  allocate ( y(1:m) )
  do i = 1, m
    y(i) = real ( 10 * i, kind = 4 )
  end do
  incy = 1

  call r4mat_print ( m, n, a, '  Matrix A:' )
  call r4vec_print ( n, x, '  Vector X:' )
  call r4vec_print ( m, y, '  Vector Y:' )

  call sgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

  call r4vec_print ( m, y, '  Result Y = alpha * A  * x + beta * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  y = alpha * A' * x + beta * y
!
  trans = 'T'
  m = 5
  n = 4
  alpha = 2.0E+00
  lda = m
  allocate ( a(1:m,1:n) )
  call r4mat_test ( trans, lda, n, m, a )
  allocate ( x(1:m) )
  do i = 1, m
    x(i) = real ( i, kind = 4 )
  end do
  incx = 1
  beta = 3.0E+00
  allocate ( y(1:n) )
  do i = 1, n
    y(i) = real ( 10 * i, kind = 4 )
  end do
  incy = 1

  call r4mat_print ( m, n, a, '  Matrix A:' )
  call r4vec_print ( m, x, '  Vector X:' )
  call r4vec_print ( n, y, '  Vector Y:' )

  call sgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

  call r4vec_print ( n, y, '  Result Y = alpha * A'' * x + beta * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SGBMV.
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
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: kl = 1
  integer ( kind = 4 ), parameter :: ku = 1
  integer ( kind = 4 ), parameter :: lda = kl + 1 + ku

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character trans
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a general band matrix A,'
  write ( *, '(a)' ) '  SGBMV computes y := alpha * A * x + beta * y'

  trans = 'N'
  alpha = 2.0E+00
  incx = 1
  beta = 3.0E+00
  incy = 1

  do i = 1, m

    jlo = max ( 1, i - kl )
    jhi = min ( n, i + ku )

    do j = jlo, jhi

      if ( i == j ) then
        a(ku+1+i-j,j) = 2.0E+00
      else if ( i == j - 1 .or. i == j + 1 ) then
        a(ku+1+i-j,j) = -1.0E+00
      else
        a(ku+1+i-j,j) = 0.0E+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 4 )
  end do

  call sgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy )

  call r4vec_print ( m, y, '  Result vector y' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SSYMV.
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
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  character uplo
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a general symmetric matrix A,'
  write ( *, '(a)' ) '  SSYMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0E+00
  incx = 1
  beta = 3.0E+00
  incy = 1

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0E+00
      else if ( i == j - 1 ) then
        a(i,j) = -1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( 10 * i, kind = 4 )
  end do

  call ssymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

  call r4vec_print ( n, y, '  Result vector y:' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SSBMV.
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
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: k = 1
  integer ( kind = 4 ), parameter :: lda = k + 1

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character uplo
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a symmetric band matrix A,'
  write ( *, '(a)' ) '  SSBMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0E+00
  incx = 1
  beta = 3.0E+00
  incy = 1

  do i = 1, m

    jhi = min ( n, i + k )

    do j = i, jhi

      if ( i == j ) then
        a(k+1+i-j,j) = 2.0E+00
      else if ( i == j - 1 ) then
        a(k+1+i-j,j) = -1.0E+00
      else
        a(k+1+i-j,j) = 0.0E+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 4 )
  end do

  call ssbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

  call r4vec_print ( m, y, '  Result vector y:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SGER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ), allocatable :: a(:,:)
  real ( kind = 4 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character trans
  real ( kind = 4 ), allocatable :: x(:)
  real ( kind = 4 ), allocatable :: y(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a general matrix A,'
  write ( *, '(a)' ) '  SGER computes A := A + alpha * x * y'''

  m = 5
  n = 4
  alpha = 2.0E+00
  trans = 'N'
  lda = m
  allocate ( a(1:m,1:n) )
  call r4mat_test ( trans, lda, m, n, a )

  allocate ( x(1:m) )
  do i = 1, m
    x(i) = real ( i, kind = 4 )
  end do
  incx = 1

  allocate ( y(1:n) )
  do i = 1, n
    y(i) = real ( 10 * i, kind = 4 )
  end do
  incy = 1

  call r4mat_print ( m, n, a, '  Matrix A:' )
  call r4vec_print ( m, x, '  Vector X:' )
  call r4vec_print ( n, y, '  Vector Y:' )

  call sger ( m, n, alpha, x, incx, y, incy, a, lda )

  call r4mat_print ( m, n, a, '  Result A = A + alpha * x * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests STRMV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  integer ( kind = 4 ), parameter :: lda = m
  integer ( kind = 4 ), parameter :: n = m

  real ( kind = 4 ) a(lda,n)
  character diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) test
  character trans
  character uplo
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a triangular matrix A,'
  write ( *, '(a)' ) '  STRMV computes y := A * x or y := A'' * x'

  do test = 1, 2

    uplo = 'U'

    if ( test == 1 ) then
      trans = 'N'
    else
      trans = 'T'
    end if

    diag = 'N'

    do j = 1, n
      do i = 1, j
        a(i,j) = real ( i + j, kind = 4 )
      end do
      do i = j + 1, m
        a(i,j) = 0.0E+00
      end do
    end do

    incx = 1
    do i = 1, n
      x(i) = real ( i, kind = 4 )
    end do

    call strmv ( uplo, trans, diag, n, a, lda, x, incx )

    if ( trans == 'N' ) then
      call r4vec_print ( n, x, '  Result y = A * x' );
    else
      call r4vec_print ( n, x, '  Result y = A'' * x' );
    end if

  end do

  return
end
