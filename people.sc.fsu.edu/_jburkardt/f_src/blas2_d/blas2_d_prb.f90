program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS2_D_PRB.
!
!  Discussion:
!
!    BLAS2_D_PRB tests the BLAS library.
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
  write ( *, '(a)' ) 'BLAS2_D_PRB'
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
  write ( *, '(a)' ) 'BLAS2_D_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DGEMV.
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
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character trans
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a general matrix A,'
  write ( *, '(a)' ) '  DGEMV computes y := alpha * A * x + beta * y'
  write ( *, '(a)' ) '  or             y := alpha * A'' * x + beta * y.'
!
!  y = alpha * A * x + beta * y
!
  trans = 'N'
  m = 5
  n = 4
  alpha = 2.0D+00
  lda = m
  allocate ( a(1:m,1:n) )
  call r8mat_test ( trans, lda, m, n, a )
  allocate ( x(1:n) )
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do
  incx = 1
  beta = 3.0D+00
  allocate ( y(1:m) )
  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do
  incy = 1

  call r8mat_print ( m, n, a, '  Matrix A:' )
  call r8vec_print ( n, x, '  Vector X:' )
  call r8vec_print ( m, y, '  Vector Y:' )

  call dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

  call r8vec_print ( m, y, '  Result Y = alpha * A  * x + beta * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  y = alpha * A' * x + beta * y
!
  trans = 'T'
  m = 5
  n = 4
  alpha = 2.0D+00
  lda = m
  allocate ( a(1:m,1:n) )
  call r8mat_test ( trans, lda, n, m, a )
  allocate ( x(1:m) )
  do i = 1, m
    x(i) = real ( i, kind = 8 )
  end do
  incx = 1
  beta = 3.0D+00
  allocate ( y(1:n) )
  do i = 1, n
    y(i) = real ( 10 * i, kind = 8 )
  end do
  incy = 1

  call r8mat_print ( m, n, a, '  Matrix A:' )
  call r8vec_print ( m, x, '  Vector X:' )
  call r8vec_print ( n, y, '  Vector Y:' )

  call dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

  call r8vec_print ( n, y, '  Result Y = alpha * A'' * x + beta * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DGBMV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
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

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character trans
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a general band matrix A,'
  write ( *, '(a)' ) '  DGBMV computes y := alpha * A * x + beta * y'

  trans = 'N'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, m

    jlo = max ( 1, i - kl )
    jhi = min ( n, i + ku )

    do j = jlo, jhi

      if ( i == j ) then
        a(ku+1+i-j,j) = 2.0D+00
      else if ( i == j - 1 .or. i == j + 1 ) then
        a(ku+1+i-j,j) = -1.0D+00
      else
        a(ku+1+i-j,j) = 0.0D+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy )

  call r8vec_print ( m, y, '  Result vector y:' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DSYMV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  character uplo
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a general symmetric matrix A,'
  write ( *, '(a)' ) '  DSYMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( i == j - 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dsymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

  call r8vec_print ( n, y, '  Result vector y := alpha * A * x + beta * y' );

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DSBMV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
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

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character uplo
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a symmetric band matrix A,'
  write ( *, '(a)' ) '  DSBMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, m

    jhi = min ( n, i + k )

    do j = i, jhi

      if ( i == j ) then
        a(k+1+i-j,j) = 2.0D+00
      else if ( i == j - 1 ) then
        a(k+1+i-j,j) = -1.0D+00
      else
        a(k+1+i-j,j) = 0.0D+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dsbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

  call r8vec_print ( m, y, '  Result Y:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DGER.
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

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character trans
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a general matrix A,'
  write ( *, '(a)' ) '  DGER computes A := A + alpha * x * y'''

  m = 5
  n = 4
  alpha = 2.0D+00
  trans = 'N'
  lda = m
  allocate ( a(1:m,1:n) )
  call r8mat_test ( trans, lda, m, n, a )

  allocate ( x(1:m) )
  do i = 1, m
    x(i) = real ( i, kind = 8 )
  end do
  incx = 1

  allocate ( y(1:n) )
  do i = 1, n
    y(i) = real ( 10 * i, kind = 8 )
  end do
  incy = 1

  call r8mat_print ( m, n, a, '  Matrix A:' )
  call r8vec_print ( m, x, '  Vector X:' )
  call r8vec_print ( n, y, '  Vector Y:' )

  call dger ( m, n, alpha, x, incx, y, incy, a, lda )

  call r8mat_print ( m, n, a, '  Result A = A + alpha * x * y' )

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DTRMV.
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

  real ( kind = 8 ) a(lda,n)
  character diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) test
  character trans
  character uplo
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a triangular matrix A,'
  write ( *, '(a)' ) '  DTRMV computes y := A * x or y := A'' * x'

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
        a(i,j) = real ( i + j, kind = 8 )
      end do
      do i = j + 1, m
        a(i,j) = 0.0D+00
      end do
    end do

    incx = 1
    do i = 1, n
      x(i) = real ( i, kind = 8 )
    end do

    call dtrmv ( uplo, trans, diag, n, a, lda, x, incx )

    if ( trans == 'N' ) then
      call r8vec_print ( n, x, '  Result y = A * x' );
    else
      call r8vec_print ( n, x, '  Result y = A'' * x' );
    end if

  end do

  return
end
