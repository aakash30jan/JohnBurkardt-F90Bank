program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPACK_Q_PRB.
!
!  Discussion:
!
!    LINPACK_D_PRB tests the LINPACK_Q library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_Q_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINPACK_Q library.'

  call qgeco_test ( )
  call qgedi_test ( )
  call qgefa_test ( )
  call qgesl_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_Q_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine qgeco_test ( )

!*****************************************************************************80
!
!! QGECO_TEST calls QGECO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    LDA defines the maximum matrix size we will use.
!
  implicit none

  integer, parameter :: qp = selected_real_kind ( 33, 4931 )

  integer ( kind = 4 ), parameter :: lda = 10

  real ( kind = qp ) a(lda,lda)
  real ( kind = qp ) b(lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(lda)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) n
  real ( kind = qp ) rcond
  real ( kind = qp ) z(lda)

  n = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QGECO_TEST'
  write ( *, '(a)' ) '  QGECO factors a general matrix and computes'
  write ( *, '(a)' ) '  its reciprocal condition number;'
  write ( *, '(a)' ) '  QGESL can be called afterwards to solve a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0Q+00
  a(1,2) = 2.0Q+00
  a(1,3) = 3.0Q+00

  a(2,1) = 4.0Q+00
  a(2,2) = 5.0Q+00
  a(2,3) = 6.0Q+00

  a(3,1) = 7.0Q+00
  a(3,2) = 8.0Q+00
  a(3,3) = 0.0Q+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call qgeco ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a,g14.6)' ) '  The reciprocal matrix condition number = ', rcond

  if ( rcond + 1.0Q+00 == 1.0Q+00 ) then
    write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
    return
  end if
!
!  Set a right hand side.
!
  b(1:3) = (/ 14.0Q+00, 32.0Q+00, 23.0Q+00 /)
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call qgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by QGESL'
  write ( *, '(a)' ) '  (Should be (1,2,3))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  A second right hand side can be solved without refactoring a.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call QGESL for a new right hand '
  write ( *, '(a)' ) '  side for the same, factored matrix.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 1.0Q+00, 4.0Q+00, 7.0Q+00 /)
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve a linear system.'

  job = 0
  call qgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by QGESL'
  write ( *, '(a)' ) '  (should be (1,0,0))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  The transposed problem A'*X = B can be solved by QGESL
!  as well, without any refactoring.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call QGESL for transposed problem.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 6.0Q+00, 6.0Q+00, -3.0Q+00 /)
!
!  Solve the transposed system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call QGESL to solve a transposed linear system.'

  job = 1
  call qgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by QGESL'
  write ( *, '(a)' ) '  (should be (-1,0,1))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do

  return
end
subroutine qgedi_test ( )

!*****************************************************************************80
!
!! QGEDI_TEST tests QGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: qp = selected_real_kind ( 33, 4931 )

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = qp ) a(lda,n)
  real ( kind = qp ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  real ( kind = qp ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QGEDI_TEST'
  write ( *, '(a)' ) '  QGEDI computes the inverse and determinant'
  write ( *, '(a)' ) '  of a matrix factored by QGEFA.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0Q+00
  a(1,2) = 2.0Q+00
  a(1,3) = 3.0Q+00

  a(2,1) = 4.0Q+00
  a(2,2) = 5.0Q+00
  a(2,3) = 6.0Q+00

  a(3,1) = 7.0Q+00
  a(3,2) = 8.0Q+00
  a(3,3) = 0.0Q+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix'

  call qgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
    return
  end if
!
!  Get the inverse and determinant.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the inverse and determinant'

  job = 11
  call qgedi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  The determinant = ', det(1), ' * 10 ** ', det(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inverse matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  return
end
subroutine qgefa_test ( )

!*****************************************************************************80
!
!! QGEFA_TEST tests QGEFA.
!
!  Discussion:
!
!    Solve A*x = b where A is a given matrix, and B a right hand side.
!
!    We will also assume that A is stored in the simplest
!    possible way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: qp = selected_real_kind ( 33, 4931 )

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = qp ) a(lda,n)
  real ( kind = qp ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QGEFA_TEST'
  write ( *, '(a)' ) '  QGEFA factors a general matrix,'
  write ( *, '(a)' ) '  which can then be solved by QGESL.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0Q+00
  a(1,2) = 2.0Q+00
  a(1,3) = 3.0Q+00

  a(2,1) = 4.0Q+00
  a(2,2) = 5.0Q+00
  a(2,3) = 6.0Q+00

  a(3,1) = 7.0Q+00
  a(3,2) = 8.0Q+00
  a(3,3) = 0.0Q+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  b(1:3) = (/ 14.0Q+00, 32.0Q+00, 23.0Q+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix'

  call qgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  QGEFA returned an error flag INFO = ', info
    return
  end if
!
!  If no error occurred, now use QGESL to solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call qgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QGESL returns the solution:'
  write ( *, '(a)' ) '  (Should be (1,2,3))'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do

  return
end
subroutine qgesl_test ( )

!*****************************************************************************80
!
!! QGESL_TEST tests QGESL.
!
!  Discussion:
!
!    In this example, we solve a relatively large linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: qp = selected_real_kind ( 33, 4931 )

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = qp ) a(lda,n)
  real ( kind = qp ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QGESL_TEST'
  write ( *, '(a)' ) '  QGESL solves a factored linear system after it'
  write ( *, '(a)' ) '  has been factored by QGEFA.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
!  The problem is just an enlarged version of the
!  problem for N = 5, which is:
!
!  Matrix A is ( n -1 -1 -1 -1)    Right hand side B is  (1)
!              (-1  n -1 -1 -1)                          (1)
!              (-1 -1  n -1 -1)                          (1)
!              (-1 -1 -1  n -1)                          (1)
!              (-1 -1 -1 -1  n)                          (1)
!
!  Solution is   (1)
!                (1)
!                (1)
!                (1)
!                (1)
!
  b(1:n) = 1.0Q+00

  a(1:n,1:n) = -1.0Q+00
  do i = 1, n
    a(i,i) = real ( n, kind = qp )
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call qgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  QGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call qgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last five solution entries:'
  write ( *, '(a)' ) '  (All of them should be 1.)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    if ( i <= 5 .or. n - 5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end