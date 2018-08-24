program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPLUS_C8_PRB.
!
!  Discussion:
!
!    LINPLUS_C8_PRB tests the LINPLUS_C8 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPLUS_C8_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the LINPLUS_C8 library.'

  call c83_cr_fa_test ( )
  call c83_cr_sls_test ( )

  call c83_np_fa_test ( )
  call c83_np_ml_test ( )
  call c83_np_sl_test ( )

  call c8ci_sl_test ( )
  call c8to_sl_test ( )
  call c8vec_unity_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPLUS_C8_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine c83_cr_fa_test ( )

!*****************************************************************************80
!
!! C83_CR_FA_TEST tests C83_CR_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: nb = 1

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) a_cr(3,0:2*n)
  complex ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n,nb)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C83_CR_FA_TEST'
  write ( *, '(a)' ) '  C83_CR_FA factors a complex tridiagonal matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  do j = 2, n
    a(1,j) = cmplx ( -1.0D+00, - real ( j - 1, kind = 8 ) )
  end do

  do j = 1, n
    a(2,j) = cmplx ( 2.0D+00, 2.0D+00 * real ( j, kind = 8 ) )
  end do

  do j = 1, n - 1
    a(3,j) = cmplx ( -1.0D+00, - real ( j + 1, kind = 8 ) )
  end do
  a(3,n) = 0.0D+00
!
!  Set the desired solution.
!
  do i = 1, n
    x(i,1) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  i = 1
  b(1,1) = cmplx ( 2.0D+00, 2.0D+00 * real ( i, kind = 8 ) ) * x(i,1) &
         + cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i+1,1)

  do i = 2, n - 1
    b(i,1) = cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i-1,1) &
           + cmplx (  2.0D+00, real ( 2*i, kind = 8  ) ) * x(i,1) &
           + cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i+1,1)
  end do

  b(n,1) = cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i-1,1) &
         + cmplx (  2.0D+00, real ( 2*i, kind = 8 ) ) * x(i,1)
!
!  Factor the matrix.
!
  call c83_cr_fa ( n, a, a_cr )
!
!  Solve the linear system.
!
  call c83_cr_sls ( n, a_cr, nb, b, x )

  call c8ge_print_some ( n, nb, x, 1, 1, 10, nb, '  Solution:' )

  return
end
subroutine c83_cr_sls_test ( )

!*****************************************************************************80
!
!! C83_CR_SLS_TEST tests C83_CR_SLS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: nb = 1

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) a_cr(3,0:2*n)
  complex ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n,nb)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C83_CR_SLS_TEST'
  write ( *, '(a)' ) '  C83_CR_SLS solves one or more linear systems'
  write ( *, '(a)' ) '  that were factored by C83_CR_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  do j = 2, n
    a(1,j) = cmplx ( -1.0D+00, - real ( j - 1, kind = 8 ) )
  end do

  do j = 1, n
    a(2,j) = cmplx ( 2.0D+00, 2.0D+00 * real ( j, kind = 8 ) )
  end do

  do j = 1, n - 1
    a(3,j) = cmplx ( -1.0D+00, - real ( j + 1, kind = 8 ) )
  end do
  a(3,n) = 0.0D+00
!
!  Set the desired solution.
!
  do i = 1, n
    x(i,1) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  i = 1
  b(1,1) = cmplx ( 2.0D+00, 2.0D+00 * real ( i, kind = 8 ) ) * x(i,1) &
         + cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i+1,1)

  do i = 2, n - 1
    b(i,1) = cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i-1,1) &
           + cmplx (  2.0D+00, real ( 2*i, kind = 8  ) ) * x(i,1) &
           + cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i+1,1)
  end do

  b(n,1) = cmplx ( -1.0D+00, real ( -i, kind = 8 ) ) * x(i-1,1) &
         + cmplx (  2.0D+00, real ( 2*i, kind = 8 ) ) * x(i,1)
!
!  Factor the matrix.
!
  call c83_cr_fa ( n, a, a_cr )
!
!  Solve the linear system.
!
  call c83_cr_sls ( n, a_cr, nb, b, x )

  call c8ge_print_some ( n, nb, x, 1, 1, 10, nb, '  Solution:' )

  return
end
subroutine c83_np_fa_test ( )

!*****************************************************************************80
!
!! C83_NP_FA_TEST tests C83_NP_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C83_NP_FA_TEST'
  write ( *, '(a)' ) '  C83_NP_FA factors a C83 matrix with no pivoting;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c83_random ( n, seed, a )

  call c83_print ( n, a, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_mv ( n, a, x, b )

  call c8vec_print ( n, b, '  The right hand side' )
!
!  Factor the matrix.
!
  call c83_np_fa ( n, a, info )
!
!  Solve the linear system.
!
  job = 0
  call c83_np_sl ( n, a, b, job )

  call c8vec_print ( n, b, '  The solution' )
!
!  Now set a SECOND desired solution.
!
  do i = 1, n
    x(i) = cmplx ( real ( 10 * i, kind = 8 ), real ( i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side, using the FACTORED matrix.
!
  call c83_np_ml ( n, a, x, b, job )

  call c8vec_print ( n, b, '  The second right hand side' )
!
!  Solve the linear system.
!
  call c83_np_sl ( n, a, b, job )

  call c8vec_print ( n, b, '  The second solution' )

  return
end
subroutine c83_np_ml_test ( )

!*****************************************************************************80
!
!! C83_NP_ML_TEST tests C83_NP_ML.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C83_NP_ML_TEST'
  write ( *, '(a)' ) '  C83_NP_ML multiplies A*X'
  write ( *, '(a)' ) '  after A has been factored by C83_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c83_random ( n, seed, a )

  call c83_print ( n, a, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_mv ( n, a, x, b )

  call c8vec_print ( n, b, '  The right hand side' )
!
!  Factor the matrix.
!
  call c83_np_fa ( n, a, info )
!
!  Solve the linear system.
!
  job = 0
  call c83_np_sl ( n, a, b, job )

  call c8vec_print ( n, b, '  The solution' )
!
!  Now set a SECOND desired solution.
!
  do i = 1, n
    x(i) = cmplx ( real ( 10 * i, kind = 8 ), real ( i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side, using the FACTORED matrix.
!
  call c83_np_ml ( n, a, x, b, job )

  call c8vec_print ( n, b, '  The second right hand side' )
!
!  Solve the linear system.
!
  call c83_np_sl ( n, a, b, job )

  call c8vec_print ( n, b, '  The second solution' )

  return
end
subroutine c83_np_sl_test ( )

!*****************************************************************************80
!
!! C83_NP_SL_TEST tests C83_NP_SL.
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

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C83_NP_SL_TEST'
  write ( *, '(a)' ) '  C83_NP_SL solves a linear system factored by C83_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will look at the TRANSPOSED linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c83_random ( n, seed, a )

  call c83_print ( n, a, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_mtv ( n, a, x, b )

  call c8vec_print ( n, b, '  The right hand side B1' )
!
!  Factor the matrix.
!
  call c83_np_fa ( n, a, info )
!
!  Solve the linear system.
!
  job = 1
  call c83_np_sl ( n, a, b, job )
 
  call c8vec_print ( n, b, '  The solution to At * X1 = B1' )

  return
end
subroutine c8ci_sl_test ( )

!*****************************************************************************80
!
!! C8CI_SL_TEST tests C8CI_SL.
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

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8CI_SL_TEST'
  write ( *, '(a)' ) '  C8CI_SL solves a complex circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c8ci_random ( n, seed, a )

  call c8ci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    do i = 1, n
      x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
    end do
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call c8ci_mv ( n, a, x, b )
    else
      call c8ci_mtv ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call c8ci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call c8vec_print ( n, x, '  Solution:' )
    else
      call c8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine c8to_sl_test ( )

!*****************************************************************************80
!
!! C8TO_SL_TEST tests C8TO_SL.
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

  integer ( kind = 4 ), parameter :: n = 4

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8TO_SL_TEST'
  write ( *, '(a)' ) '  C8TO_SL solves a complex Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c8to_random ( n, seed, a )

  call c8to_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call c8vec_indicator ( n, x )

    call c8vec_print_some ( n, x, 10, '  Desired solution:' )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call c8to_mv ( n, a, x, b )
    else
      call c8to_mtv ( n, a, x, b )
    end if

    call c8vec_print_some ( n, b, 10, '  Right hand side:' )
!
!  Solve the linear system.
!
    call c8to_sl ( n, a, b, x, job )

    call c8vec_print_some ( n, x, 10, '  Computed solution:' )

  end do

  return
end
subroutine c8vec_unity_test ( )

!*****************************************************************************80
!
!! C8VEC_UNITY_TEST tests C8VEC_UNITY.
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

  integer ( kind = 4 ), parameter :: n_max = 5

  integer ( kind = 4 ) n
  complex ( kind = 8 ) x(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8VEC_UNITY_TEST'
  write ( *, '(a)' ) '  C8VEC_UNITY returns the N complex roots of unity.'

  do n = 1, n_max

    call c8vec_unity ( n, x )

    call c8vec_print ( n, x, '  Roots of unity:' )

  end do

  return
end


