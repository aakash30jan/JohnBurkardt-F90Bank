program main

!*****************************************************************************80
!
!! MAIN is the main program for R8GE_NP_PRB.
!
!  Discussion:
!
!    R8GE_NP_PRB tests the R8GE_NP library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8GE_NP library.'

  call r8ge_np_det_test ( )
  call r8ge_np_fa_test ( )
  call r8ge_np_inverse_test ( )
  call r8ge_np_ml_test ( )
  call r8ge_np_sl_test ( )
  call r8ge_np_trf_test ( )
  call r8ge_np_trm_test ( )
  call r8ge_np_trs_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8ge_np_det_test ( )

!*****************************************************************************80
!
!! R8GE_NP_DET_TEST tests R8GE_NP_DET.
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

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_DET_TEST'
  write ( *, '(a)' ) '  R8GE_NP_DET computes the determinant of a matrix'
  write ( *, '(a)' ) '  that was factored by R8GE_NP_FA,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_dif2 ( n, n, a )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_DET_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Get the determinant.
!
  call r8ge_np_det ( n, a, det )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant of -1, 2, -1 matrix is ', det
  write ( *, '(a,g14.6)' ) '  Exact value is ', real ( n + 1, kind = 8 )

  return
end
subroutine r8ge_np_fa_test ( )

!*****************************************************************************80
!
!! R8GE_NP_FA_TEST tests R8GE_NP_FA.
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

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_FA_TEST'
  write ( *, '(a)' ) '  R8GE_NP_FA computes the LU factors of a general'
  write ( *, '(a)' ) '  storage matrix without pivoting,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_inverse_test ( )

!*****************************************************************************80
!
!! R8GE_NP_INVERSE_TEST tests R8GE_NP_INVERSE.
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
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_INVERSE'
  write ( *, '(a)' ) '  R8GE_NP_INVERSE computes the inverse of an R8GE matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )

  call r8ge_print ( n, n, a, '  The random matrix:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call r8ge_np_fa ( n, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_INVERSE - Warning!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if

  call r8ge_np_inverse ( n, b )

  call r8ge_print ( n, n, b, '  The inverse matrix:' )
!
!  Compute A * B = C.
!
  call r8ge_mm ( n, n, n, a, b, c )

  call r8ge_print ( n, n, c, '  The product:' )

  return
end
subroutine r8ge_np_ml_test ( )

!*****************************************************************************80
!
!! R8GE_NP_ML_TEST tests R8GE_NP_ML.
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

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_ML_TEST'
  write ( *, '(a)' ) '  R8GE_NP_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  where A has been factored by R8GE_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8ge_mv ( n, n, a, x, b )
    else
      call r8ge_mtv ( n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8ge_np_fa ( n, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_NP_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8ge_np_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8ge_np_sl_test ( )

!*****************************************************************************80
!
!! R8GE_NP_SL_TEST tests R8GE_NP_SL.
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

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_SL_TEST'
  write ( *, '(a)' ) '  R8GE_NP_SL solves a linear system that was factored'
  write ( *, '(a)' ) '  by R8GE_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trf_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRF_TEST tests R8GE_NP_TRF.
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

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: nrhs = 1

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRF_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRF factors an R8GE matrix without pivoting,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trm_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRM_TEST tests R8GE_NP_TRM.
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

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: nrhs = 1

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRM_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRM computes b=A*x after A has'
  write ( *, '(a)' ) '  been factored by R8GE_NP_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trs_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRS_TEST tests R8GE_NP_TRS.
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

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: nrhs = 1

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRS_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRS solves a linear system A*x=b'
  write ( *, '(a)' ) '  which has been factored by R8GE_NP_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end

