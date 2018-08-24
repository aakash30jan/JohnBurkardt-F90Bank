program main

!*****************************************************************************80
!
!! MAIN is the main program for R8PP_PRB.
!
!  Discussion:
!
!    R8PP_PRB tests the R8PP library.
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

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8PP library.'

  call r8pp_det_test ( )
  call r8pp_dif2_test ( )
  call r8pp_fa_test ( )
  call r8pp_indicator_test ( )
  call r8pp_mv_test ( )
  call r8pp_print_test ( )
  call r8pp_print_some_test ( )
  call r8pp_random_test ( )
  call r8pp_sl_test ( )
  call r8pp_to_r8ge_test ( )
  call r8pp_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8pp_det_test ( )

!*****************************************************************************80
!
!! R8PP_DET_TEST tests R8PP_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) det
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_DET_TEST'
  write ( *, '(a)' ) '  R8PP_DET computes the determinant of an R8PP matrix'
  write ( *, '(a)' ) '  factored by R8PP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_dif2 ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8PP_DET_TEST - Warning!'
    write ( *, '(a)' ) '  R8PP_FA failed to factor the matrix.'
    return
  end if
!
!  Compute the determinant.
!
  call r8pp_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant = ', det
  write ( *, '(a,g14.6)' ) '  Exact determinant = ', real ( n + 1, kind = 8 )

  return
end
subroutine r8pp_dif2_test ( )

!*****************************************************************************80
!
!! R8PP_DIF2_TEST tests R8PP_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_DIF2_TEST'
  write ( *, '(a)' ) '  R8PP_DIF2 sets up an R8PP second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_dif2 ( n, a )

  call r8pp_print ( n, a, '  The R8PP second difference matrix:' )
 
  return
end
subroutine r8pp_fa_test ( )

!*****************************************************************************80
!
!! R8PP_FA_TEST tests R8PP_FA.
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

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_FA_TEST'
  write ( *, '(a)' ) '  R8PP_FA factors an R8PP system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The desired solution:' )
!
!  Compute the corresponding right hand side.
!
  call r8pp_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PP_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8PP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The R8PP matrix has been factored.'
  end if
!
!  Solve the linear system.
!
  call r8pp_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8pp_indicator_test ( )

!*****************************************************************************80
!
!! R8PP_INDICATOR_TEST tests R8PP_INDICATOR.
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

  real ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8PP_INDICATOR sets up an R8PP indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP indicator matrix:' )
 
  return
end
subroutine r8pp_mv_test ( )

!*****************************************************************************80
!
!! R8PP_MV_TEST tests R8PP_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_MV_TEST'
  write ( *, '(a)' ) '  R8PP_MV computes b=A*x, where A is an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8pp_indicator ( n, a )
  call r8pp_print ( n, a, '  The R8PP indicator matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )

  call r8pp_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Result b=A*x:' )
 
  return
end
subroutine r8pp_print_test ( )

!*****************************************************************************80
!
!! R8PP_PRINT_TEST tests R8PP_PRINT.
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

  real ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRINT_TEST'
  write ( *, '(a)' ) '  R8PP_PRINT prints an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
 
  return
end
subroutine r8pp_print_some_test ( )

!*****************************************************************************80
!
!! R8PP_PRINT_SOME_TEST tests R8PP_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8PP_PRINT_SOME prints some of an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print_some ( n, a, 2, 3, 6, 5, '  Rows 2-6, Cols 3-5:' )
 
  return
end
subroutine r8pp_random_test ( )

!*****************************************************************************80
!
!! R8PP_RANDOM_TEST tests R8PP_RANDOM.
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

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_RANDOM_TEST'
  write ( *, '(a)' ) '  R8PP_RANDOM, compute a random positive definite'
  write ( *, '(a)' ) '  symmetric packed matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The random R8PP matrix:' )
 
  return
end
subroutine r8pp_sl_test ( )

!*****************************************************************************80
!
!! R8PP_SL_TEST tests R8PP_SL.
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

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_SL_TEST'
  write ( *, '(a)' ) '  R8PP_SL solves a linear system factored by R8PP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The desired solution:' )
!
!  Compute the corresponding right hand side.
!
  call r8pp_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PP_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8PP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The R8PP matrix has been factored.'
  end if
!
!  Solve the linear system.
!
  call r8pp_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8pp_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8PP_TO_R8GE_TEST tests R8PP_TO_R8GE.
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

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8PP_TO_R8GE converts an R8PP matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )

  call r8pp_to_r8ge ( n, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8pp_zeros_test ( )

!*****************************************************************************80
!
!! R8PP_ZEROS_TEST tests R8PP_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_ZEROS_TEST'
  write ( *, '(a)' ) '  R8PP_ZEROS sets a zero R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8pp_zeros ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
 
  return
end
