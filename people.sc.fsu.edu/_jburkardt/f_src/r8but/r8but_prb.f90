program main

!*****************************************************************************80
!
!! MAIN is the main program for R8BUT_PRB.
!
!  Discussion:
!
!    R8BUT_PRB tests the R8BUT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8BUT library.'

  call r8but_det_test ( )
  call r8but_indicator_test ( )
  call r8but_mtv_test ( )
  call r8but_mv_test ( )
  call r8but_print_test ( )
  call r8but_print_some_test ( )
  call r8but_random_test ( )
  call r8but_sl_test ( )
  call r8but_slt_test ( )
  call r8but_to_r8ge_test ( )
  call r8but_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8but_det_test ( )

!*****************************************************************************80
!
!! R8BUT_DET_TEST tests R8BUT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_DET_TEST'
  write ( *, '(a)' ) '  R8BUT_DET computes the determinant of an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8but_det ( n, mu, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant = ', det

  return
end
subroutine r8but_indicator_test ( )

!*****************************************************************************80
!
!! R8BUT_INDICATOR_TEST tests R8BUT_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8BUT_INDICATOR sets up an R8BUT indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8but_indicator ( n, mu, a )

  call r8but_print ( n, mu, a, '  The R8BUT indicator matrix:' )

  return
end
subroutine r8but_mtv_test ( )

!*****************************************************************************80
!
!! R8BUT_MTV_TEST tests R8BUT_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_MTV_TEST'
  write ( *, '(a)' ) '  R8BUT_MTV computes b=A''*x, where A is an R8BUT matrix.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8but_mtv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8but_mv_test ( )

!*****************************************************************************80
!
!! R8BUT_MV_TEST tests R8BUT_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_MV_TEST'
  write ( *, '(a)' ) '  R8BUT_MV computes b=A*x, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8but_mv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8but_print_test ( )

!*****************************************************************************80
!
!! R8BUT_PRINT_TEST tests R8BUT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) :: seed = 123456789
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRINT_TEST'
  write ( *, '(a)' ) '  R8BUT_PRINT prints an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  return
end
subroutine r8but_print_some_test ( )

!*****************************************************************************80
!
!! R8BUT_PRINT_SOME_TEST tests R8BUT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8BUT_PRINT_SOME prints some of an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_indicator ( n, mu, a )

  call r8but_print_some ( n, mu, a, 1, 2, 4, 4, '  Rows 1:4, Cols 2:4:' )

  return
end
subroutine r8but_random_test ( )

!*****************************************************************************80
!
!! R8BUT_RANDOM_TEST tests R8BUT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) :: seed = 123456789
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_RANDOM_TEST'
  write ( *, '(a)' ) '  R8BUT_RANDOM randomizes an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  return
end
subroutine r8but_sl_test ( )

!*****************************************************************************80
!
!! R8BUT_SL_TEST tests R8BUT_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_SL_TEST'
  write ( *, '(a)' ) '  R8BUT_SL solves A*x=b, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8but_mv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b:' )
!
!  Solve the linear system.
!
  call r8but_sl ( n, mu, a, b )
  call r8vec_print ( n, b, '  x:' )

  return
end
subroutine r8but_slt_test ( )

!*****************************************************************************80
!
!! R8BUT_SLT_TEST tests R8BUT_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_SLT_TEST'
  write ( *, '(a)' ) '  R8BUT_SLT solves A''*x=b, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8but_mtv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b:' )
!
!  Solve the linear system.
!
  call r8but_slt ( n, mu, a, b )
  call r8vec_print ( n, b, '  x:' )

  return
end
subroutine r8but_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8BUT_TO_R8GE_TEST tests R8BUT_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a_r8but(mu+1,n)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8BUT_TO_R8GE converts a matrix from R8BUT to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a_r8but )

  call r8but_print ( n, mu, a_r8but, '  The R8BUT matrix:' )

  call r8but_to_r8ge ( n, mu, a_r8but, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix' )

  return
end
subroutine r8but_zeros_test ( )

!*****************************************************************************80
!
!! R8BUT_ZEROS_TEST tests R8BUT_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_ZEROS_TEST'
  write ( *, '(a)' ) '  R8BUT_ZEROS zeros an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_zeros ( n, mu, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  return
end

