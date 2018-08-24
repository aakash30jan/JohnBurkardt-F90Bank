program main

!*****************************************************************************80
!
!! MAIN is the main program for R8BLT_PRB.
!
!  Discussion:
!
!    R8BLT_PRB tests the R8BLT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8BLT library.'

  call r8blt_det_test ( )
  call r8blt_indicator_test ( )
  call r8blt_mtv_test ( )
  call r8blt_mv_test ( )
  call r8blt_print_test ( )
  call r8blt_print_some_test ( )
  call r8blt_random_test ( )
  call r8blt_sl_test ( )
  call r8blt_slt_test ( )
  call r8blt_to_r8ge_test ( )
  call r8blt_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8blt_det_test ( )

!*****************************************************************************80
!
!! R8BLT_DET_TEST tests R8BLT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_DET_TEST'
  write ( *, '(a)' ) '  R8BLT_DET computes the determinant of an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  call r8blt_det ( n, ml, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant = ', det

  return
end
subroutine r8blt_indicator_test ( )

!*****************************************************************************80
!
!! R8BLT_INDICATOR_TEST tests R8BLT_INDICATOR.
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
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8BLT_INDICATOR sets up an R8BLT indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
!
!  Set the matrix.
!
  call r8blt_indicator ( n, ml, a )

  call r8blt_print ( n, ml, a, '  The R8BLT indicator matrix:' )

  return
end
subroutine r8blt_mtv_test ( )

!*****************************************************************************80
!
!! R8BLT_MTV_TEST tests R8BLT_MTV.
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

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_MTV_TEST'
  write ( *, '(a)' ) '  R8BLT_MTV computes b=A''*x, where A is an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )
  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8blt_mtv ( n, ml, a, x, b )
  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8blt_mv_test ( )

!*****************************************************************************80
!
!! R8BLT_MV_TEST tests R8BLT_MV.
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

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_MV_TEST'
  write ( *, '(a)' ) '  R8BLT_MV computes b=A*x, where A is an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )
  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8blt_mv ( n, ml, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8blt_print_test ( )

!*****************************************************************************80
!
!! R8BLT_PRINT_TEST tests R8BLT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_PRINT_TEST'
  write ( *, '(a)' ) '  R8BLT_PRINT prints an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  return
end
subroutine r8blt_print_some_test ( )

!*****************************************************************************80
!
!! R8BLT_PRINT_SOME_TEST tests R8BLT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8BLT_PRINT_SOME prints some of an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_indicator ( n, ml, a )

  call r8blt_print_some ( n, ml, a, 2, 3, 5, 5, '  Rows 2:5, Cols 3:5:' )

  return
end
subroutine r8blt_random_test ( )

!*****************************************************************************80
!
!! R8BLT_RANDOM_TEST tests R8BLT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_RANDOM_TEST'
  write ( *, '(a)' ) '  R8BLT_RANDOM randomizes an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  return
end
subroutine r8blt_sl_test ( )

!*****************************************************************************80
!
!! R8BLT_SL_TEST tests R8BLT_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_SL_TEST'
  write ( *, '(a)' ) '  R8BLT_SL solves A*x=b where A is an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8blt_mv ( n, ml, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Solve the linear system.
!
  call r8blt_sl ( n, ml, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8blt_slt_test ( )

!*****************************************************************************80
!
!! R8BLT_SLT_TEST tests R8BLT_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_SLT_TEST'
  write ( *, '(a)' ) '  R8BLT_SLT solves A''*x=b where A is an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8blt_mtv ( n, ml, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Solve the linear system.
!
  call r8blt_slt ( n, ml, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8blt_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8BLT_TO_R8GE_TEST tests R8BLT_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a_r8blt(ml+1,n)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8BLT_TO_R8GE converts a matrix from R8BLT to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a_r8blt )

  call r8blt_print ( n, ml, a_r8blt, '  The R8BLT matrix:' )

  call r8blt_to_r8ge ( n, ml, a_r8blt, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix' )

  return
end
subroutine r8blt_zeros_test ( )

!*****************************************************************************80
!
!! R8BLT_ZEROS_TEST tests R8BLT_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BLT_ZEROS_TEST'
  write ( *, '(a)' ) '  R8BLT_ZEROS zeros an R8BLT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_zeros ( n, ml, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  return
end

