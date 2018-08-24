program main

!*****************************************************************************80
!
!! MAIN is the main program for R8PBL_PRB.
!
!  Discussion:
!
!    R8PBL_PRB tests the R8PBL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8PBL library.'

  call r8pbl_dif2_test ( )
  call r8pbl_indicator_test ( )
  call r8pbl_mv_test ( )
  call r8pbl_print_test ( )
  call r8pbl_print_some_test ( )
  call r8pbl_random_test ( )
  call r8pbl_to_r8ge_test ( )
  call r8pbl_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8pbl_dif2_test ( )

!*****************************************************************************80
!
!! R8PBL_DIF2_TEST tests R8PBL_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ml = 1

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_DIF2_TEST'
  write ( *, '(a)' ) '  R8PBL_DIF2 sets up an R8PBL indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  call r8pbl_dif2 ( n, ml, a )

  call r8pbl_print ( n, ml, a, '  The R8PBL second difference matrix:' )

  return
end
subroutine r8pbl_indicator_test ( )

!*****************************************************************************80
!
!! R8PBL_INDICATOR_TEST tests R8PBL_INDICATOR.
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
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8PBL_INDICATOR sets up an R8PBL indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  call r8pbl_indicator ( n, ml, a )

  call r8pbl_print ( n, ml, a, '  The R8PBL indicator matrix:' )

  return
end
subroutine r8pbl_mv_test ( )

!*****************************************************************************80
!
!! R8PBL_MV_TEST tests R8PBL_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_MV_TEST'
  write ( *, '(a)' ) '  R8PBL_MV computes A*x, where A is an R8PBL matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
!
!  Set the matrix.
!
  call r8pbl_random ( n, ml, seed, a )
  call r8pbl_print ( n, ml, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )
!
!  Compute the corresponding right hand side.
!
  call r8pbl_mv ( n, ml, a, x, b )
  call r8vec_print ( n, b, '  Product b=A*x' )

  return
end
subroutine r8pbl_print_test ( )

!*****************************************************************************80
!
!! R8PBL_PRINT_TEST tests R8PBL_PRINT.
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
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_PRINT_TEST'
  write ( *, '(a)' ) '  R8PBL_PRINT prints an R8PBL matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  call r8pbl_indicator ( n, ml, a )

  call r8pbl_print ( n, ml, a, '  The R8PBL matrix:' )

  return
end
subroutine r8pbl_print_some_test ( )

!*****************************************************************************80
!
!! R8PBL_PRINT_SOME_TEST tests R8PBL_PRINT_SOME.
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

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: ml = 4

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8PBL_PRINT_SOME prints some of an R8PBL matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  call r8pbl_indicator ( n, ml, a )

  call r8pbl_print_some ( n, ml, a, 4, 5, 8, 9, '  Row(4:8), Col(5:9):' )

  return
end
subroutine r8pbl_random_test ( )

!*****************************************************************************80
!
!! R8PBL_RANDOM_TEST tests R8PBL_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_RANDOM_TEST'
  write ( *, '(a)' ) '  R8PBL_RANDOM sets up a random R8PBL matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  seed = 123456789
  call r8pbl_random ( n, ml, seed, a )

  call r8pbl_print ( n, ml, a, '  The random R8PBL matrix:' )

  return
end
subroutine r8pbl_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8PBL_TO_R8GE_TEST tests R8PBL_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8PBL_TO_R8GE converts an R8PBL matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  seed = 123456789
  call r8pbl_random ( n, ml, seed, a )

  call r8pbl_print ( n, ml, a, '  The R8PBL matrix:' )

  call r8pbl_to_r8ge ( n, ml, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8pbl_zeros_test ( )

!*****************************************************************************80
!
!! R8PBL_ZEROS_TEST tests R8PBL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ml = 2

  real ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBL_ZEROS_TEST'
  write ( *, '(a)' ) '  R8PBL_ZEROS sets up an R8PBL zero matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth ML =   ', ml

  call r8pbl_zeros ( n, ml, a )

  call r8pbl_print ( n, ml, a, '  The R8PBL indicator matrix:' )

  return
end
