program main

!*****************************************************************************80
!
!! MAIN is the main program for R8GD_PRB.
!
!  Discussion:
!
!    R8GD_PRB tests the R8GD library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8GD library.'

  call r8gd_dif2_test ( )
  call r8gd_indicator_test ( )
  call r8gd_mtv_test ( )
  call r8gd_mv_test ( )
  call r8gd_print_test ( )
  call r8gd_print_some_test ( )
  call r8gd_random_test ( )
  call r8gd_to_r8ge_test ( )
  call r8gd_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8gd_dif2_test ( )

!*****************************************************************************80
!
!! R8GD_DIF2_TEST tests R8GD_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: ndiag = 3

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_DIF2_TEST'
  write ( *, '(a)' ) '  R8GD_DIF2 sets up an R8GD second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag

  call r8gd_dif2 ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The R8GS second difference matrix:' )

  return
end
subroutine r8gd_indicator_test ( )

!*****************************************************************************80
!
!! R8GD_INDICATOR_TEST tests R8GD_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8GD_INDICATOR sets up an indicator matrix'
  write ( *, '(a)' ) '  for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_mtv_test ( )

!*****************************************************************************80
!
!! R8GD_MTV_TEST tests R8GD_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_MTV_TEST'
  write ( *, '(a)' ) '  R8GD_MTV computes A'' * x for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8gd_mtv ( n, ndiag, offset, a, x, b )
  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8gd_mv_test ( )

!*****************************************************************************80
!
!! R8GD_MV_TEST tests R8GD_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_MV_TEST'
  write ( *, '(a)' ) '  R8GD_MV computes A * x for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8gd_mv ( n, ndiag, offset, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8gd_print_test ( )

!*****************************************************************************80
!
!! R8GD_PRINT_TEST tests R8GD_PRINT.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRINT_TEST'
  write ( *, '(a)' ) '  R8GD_PRINT prints a general diagonal matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_print_some_test ( )

!*****************************************************************************80
!
!! R8GD_PRINT_SOME_TEST tests R8GD_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8GD_PRINT_SOME prints some of an R8GD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print_some ( n, ndiag, offset, a, 3, 3, 6, 6, '  Rows 3-6, Cols 3-6:' )

  return
end
subroutine r8gd_random_test ( )

!*****************************************************************************80
!
!! R8GD_RANDOM_TEST tests R8GD_RANDOM.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_RANDOM_TEST'
  write ( *, '(a)' ) '  R8GD_RANDOM randomly generates a general diagonal matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8GD_TO_R8GE_TEST tests R8GD_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8GD_TO_R8GE converts an R8GD matrix to R8GD format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The R8GD matrix:' )

  call r8gd_to_r8ge ( n, ndiag, offset, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8gd_zeros_test ( )

!*****************************************************************************80
!
!! R8GD_ZEROS_TEST tests R8GD_ZEROS.
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
  integer ( kind = 4 ), parameter :: ndiag = 4

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_ZEROS_TEST'
  write ( *, '(a)' ) '  R8GD_ZEROS returns a zero R8GD matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_zeros ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The zero R8GD matrix:' )

  return
end

