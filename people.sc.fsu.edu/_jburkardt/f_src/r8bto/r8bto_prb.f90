program main

!*****************************************************************************80
!
!! MAIN is the main program for R8BTO_PRB.
!
!  Discussion:
!
!    R8BTO_PRB tests the R8BTO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8BTO library.'

  call r8bto_dif2_test ( )
  call r8bto_indicator_test ( )
  call r8bto_mtv_test ( )
  call r8bto_mv_test ( )
  call r8bto_print_test ( )
  call r8bto_print_some_test ( )
  call r8bto_random_test ( )
  call r8bto_sl_test ( )
  call r8bto_to_r8ge_test ( )
  call r8bto_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8bto_dif2_test ( )

!*****************************************************************************80
!
!! R8BTO_DIF2_TEST tests R8BTO_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 5
  integer ( kind = 4 ), parameter :: m = 1

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_DIF2_TEST'
  write ( *, '(a)' ) '  R8BTO_DIF2 sets up an R8BTO version of'
  write ( *, '(a)' ) '  the second difference matrix (assuming M = 1 ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_dif2 ( m, l, a )

  call r8bto_print ( m, l, a, '  The R8BTO second difference matrix:' )

  return
end
subroutine r8bto_indicator_test ( )

!*****************************************************************************80
!
!! R8BTO_INDICATOR_TEST tests R8BTO_INDICATOR.
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

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8BTO_INDICATOR sets up an indicator matrix'
  write ( *, '(a)' ) '  for a real block Toeplitz matrix,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_indicator ( m, l, a )

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  return
end
subroutine r8bto_mtv_test ( )

!*****************************************************************************80
!
!! R8BTO_MTV_TEST tests R8BTO_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ), parameter :: n = m * l
  integer ( kind = 4 ), parameter :: p = 2 * l - 1

  real ( kind = 8 ), dimension ( m, m, p ) ::  a = reshape ( (/ &
    1.0D+00, 5.0D+00, 2.0D+00, 5.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00, &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, p /) )
  real ( kind = 8 ) b(m,l)
  real ( kind = 8 ) x(m,l)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_MTV_TEST'
  write ( *, '(a)' ) '  R8BTO_MTV computes A''* x for a block Toeplitz matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  call r8ge_indicator ( m, l, x )

  call r8ge_print ( m, l, x, '  The "vector" x:' )

  call r8bto_mtv ( m, l, a, x, b )

  call r8ge_print ( m, l, b, '  The product A''*x:' )

  return
end
subroutine r8bto_mv_test ( )

!*****************************************************************************80
!
!! R8BTO_MV_TEST tests R8BTO_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ), parameter :: n = m * l
  integer ( kind = 4 ), parameter :: p = 2 * l - 1

  real ( kind = 8 ), dimension ( m, m, p ) ::  a = reshape ( (/ &
    1.0D+00, 5.0D+00, 2.0D+00, 5.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00, &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, p /) )
  real ( kind = 8 ) b(m,l)
  real ( kind = 8 ) x(m,l)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_MV_TEST'
  write ( *, '(a)' ) '  R8BTO_MV computes A * x for a block Toeplitz matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  call r8ge_indicator ( m, l, x )

  call r8ge_print ( m, l, x, '  The "vector" x:' )

  call r8bto_mv ( m, l, a, x, b )

  call r8ge_print ( m, l, b, '  The product A*x:' )

  return
end
subroutine r8bto_print_test ( )

!*****************************************************************************80
!
!! R8BTO_PRINT_TEST tests R8BTO_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_PRINT_TEST'
  write ( *, '(a)' ) '  R8BTO_PRINT prints an R8BTO matrix,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_indicator ( m, l, a )

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  return
end
subroutine r8bto_print_some_test ( )

!*****************************************************************************80
!
!! R8BTO_PRINT_SOME_TEST tests R8BTO_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8BTO_PRINT_SOME prints some of an R8BTO matrix,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_indicator ( m, l, a )

  call r8bto_print_some ( m, l, a, 1, 3, 6, 4, '  Row (1:6), Cols (3:4):' )

  return
end
subroutine r8bto_random_test ( )

!*****************************************************************************80
!
!! R8BTO_RANDOM_TEST tests R8BTO_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) :: a
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_RANDOM_TEST'
  write ( *, '(a)' ) '  R8BTO_RANDOM returns a random R8BTO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  seed = 123456789
  call r8bto_random ( m, l, seed, a )

  call r8bto_print ( m, l, a, '  The random R8BTO matrix:' )

  return
end
subroutine r8bto_sl_test ( )

!*****************************************************************************80
!
!! R8BTO_SL_TEST tests R8BTO_SL.
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

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: l = 3

  integer ( kind = 4 ), parameter :: n = m * l
  integer ( kind = 4 ), parameter :: p = 2 * l - 1

  real ( kind = 8 ), dimension ( m, m, p ) ::  a = reshape ( (/ &
    9.0D+00, 2.0D+00, 1.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00, &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, p /) )
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_SL_TEST'
  write ( *, '(a)' ) '  R8BTO_SL solves a block Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  call r8bto_mv ( m, l, a, x, b )

  call r8vec_print ( n, b, '  The right hand side B:' )

  call r8bto_sl ( m, l, a, b, x )

  call r8vec_print ( n, x, '  The computed solution X:' )

  return
end
subroutine r8bto_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8BTO_TO_R8GE_TEST tests R8BTO_TO_R8GE.
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

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ), parameter :: n = l * m

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) :: a
  real ( kind = 8 ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8BTO_TO_R8GE converts an R8BTO matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_indicator ( m, l, a )

  call r8bto_print ( m, l, a, '  The R8BTO matrix:' )

  call r8bto_to_r8ge ( m, l, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8bto_zeros_test ( )

!*****************************************************************************80
!
!! R8BTO_ZEROS_TEST tests R8BTO_ZEROS.
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

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BTO_ZEROS_TEST'
  write ( *, '(a)' ) '  R8BTO_ZEROS zeros an R8BTO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_zeros ( m, l, a )

  call r8bto_print ( m, l, a, '  The zero R8BTO matrix:' )

  return
end
