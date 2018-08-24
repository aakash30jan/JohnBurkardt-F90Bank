program main

!*****************************************************************************80
!
!! MAIN is the main program for R8CI_PRB.
!
!  Discussion:
!
!    R8CI_PRB tests the R8CI library.
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
  write ( *, '(a)' ) 'R8CI_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8CI library.'

  call r8ci_det_test ( )
  call r8ci_dif2_test ( )
  call r8ci_eval_test ( )
  call r8ci_indicator_test ( )
  call r8ci_mtv_test ( )
  call r8ci_mv_test ( )
  call r8ci_print_test ( )
  call r8ci_print_some_test ( )
  call r8ci_random_test ( )
  call r8ci_sl_test ( )
  call r8ci_to_r8ge_test ( )
  call r8ci_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8ci_det_test ( )

!*****************************************************************************80
!
!! R8CI_DET_TEST tests R8CI_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_DET_TEST'
  write ( *, '(a)' ) '  R8CI_DET finds the determinant of '
  write ( *, '(a)' ) '  a real circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  call r8ci_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Computed determinant = ', det

  return
end
subroutine r8ci_dif2_test ( )

!*****************************************************************************80
!
!! R8CI_DIF2_TEST tests R8CI_DIF2.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_DIF2_TEST'
  write ( *, '(a)' ) '  R8CI_DIF2 sets up an R8CI periodic second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_dif2 ( n, a )

  call r8ci_print ( n, a, '  The R8CI second difference matrix:' )

  return
end
subroutine r8ci_eval_test ( )

!*****************************************************************************80
!
!! R8CI_EVAL_TEST tests R8CI_EVAL.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  complex ( kind = 8 ) lambda(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_EVAL_TEST'
  write ( *, '(a)' ) '  R8CI_EVAL finds the eigenvalues of '
  write ( *, '(a)' ) '  a real circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  call r8ci_eval ( n, a, lambda )

  call c8vec_print ( n, lambda, '  The eigenvalues:' )

  return
end
subroutine r8ci_indicator_test ( )

!*****************************************************************************80
!
!! R8CI_INDICATOR_TEST tests R8CI_INDICATOR.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CI_INDICATOR sets up an R8CI indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_mtv_test ( )

!*****************************************************************************80
!
!! R8CI_MTV_TEST tests R8CI_MTV.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_MTV_TEST'
  write ( *, '(a)' ) '  R8CI_MTV computes b=A''*x, where A is an R8CI matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )
  call r8ci_print ( n, a, '  The circulant matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector x:' );

  call r8ci_mtv ( n, a, x, b )
  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8ci_mv_test ( )

!*****************************************************************************80
!
!! R8CI_MV_TEST tests R8CI_MV.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_MV_TEST'
  write ( *, '(a)' ) '  R8CI_MV computes b=A*x, where A is an R8CI matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )
  call r8ci_print ( n, a, '  The circulant matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector x:' );

  call r8ci_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  The product b=A*x:' )

  return
end
subroutine r8ci_print_test ( )

!*****************************************************************************80
!
!! R8CI_PRINT_TEST tests R8CI_PRINT.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_PRINT_TEST'
  write ( *, '(a)' ) '  R8CI_PRINT prints an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_print_some_test ( )

!*****************************************************************************80
!
!! R8CI_PRINT_SOME_TEST tests R8CI_PRINT_SOME.
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

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CI_PRINT_SOME prints some of an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print_some ( n, a, 2, 3, 6, 5, '  Rows 2-6, Cols 3-5:' )

  return
end
subroutine r8ci_random_test ( )

!*****************************************************************************80
!
!! R8CI_RANDOM_TEST tests R8CI_RANDOM.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CI_RANDOM sets a random R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_sl_test ( )

!*****************************************************************************80
!
!! R8CI_SL_TEST tests R8CI_SL.
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

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_SL_TEST'
  write ( *, '(a)' ) '  R8CI_SL solves a circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8ci_mv ( n, a, x, b )
    else
      call r8ci_mtv ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call r8ci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine r8ci_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CI_TO_R8GE_TEST tests R8CI_TO_R8GE.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CI_TO_R8GE converts an R8CI matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  call r8ci_to_r8ge ( n, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8ci_zeros_test ( )

!*****************************************************************************80
!
!! R8CI_ZEROS_TEST tests R8CI_ZEROS.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CI_ZEROS zeros an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_zeros ( n, a )

  call r8ci_print ( n, a, '  The zero R8CI matrix:' )

  return
end

