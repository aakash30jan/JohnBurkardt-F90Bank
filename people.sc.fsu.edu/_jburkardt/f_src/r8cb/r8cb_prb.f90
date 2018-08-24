program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPLUS_PRB.
!
!  Discussion:
!
!    LINPLUS_PRB tests the LINPLUS library.
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
  write ( *, '(a)' ) 'LINPLUS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the LINPLUS library.'

  call r8cb_det_test ( )
  call r8cb_dif2_test ( )
  call r8cb_indicator_test ( )
  call r8cb_ml_test ( )
  call r8cb_mtv_test ( )
  call r8cb_mv_test ( )
  call r8cb_np_fa_test ( )
  call r8cb_np_sl_test ( )
  call r8cb_print_test ( )
  call r8cb_print_some_test ( )
  call r8cb_random_test ( )
  call r8cb_to_r8ge_test ( )
  call r8cb_to_r8vec_test ( )
  call r8cb_zeros_test ( )
  call r8vec_to_r8cb_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPLUS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8cb_det_test ( )

!*****************************************************************************80
!
!! R8CB_DET_TEST tests R8CB_DET.
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
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) a_r8ge(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_DET_TEST'
  write ( *, '(a)' ) '  R8CB_DET computes the determinant of a matrix factored'
  write ( *, '(a)' ) '  by R8CB_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cb_random ( n, n, ml, mu, seed, a )

  call r8cb_print ( n, n, ml, mu, a, '  The compact band matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8cb_to_r8ge ( n, n, ml, mu, a, a_r8ge )
!
!  Factor the matrix.
!
  call r8cb_np_fa ( n, ml, mu, a, info )
!
!  Compute the determinant.
!
  call r8cb_det ( n, ml, mu, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8CB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call r8ge_fa ( n, a_r8ge, pivot, info )
!
!  Compute the determinant.
!
  call r8ge_det ( n, a_r8ge, pivot, det )

  write ( *, '(a,g14.6)' ) '  R8GE_DET computes the determinant = ', det

  return
end
subroutine r8cb_dif2_test ( )

!*****************************************************************************80
!
!! R8CB_DIF2_TEST tests R8CB_DIF2.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_DIF2_TEST'
  write ( *, '(a)' ) '  R8CB_DIF2 computes an R8CB second difference matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_dif2 ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB second difference matrix:' )

  return
end
subroutine r8cb_indicator_test ( )

!*****************************************************************************80
!
!! R8CB_INDICATOR_TEST tests R8CB_INDICATOR.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CB_INDICATOR computes an R8CB indicator matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  return
end
subroutine r8cb_ml_test ( )

!*****************************************************************************80
!
!! R8CB_ML_TEST tests R8CB_ML.
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
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) b2(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_ML_TEST'
  write ( *, '(a)' ) '  R8CB_ML computes A*x or A''*x'
  write ( *, '(a)' ) '  for a compact band matrix A which has'
  write ( *, '(a)' ) '  been factored by R8CB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8CB_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8cb_ml ( n, ml, mu, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print ( n, b, b2, '  A*x and PLU*x' )
    else
      call r8vec2_print ( n, b, b2, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8cb_mtv_test ( )

!*****************************************************************************80
!
!! R8CB_MTV_TEST tests R8CB_MTV.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_MTV_TEST'
  write ( *, '(a)' ) '  R8CB_MTV computes b=A''*x, where A is an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r8cb_mtv ( m, n, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8cb_mv_test ( )

!*****************************************************************************80
!
!! R8CB_MV_TEST tests R8CB_MV.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_MV_TEST'
  write ( *, '(a)' ) '  R8CB_MV computes b=A*x, where A is an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8cb_mv ( m, n, ml, mu, a, x, b )

  call r8vec_print ( m, b, '  The product b=A*x:' )

  return
end
subroutine r8cb_np_fa_test ( )

!*****************************************************************************80
!
!! R8CB_NP_FA_TEST tests R8CB_NP_FA.
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
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_NP_FA_TEST'
  write ( *, '(a)' ) '  R8CB_NP_FA factors an R8CB matrix, no pivoting:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_FA_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8CB_NP_FA claims the matrix is singular.'
      write ( *, '(a,i8)' ) '  The value of info is ', info
      return
    end if
!
!  Solve the system.
!
    call r8cb_np_sl ( n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine r8cb_np_sl_test ( )

!*****************************************************************************80
!
!! R8CB_NP_SL_TEST tests R8CB_NP_SL.
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
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_NP_SL_TEST'
  write ( *, '(a)' ) '  R8CB_NP_SL solves a linear system factored by R8CB_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_SL_TEST Fatal error!'
      write ( *, '(a)' ) '  R8CB_NP_FA claims the matrix is singular.'
      write ( *, '(a,i8)' ) '  The value of info is ', info
      return
    end if
!
!  Solve the system.
!
    call r8cb_np_sl ( n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine r8cb_print_test ( )

!*****************************************************************************80
!
!! R8CB_PRINT_TEST tests R8CB_PRINT.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_PRINT_TEST'
  write ( *, '(a)' ) '  R8CB_PRINT prints an R8CB matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  return
end
subroutine r8cb_print_some_test ( )

!*****************************************************************************80
!
!! R8CB_PRINT_SOME_TEST tests R8CB_PRINT_SOME.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CB_PRINT_SOME prints some of an R8CB matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print_some ( m, n, ml, mu, a, 3, 3, 6, 6, '  Rows 3-6, Cols 3-6:' )

  return
end
subroutine r8cb_random_test ( )

!*****************************************************************************80
!
!! R8CB_RANDOM_TEST tests R8CB_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CB_RANDOM randomizes an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  seed = 123456789
  call r8cb_random ( m, n, ml, mu, seed, a )

  call r8cb_print ( m, n, ml, mu, a, '  The random R8CB matrix:' )

  return
end
subroutine r8cb_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CB_TO_R8GE_TEST tests R8CB_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) a_r8ge(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CB_TO_R8GE converts an R8CB matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8cb_to_r8ge ( m, n, ml, mu, a, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8cb_to_r8vec_test ( )

!*****************************************************************************80
!
!! R8CB_TO_R8VEC_TEST tests R8CB_TO_R8VEC.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x((ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_TO_R8VEC_TEST'
  write ( *, '(a)' ) '  R8CB_TO_R8VEC converts an R8CB matrix to an R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8cb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, ml + mu + 1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8cb ( m, n, ml, mu, x, a )

  call r8cb_print ( m, n, ml, mu, a, '  The recovered R8CB matrix:' )

  return
end
subroutine r8cb_zeros_test ( )

!*****************************************************************************80
!
!! R8CB_ZEROS_TEST tests R8CB_ZEROS.
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

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CB_ZEROS zeros an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_zeros ( n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB zero matrix:' )

  return
end
subroutine r8vec_to_r8cb_test ( )

!*****************************************************************************80
!
!! R8VEC_TO_R8CB_TEST tests R8VEC_TO_R8CB.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x((ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_TO_R8CB_TEST'
  write ( *, '(a)' ) '  R8VEC_TO_R8CB converts an R8VEC to an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  call r8cb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, ml + mu + 1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8cb ( m, n, ml, mu, x, a )

  call r8cb_print ( m, n, ml, mu, a, '  The recovered R8CB indicator matrix:' )

  return
end

