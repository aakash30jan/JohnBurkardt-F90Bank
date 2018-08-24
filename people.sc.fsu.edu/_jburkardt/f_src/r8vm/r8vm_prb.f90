program main

!*****************************************************************************80
!
!! MAIN is the main program for R8VM_PRB.
!
!  Discussion:
!
!    R8VM_PRB tests the R8VM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8VM library.'

  call r8ge_to_r8vm_test ( )

  call r8vm_det_test ( )
  call r8vm_indicator_test ( )
  call r8vm_mtv_test ( )
  call r8vm_mv_test ( )
  call r8vm_print_test ( )
  call r8vm_print_some_test ( )
  call r8vm_random_test ( )
  call r8vm_sl_test ( )
  call r8vm_slt_test ( )
  call r8vm_to_r8ge_test ( )
  call r8vm_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8ge_to_r8vm_test ( )

!*****************************************************************************80
!
!! R8GE_TO_R8VM_TEST tests R8GE_TO_R8VM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a_ge(m,n)
  real ( kind = 8 ) a_vm(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_TO_R8VM_TEST'
  write ( *, '(a)' ) '  R8GE_TO_R8VM converts an R8GE matrix to R8VM format.'

  seed = 123456789
  call r8ge_random ( m, n, seed, a_ge )

  call r8ge_print ( m, n, a_ge, '  The random R8GE matrix:' )

  call r8ge_to_r8vm ( m, n, a_ge, a_vm )

  call r8vm_print ( m, n, a_vm, '  The R8VM matrix' )

  return
end
subroutine r8vm_det_test ( )

!*****************************************************************************80
!
!! R8VM_DET_TEST tests R8VM_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) det1
  real ( kind = 8 ) det2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_DET_TEST'
  write ( *, '(a)' ) '  R8VM_DET, determinant of a Vandermonde matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8vm_indicator ( n, n, a )

  call r8vm_print ( n, n, a, '  The Vandermonde matrix:' )
!
!  Compute the determinant.
!
  call r8vm_det ( n, a, det1 )
!
!  Compare to exact value.
!
  call r8vm_indicator_det ( n, det2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8VM_DET    = ', det1
  write ( *, '(a,g14.6)' ) '  Exact value = ', det2

  return
end
subroutine r8vm_indicator_test ( )

!*****************************************************************************80
!
!! R8VM_INDICATOR_TEST tests R8VM_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8VM_INDICATOR sets an indicator matrix in R8VM format.'

  call r8vm_indicator ( m, n, a )

  call r8vm_print ( m, n, a, '  The R8VM indicator matrix:' )
 
  return
end
subroutine r8vm_mtv_test ( )

!*****************************************************************************80
!
!! R8VM_MTV_TEST tests R8VM_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_MTV_TEST'
  write ( *, '(a)' ) '  R8VM_MTV computes b=A''*x, where A is an R8VM matrix.'
!
!  Set A.
!
  seed = 123456789
  call r8vm_random ( n, n, seed, a )
!
!  Set X.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute b=A'*x.
!
  call r8vm_mtv ( n, n, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x' )
 
  return
end
subroutine r8vm_mv_test ( )

!*****************************************************************************80
!
!! R8VM_MV_TEST tests R8VM_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_MV_TEST'
  write ( *, '(a)' ) '  R8VM_MV computes b=A*x, where A is an R8VM matrix.'
!
!  Set A.
!
  seed = 123456789
  call r8vm_random ( n, n, seed, a )
!
!  Set X.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute b=A*x.
!
  call r8vm_mv ( n, n, a, x, b )

  call r8vec_print ( n, b, '  b=A*x' )
 
  return
end
subroutine r8vm_print_test ( )

!*****************************************************************************80
!
!! R8VM_PRINT_TEST tests R8VM_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_PRINT_TEST'
  write ( *, '(a)' ) '  R8VM_PRINT prints an R8VM format.'

  call r8vm_indicator ( m, n, a )

  call r8vm_print ( m, n, a, '  The R8VM matrix:' )
 
  return
end
subroutine r8vm_print_some_test ( )

!*****************************************************************************80
!
!! R8VM_PRINT_SOME_TEST tests R8VM_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8VM_PRINTSOME prints some of an R8VM format.'

  call r8vm_indicator ( m, n, a )

  call r8vm_print_some ( m, n, a, 2, 2, 5, 4, '  Rows 2-5, Cols 2:4' )
 
  return
end
subroutine r8vm_random_test ( )

!*****************************************************************************80
!
!! R8VM_RANDOM_TEST tests R8VM_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_RANDOM_TEST'
  write ( *, '(a)' ) '  R8VM_RANDOM randomizes an R8VM matrix.'

  seed = 123456789
  call r8vm_random ( m, n, seed, a )

  call r8vm_print ( m, n, a, '  The R8VM matrix:' )

  return
end
subroutine r8vm_sl_test ( )

!*****************************************************************************80
!
!! R8VM_SL_TEST tests R8VM_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_SL_TEST'
  write ( *, '(a)' ) '  R8VM_SL solves a Vandermonde system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8vm_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8vm_mv ( n, n, a, x, b )
!
!  Solve the linear system.
!
  call r8vm_sl ( n, a, b, x, info )

  call r8vec_print_some ( n, x, 10, '  Solution:' )
 
  return
end
subroutine r8vm_slt_test ( )

!*****************************************************************************80
!
!! R8VM_SLT_TEST tests R8VM_SLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_SLT_TEST'
  write ( *, '(a)' ) '  R8VM_SLT solves a transposed Vandermonde system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8vm_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8vm_mtv ( n, n, a, x, b )
!
!  Solve the linear system.
!
  call r8vm_slt ( n, a, b, x, info )

  call r8vec_print_some ( n, x, 10, '  Solution to transposed system:' )

  return
end
subroutine r8vm_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8VM_TO_R8GE_TEST tests R8VM_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a_ge(m,n)
  real ( kind = 8 ) a_vm(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8VM_TO_R8GE converts an R8VM matrix to R8GE format.'

  call r8vm_random ( m, n, seed, a_vm )

  call r8vm_print ( m, n, a_vm, '  The random R8VM matrix:' )

  call r8vm_to_r8ge ( m, n, a_vm, a_ge )

  call r8ge_print ( m, n, a_ge, '  The R8GE matrix' )

  return
end
subroutine r8vm_zeros_test ( )

!*****************************************************************************80
!
!! R8VM_ZEROS_TEST tests R8VM_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VM_ZEROS_TEST'
  write ( *, '(a)' ) '  R8VM_ZEROS zeros out a matrix in R8VM format.'

  call r8vm_zeros ( m, n, a )

  call r8vm_print ( m, n, a, '  The zero R8VM matrix:' )
 
  return
end
