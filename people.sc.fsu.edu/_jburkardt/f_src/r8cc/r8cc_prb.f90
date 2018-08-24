program main

!*****************************************************************************80
!
!! MAIN is the main program for R8CC_PRB.
!
!  Discussion:
!
!    R8CC_PRB tests the R8CC library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the R8CC library.'

  call r8cc_dif2_test ( )
  call r8cc_get_test ( )
  call r8cc_ijk_test ( )
  call r8cc_inc_test ( )
  call r8cc_indicator_test ( )
  call r8cc_kij_test ( )
  call r8cc_mtv_test ( )
  call r8cc_mv_test ( )
  call r8cc_print_test ( )
  call r8cc_print_some_test ( )
  call r8cc_random_test ( )
  call r8cc_read_test ( )
  call r8cc_set_test ( )
  call r8cc_to_r8ge_test ( )
  call r8cc_write_test ( )
  call r8cc_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8cc_dif2_test ( )

!*****************************************************************************80
!
!! R8CC_DIF2_TEST tests R8CC_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  integer ( kind = 4 ), allocatable :: colptr(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable :: rowind(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_DIF2_TEST'
  write ( *, '(a)' ) '  R8CC_DIF2 sets the second difference as an R8CC matrix;'

  m = 5
  n = 5
  if ( m == n ) then
    nz_num = 3 * m - 2
  else
    nz_num = 3 * min ( m, n ) - 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  allocate ( a(1:nz_num) )
  allocate ( rowind(1:nz_num) )
  allocate ( colptr(1:n+1) )

  call r8cc_dif2 ( m, n, nz_num, colptr, rowind, a )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The R8CC matrix:' )

  return
end
subroutine r8cc_get_test ( )

!*****************************************************************************80
!
!! R8CC_GET_TEST tests R8CC_GET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_GET_TEST'
  write ( *, '(a)' ) '  R8CC_GET gets an entry of a matrix in the R8CC format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The  R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_GET retrieves 10 entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    call r8cc_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  return
end
subroutine r8cc_ijk_test ( )

!*****************************************************************************80
!
!! R8CC_IJK_TEST tests R8CC_IJK.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_IJK_TEST'
  write ( *, '(a)' ) '  R8CC_IJK gets K from (I,J)'
  write ( *, '(a)' ) '  for a matrix in the R8CC format,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_IJK locates some (I,J) entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K'
  write ( *, '(a)' ) ' '

  do test = 1, 20
    i = i4_uniform_ab ( 1, m, seed )
    j = i4_uniform_ab ( 1, n, seed )
    call r8cc_ijk ( m, n, nz_num, colptr, rowind, i, j, k )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine r8cc_inc_test ( )

!*****************************************************************************80
!
!! R8CC_INC_TEST tests R8CC_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_INC_TEST'
  write ( *, '(a)' ) '  R8CC_INC increments entries in an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_INC increments 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K       NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 20.0D+00 + real ( test, kind = 8 )
    call r8cc_inc ( m, n, nz_num, colptr, rowind, a, i, j, value )
    call r8cc_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The final R8CC matrix:' )

  return
end
subroutine r8cc_indicator_test ( )

!*****************************************************************************80
!
!! R8CC_INDICATOR_TEST tests R8CC_INDICATOR.
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
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CC_INDICATOR sets an indicator R8CC matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8cc_indicator ( m, n, nz_num, colptr, rowind, a )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The R8CC indicator matrix:' )

  return
end
subroutine r8cc_kij_test ( )

!*****************************************************************************80
!
!! R8CC_KIJ_TEST tests R8CC_KIJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_KIJ_TEST'
  write ( *, '(a)' ) '  R8CC_KIJ gets (I,J) from K'
  write ( *, '(a)' ) '  for a matrix in the R8CC format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_KIJ locates some K entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         K         I         J'
  write ( *, '(a)' ) ' '

  do test = 1, 20  
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) k, i, j
  end do

  return
end
subroutine r8cc_mtv_test ( )

!*****************************************************************************80
!
!! R8CC_MTV_TEST tests R8CC_MTV.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_MTV_TEST'
  write ( *, '(a)' ) '  R8CC_MTV compute b=A''*x, where A is an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Compute BN = XM * A(MxN)
!
  x(1) = 2.0D+00
  x(2:m-1) = 0.0D+00
  x(m) = -3.0D+00

  call r8vec_print ( m, x, '  x:' )

  call r8cc_mtv ( m, n, nz_num, colptr, rowind, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8cc_mv_test ( )

!*****************************************************************************80
!
!! R8CC_MV_TEST tests R8CC_MV.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_MV_TEST'
  write ( *, '(a)' ) '  R8CC_MV computes b=A*x, where A is an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Compute B = A(MxN) * X
!
  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  x:' )

  call r8cc_mv ( m, n, nz_num, colptr, rowind, a, x, b )

  call r8vec_print ( m, b, '  b=A*x:' )

  return
end
subroutine r8cc_print_test ( )

!*****************************************************************************80
!
!! R8CC_PRINT_TEST tests R8CC_PRINT.
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
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_PRINT_TEST'
  write ( *, '(a)' ) '  R8CC_PRINT prints an R8CC matrix.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Print the matrix.
!
  call r8cc_print ( m, n, nz_num, colptr, rowind, a, '  The R8CC matrix:' )

  return
end
subroutine r8cc_print_some_test ( )

!*****************************************************************************80
!
!! R8CC_PRINT_SOME_TEST tests R8CC_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: nz_num = 28

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ &
    1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 29 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1,  2,  &
    1,  2,  3, &
    2,  3,  4, &
    3,  4,  5, &
    4,  5,  6, &
    5,  6,  7, &
    6,  7,  8, &
    7,  8,  9, &
    8,  9, 10, &
    9, 10 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CC_PRINT_SOME prints some of an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8cc_indicator ( m, n, nz_num, colptr, rowind, a )

  call r8cc_print_some ( m, n, nz_num, colptr, rowind, a, &
    2, 5, 6, 8, '  Rows 2-6, Cols 5-8:' )

  return
end
subroutine r8cc_random_test ( )

!*****************************************************************************80
!
!! R8CC_RANDOM_TEST tests R8CC_RANDOM.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CC_RANDOM randomizes an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  seed = 123456789
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The R8CC matrix:' )

  return
end
subroutine r8cc_read_test ( )

!*****************************************************************************80
!
!! R8CC_READ_TEST tests R8CC_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  character ( len = 127 ) :: a_file = 'r8cc_a.txt'
  integer ( kind = 4 ) base
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  character ( len = 127 ) :: col_file = 'r8cc_col.txt'
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row
  character ( len = 127 ) :: row_file = 'r8cc_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_READ_TEST'
  write ( *, '(a)' ) '  R8CC_READ reads an R8CC matrix from 3 files.'

  call r8cc_read_size ( col_file, row_file, m, n, nz_num, base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
  write ( *, '(a,i8)' ) '  Index base (0/1)  = ', base

  allocate ( a(1:nz_num) )
  allocate ( col(1:n+1) )
  allocate ( row(1:nz_num) )

  call r8cc_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  call i4vec_print ( n + 1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call r8cc_print ( m, n, nz_num, col, row, a, '  The R8CC matrix:' )

  deallocate ( a )
  deallocate ( col )
  deallocate ( row )

  return
end
subroutine r8cc_set_test ( )

!*****************************************************************************80
!
!! R8CC_SET_TEST tests R8CC_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_SET_TEST'
  write ( *, '(a)' ) '  R8CC_SET sets entries in an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_SET sets 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 100.0D+00 + real ( test, kind = 8 )
    call r8cc_set ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  return
end
subroutine r8cc_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CC_TO_R8GE_TEST tests R8CC_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a_r8cc(nz_num)
  real ( kind = 8 ) a_r8ge(m,n)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CC_TO_R8GE converts a matrix from R8CC to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8cc_indicator ( m, n, nz_num, colptr, rowind, a_r8cc )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a_r8cc, &
    '  The R8CC matrix:' )

  call r8cc_to_r8ge ( m, n, nz_num, colptr, rowind, a_r8cc, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8cc_write_test ( )

!*****************************************************************************80
!
!! R8CC_WRITE_TEST tests R8CC_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  character ( len = 127 ) :: a_file = 'r8cc_a.txt'
  integer ( kind = 4 ), dimension (n+1) :: col = (/ 1, 4, 6, 8, 10, 13 /)
  character ( len = 127 ) :: col_file = 'r8cc_col.txt'
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  character ( len = 127 ) :: row_file = 'r8cc_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_WRITE_TEST'
  write ( *, '(a)' ) '  R8CC_WRITE writes an R8CC matrix to 3 files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call r8cc_indicator ( m, n, nz_num, col, row, a )

  call r8cc_print ( m, n, nz_num, col, row, a, '  The R8CC matrix:' )

  call r8cc_write ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  return
end
subroutine r8cc_zeros_test ( )

!*****************************************************************************80
!
!! R8CC_ZEROS_TEST tests R8CC_ZEROS.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CC_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CC_ZEROS zeros an R8CC matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8cc_zeros ( m, n, nz_num, colptr, rowind, a )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The R8CC matrix:' )

  return
end

