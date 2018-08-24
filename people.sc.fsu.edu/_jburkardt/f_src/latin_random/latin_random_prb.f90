program main

!*****************************************************************************80
!
!! MAIN is the main program for LATIN_RANDOM_PRB.
!
!  Discussion:
!
!    LATIN_RANDOM_PRB tests the LATIN_RANDOM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_RANDOM_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATIN_RANDOM library.'

  seed = 123456789

  do test = 1, 3

    call test01 ( seed )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_RANDOM_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( seed )

!*****************************************************************************80
!
!! TEST01 tests LATIN_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LATIN_RANDOM chooses a random Latin Square'
  write ( *, '(a)' ) '  cell arrangement, and then returns'
  write ( *, '(a)' ) '  a random point from each cell.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i6)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed

  call latin_random ( dim_num, point_num, seed, x )

  call r8mat_transpose_print ( dim_num, point_num, x, &
    '  The Latin Random points:' )

  return
end
