program main

!*****************************************************************************80
!
!! MAIN is the main program for GRID_PRB.
!
!  Discussion:
!
!    GRID_PRB tests the GRID library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) center
  integer ( kind = 4 ) seed

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GRID library.'

  center = 1
  seed = 123456789

  call grid_generate_test ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat with a different seed from the first run.'

  seed = 987654321
  call grid_generate_test ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat with the same seed as the first run.'

  seed = 123456789
  call grid_generate_test ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat with different centering values.'

  do center = 1, 5

    seed = 123456789
    call grid_generate_test ( center, seed )

  end do

  call grid_side_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine grid_generate_test ( center, seed )

!*****************************************************************************80
!
!! GRID_GENERATE_TEST tests GRID_GENERATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) center
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,point_num)

  write ( *, '(a)'     ) ' '
  write ( *, '(a)'     ) 'GRID_GENERATE_TEST'
  write ( *, '(a)'     ) '  GRID_GENERATE randomly chooses a given number of'
  write ( *, '(a)'     ) '  points on a uniform grid.'
  write ( *, '(a)'     ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed
  write ( *, '(a,i8)'  ) '  Centering option =   ', center

  call grid_generate ( dim_num, point_num, center, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The grid points:'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2f10.4)' ) x(1:dim_num,j)
  end do

  return
end
subroutine grid_side_test ( )

!*****************************************************************************80
!
!! GRID_SIDE_TEST tests GRID_SIDE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log10
  integer ( kind = 4 ) n_side

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GRID_SIDE_TEST'
  write ( *, '(a)' ) '  GRID_SIDE returns the smallest N_SIDE, such that'
  write ( *, '(a)' ) '  N <= NSIDE^M'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  M      N  NSIDE  NSIDE^M'

  do m = 2, 4
    write ( *, '(a)' ) ''
    do n_log10 = 1, 4
      n = 10 ** n_log10
      call grid_side ( m, n, n_side )
      write ( *, '(2x,i1,2x,i5,2x,i5,2x,i5)' ) m, n, n_side, n_side ** m
    end do
  end do

  return
end

