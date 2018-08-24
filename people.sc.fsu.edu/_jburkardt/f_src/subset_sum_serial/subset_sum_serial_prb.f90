program main

!*****************************************************************************80
!
!! MAIN is the main program for SUBSET_SUM_SERIAL_PRB.
!
!  Discussion:
!
!    SUBSET_SUM_SERIAL_PRB tests the SUBSET_SUM_SERIAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SUBSET_SUM_SERIAL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the SUBSET_SUM_SERIAL library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SUBSET_SUM_SERIAL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! SUBSET_SUM_SERIAL_TEST01 tests the SUBSET_SUM_SERIAL program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21

  integer ( kind = 4 ) choice(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) target
  integer ( kind = 4 ), dimension ( n ) :: weight = (/ &
    518533, 1037066, 2074132, 1648264, 796528, &
   1593056,  686112, 1372224,  244448, 488896, &
    977792, 1955584, 1411168,  322336, 644672, &
   1289344,   78688,  157376,  314752, 629504, &
   1259008 /)
  integer ( kind = 4 ) w_sum

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test the SUBSET_SUM_SERIAL function, which looks for a selection'
  write ( *, '(a)' ) '  from a set of weights that adds up to a given target.'
!
!  Define the problem data.
!
  target = 2463098
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Target value:'
  write ( *, '(2x,i8)' ) target
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I      W(I)'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i2,2x,i8)' ) i, weight(i)
  end do

  call subset_sum_serial ( n, weight, target, choice )

  if ( choice(1) == -1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '  No solution was found.'
  else
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '   I*     W*'
    write ( *, '(a)' ) ''
    w_sum = 0
    do i = 1, n
      if ( choice(i) == 1 ) then
        w_sum = w_sum + weight(i)
        write ( *, '(2x,i2,2x,i8)' ) i, weight(i)
      end if
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a,i12)' ) '  Sum:    ', w_sum
    write ( *, '(a,i12)' ) '  Target: ', target
  end if

  return
end

