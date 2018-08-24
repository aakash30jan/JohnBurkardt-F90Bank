program main

!*****************************************************************************80
!
!! MAIN is the main program for PARTIAL_DIGEST_PRB.
!
!  Discussion:
!
!    PARTIAL_DIGEST_PRB tests the PARTIAL_DIGEST library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTIAL_DIGEST_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PARTIAL_DIGEST library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTIAL_DIGEST_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PARTIAL_DIGEST_RECUR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nn2 = ( n * ( n - 1 ) ) / 2
!
!  Set the distance array.
!
  integer ( kind = 4 ), dimension ( ((n-1)*n)/2 ) :: dist = (/ &
    2, 2, 3, 3, 4, 5, 6, 7, 8, 10 /)
  integer ( kind = 4 ) i

  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  PARTIAL_DIGEST_RECUR generates solutions to the partial'
  write ( *, * ) '  digest problem, using recursion'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of objects to place is N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The original placement was 0,3,6,8,10.'
  write ( *, '(a)' ) '  These placements generate the following distances:'

  call i4vec_print ( nn2, dist, '  Distance array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PARTIAL_DIGEST_RECUR may recover the original placements'
  write ( *, '(a)' ) '  from the pairwise distances.  It may also find other'
  write ( *, '(a)' ) '  placements that have the same distance array.'

  call partial_digest_recur ( n, dist )

  return
end
