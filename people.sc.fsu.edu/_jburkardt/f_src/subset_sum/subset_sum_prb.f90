program main

!*****************************************************************************80
!
!! MAIN is the main program for SUBSET_SUM_PRB.
!
!  Discussion:
!
!    SUBSET_SUM_PRB tests the SUBSET_SUM library.
!
!  Licensing:
!
!    I don't care what you do with this code.
!
!  Modified:
!
!    15 July 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBSET_SUM_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SUBSET_SUM library.'

  call subset_sum_count_tests ( )
  call subset_sum_find_tests ( )
  call subset_sum_next_tests ( )
  call subset_sum_table_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBSET_SUM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end

