program main

!*****************************************************************************80
!
!! MAIN is the main program for DUCK_PRB.
!
  use duck

  integer boop
  integer toop
  integer grink

  call flink ( toop )
  write ( *, * ) 'Toop = ', toop
  boop = grink ( 11 )
  write ( *, * ) 'Boop = ', boop

  stop
end
