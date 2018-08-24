program main

!*******************************************************************************
!
!! MAIN is the main program for TOMS659_PRB.
!
!  Discussion:
!
!    TOMS659_PRB tests the TOMS659 library.
!
!  Modified:
!
!    06 October 2013
!
  integer, parameter :: maxdim = 1111

  integer atmost
  integer dimen
  double precision f
  logical flag(2)
  integer i
  integer j
  double precision quasi(maxdim)
  double precision sum
  integer taus

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS659_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS659 library.'

  do

    read (*,*) dimen,atmost

    if (dimen.eq.0) then
      exit
    end if

    write ( *, * ) ''
    write ( *, * ) 'Test SOBOL'
    write ( *, * ) 'DIMENSION = ',DIMEN
    write ( *, * ) 'ATMOST = ',ATMOST

    call cpu_time ( t1 )

    call insobl ( dimen, atmost )

    write ( *,*) 'Start time = ',T1
    write ( *,*) 'I = Iteration number'
    write ( *,*) 'EI = Estimate of integral'
    write ( *,*) ''

    sum = 0.0D+00
    do i = 1, atmost

      call i4_sobol ( dimen, quasi )

      f = 1.0
      do j = 1, dimen
        f = f * abs ( 4.0D+00 * quasi(j) - 2.0D+00 )
      end do
      sum = sum + f

      if ( mod ( i, 5000 ) .eq. 0 ) then
        write ( *,*) 'i = ', i
        write ( *,*) 'ei = ',sum / i
        call cpu_time ( t2 )
        write ( *,*) 'TIMEI = ', t2 - t1
        write ( *,'(1h )')
      end if

    end do

    write ( *, * ) ' EI = ',sum / atmost

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS659_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
  end
