program main

!*****************************************************************************80
!
!! MAIN is the main program for EP_SERIAL.
!
!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                      S E R I A L     V E R S I O N                      !
!                                                                         !
!                                   E P                                   !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial version of the NPB EP code.               !
!    Refer to NAS Technical Reports 95-020 for details.                   !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.3. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.3, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

!  Discussion:
!
!    This is the serial version of the APP Benchmark 1,
!    the "embarassingly parallel" benchmark.
!
!    M is the Log_2 of the number of complex pairs of uniform (0, 1) random
!    numbers.  MK is the Log_2 of the size of each batch of uniform random
!    numbers.  MK can be set for convenience on a given system, since it does
!    not affect the results.
!
!  Author:
!
!    P. O. Frederickson, D. H. Bailey, A. C. Woo
!
  implicit none

  include 'npbparams.h'

  double precision Mops, epsilon, a, s, t1, t2, t3, t4, x, x1, &
                   x2, q, sx, sy, tm, an, tt, gc, dum(3)
  double precision sx_verify_value, sy_verify_value, sx_err, sy_err
  integer          mk, mm, nn, nk, nq, np, &
                   i, ik, kk, l, k, nit, &
                   k_offset, j, fstatus
  logical          verified, timers_enabled
  external         randlc, timer_read
  double precision randlc, timer_read
  character*15     size

  parameter (mk = 16, mm = m - mk, nn = 2 ** mm, &
             nk = 2 ** mk, nq = 10, epsilon=1.d-8, &
             a = 1220703125.d0, s = 271828183.d0)

  common/storage/ x(2*nk), q(0:nq-1)
  data             dum /1.d0, 1.d0, 1.d0/

  call timestamp ( )

  open(unit=2, file='timer.flag', status='old', iostat=fstatus)
  if (fstatus == 0) then
     timers_enabled = .true.
     close(2)
  else
     timers_enabled = .false.
  endif
!
!   Because the size of the problem is too large to store in a 32-bit
!   integer for some classes, we put it into a string (for printing).
!   Have to strip off the decimal point put in there by the floating
!   point print statement (internal file)
!
  write(*, 1000)
  write(size, '(f15.0)' ) 2.d0**(m+1)
  j = 15
  if (size(j:j) == '.') j = j - 1
  write (*,1001) size(1:j)
  write (*,*)

 1000 format(//,' NAS Parallel Benchmarks (NPB3.3-SER)', &
            ' - EP Benchmark', /)
 1001 format(' Number of random numbers generated: ', a15)

  verified = .false.
!
!   Compute the number of "batches" of random number pairs generated 
!   per processor. Adjust if the number of processors does not evenly 
!   divide the total number
!
  np = nn 
!
!   Call the random number generator functions and initialize
!   the x-array to reduce the effects of paging on the timings.
!   Also, call all mathematical functions that are used. Make
!   sure these initializations cannot be eliminated as dead code.
!
  call vranlc(0, dum(1), dum(2), dum(3))
  dum(1) = randlc(dum(2), dum(3))

  do i = 1, 2*nk
     x(i) = -1.d99
  end do

  Mops = log(sqrt(abs(max(1.d0,1.d0))))

  call timer_clear(1)
  call timer_clear(2)
  call timer_clear(3)
  call timer_start(1)

  t1 = a
  call vranlc(0, t1, a, x)
!
!   Compute AN = A ^ (2 * NK) (mod 2^46).
!
  t1 = a

  do i = 1, mk + 1
     t2 = randlc(t1, t1)
  end do

  an = t1
  tt = s
  gc = 0.d0
  sx = 0.d0
  sy = 0.d0

  do i = 0, nq - 1
     q(i) = 0.d0
  end do
!
!   Each instance of this loop may be performed independently. We compute
!   the k offsets separately to take into account the fact that some nodes
!   have more numbers to generate than others
!
  k_offset = -1

  do k = 1, np

     kk = k_offset + k 
     t1 = s
     t2 = an
!
!  Find starting seed t1 for this kk.
!
     do i = 1, 100
        ik = kk / 2
        if (2 * ik .ne. kk) t3 = randlc(t1, t2)
        if (ik == 0) goto 130
        t3 = randlc(t2, t2)
        kk = ik
     end do
!
!  Compute uniform pseudorandom numbers.
!
 130     continue

     if (timers_enabled) call timer_start(3)
     call vranlc(2 * nk, t1, a, x)
     if (timers_enabled) call timer_stop(3)
!
!  Compute Gaussian deviates by acceptance-rejection method and 
!  tally counts in concentric square annuli.  This loop is not 
!  vectorizable. 
!
     if (timers_enabled) call timer_start(2)

     do i = 1, nk
        x1 = 2.d0 * x(2*i-1) - 1.d0
        x2 = 2.d0 * x(2*i) - 1.d0
        t1 = x1 ** 2 + x2 ** 2
        if (t1 .le. 1.d0) then
           t2   = sqrt(-2.d0 * log(t1) / t1)
           t3   = (x1 * t2)
           t4   = (x2 * t2)
           l    = max(abs(t3), abs(t4))
           q(l) = q(l) + 1.d0
           sx   = sx + t3
           sy   = sy + t4
        endif
     end do

     if (timers_enabled) call timer_stop(2)

  end do

  do i = 0, nq - 1
    gc = gc + q(i)
  end do

  call timer_stop(1)
  tm  = timer_read(1)

  nit=0
  verified = .true.
  if (m==24) then
     sx_verify_value = -3.247834652034740D+3
     sy_verify_value = -6.958407078382297D+3
  elseif (m==25) then
     sx_verify_value = -2.863319731645753D+3
     sy_verify_value = -6.320053679109499D+3
  elseif (m==28) then
     sx_verify_value = -4.295875165629892D+3
     sy_verify_value = -1.580732573678431D+4
  elseif (m==30) then
     sx_verify_value =  4.033815542441498D+4
     sy_verify_value = -2.660669192809235D+4
  elseif (m==32) then
     sx_verify_value =  4.764367927995374D+4
     sy_verify_value = -8.084072988043731D+4
  elseif (m==36) then
     sx_verify_value =  1.982481200946593D+5
     sy_verify_value = -1.020596636361769D+5
  elseif (m==40) then
     sx_verify_value = -5.319717441530D+05
     sy_verify_value = -3.688834557731D+05
  else
     verified = .false.
  endif

  if (verified) then
     sx_err = abs((sx - sx_verify_value)/sx_verify_value)
     sy_err = abs((sy - sy_verify_value)/sy_verify_value)
     verified = ((sx_err.le.epsilon) .and. (sy_err.le.epsilon))
  endif
  Mops = 2.d0**(m+1)/tm/1000000.d0

  write (6,11) tm, m, gc, sx, sy, (i, q(i), i = 0, nq - 1)
 11   format ('EP Benchmark Results:'//'CPU Time =',f10.4/'N = 2^', &
          i5/'No. Gaussian Pairs =',f15.0/'Sums = ',1p,2d25.15/ &
          'Counts:'/(i3,0p,f15.0))

  call print_results('EP', class, m+1, 0, 0, nit, &
                     tm, Mops,  &
                     'Random numbers generated',  &
                     verified, npbversion, compiletime, cs1, &
                     cs2, cs3, cs4, cs5, cs6, cs7)

  if (timers_enabled) then
     if (tm .le. 0.d0) tm = 1.0
     tt = timer_read(1)
     print 810, 'Total time:    ', tt, tt*100./tm
     tt = timer_read(2)
     print 810, 'Gaussian pairs:', tt, tt*100./tm
     tt = timer_read(3)
     print 810, 'Random numbers:', tt, tt*100./tm
810      format(1x,a,f9.3,' (',f6.2,'%)')
  endif
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EP_SERIAL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

 stop
end
subroutine print_results ( name, class, n1, n2, n3, niter,  &
                t, mops, optype, verified, npbversion,  &
                compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)

!*****************************************************************************80
!
!! PRINT_RESULTS prints the results.
!
  implicit none

  character name*(*)
  character class*1
  integer   n1, n2, n3, niter, j
  double precision t, mops
  character optype*24, size*15
  logical   verified
  character*(*) npbversion, compiletime,  &
                cs1, cs2, cs3, cs4, cs5, cs6, cs7

     write (*, 2) name 
 2       format(//, ' ', A, ' Benchmark Completed.')

     write (*, 3) Class
 3       format(' Class           = ', 12x, a12)
!
!   If this is not a grid-based problem (EP, FT, CG), then
!   we only print n1, which contains some measure of the
!   problem size. In that case, n2 and n3 are both zero.
!   Otherwise, we print the grid size n1xn2xn3
!
     if ((n2 == 0) .and. (n3 == 0)) then
        if (name(1:2) == 'EP') then
           write(size, '(f15.0)' ) 2.d0**n1
           j = 15
           if (size(j:j) == '.') then
              size(j:j) = ' '
              j = j - 1
           endif
           write (*,42) size(1:j)
 42            format(' Size            = ',9x, a15)
        else
           write (*,44) n1
 44            format(' Size            = ',12x, i12)
        endif
     else
        write (*, 4) n1,n2,n3
 4          format(' Size            =  ',9x, i4,'x',i4,'x',i4)
     endif

     write (*, 5) niter
 5       format(' Iterations      = ', 12x, i12)
     
     write (*, 6) t
 6       format(' Time in seconds = ',12x, f12.2)
     
     write (*,9) mops
 9       format(' Mop/s total     = ',12x, f12.2)

     write(*, 11) optype
 11      format(' Operation type  = ', a24)

     if (verified) then 
        write(*,12) '  SUCCESSFUL'
     else
        write(*,12) 'UNSUCCESSFUL'
     endif
 12      format(' Verification    = ', 12x, a)

     write(*,13) npbversion
 13      format(' Version         = ', 12x, a12)

     write(*,14) compiletime
 14      format(' Compile date    = ', 12x, a12)


     write (*,121) cs1
 121     format(/, ' Compile options:', /,  &
            '    F77          = ', A)

     write (*,122) cs2
 122     format('    FLINK        = ', A)

     write (*,123) cs3
 123     format('    F_LIB        = ', A)

     write (*,124) cs4
 124     format('    F_INC        = ', A)

     write (*,125) cs5
 125     format('    FFLAGS       = ', A)

     write (*,126) cs6
 126     format('    FLINKFLAGS   = ', A)

     write(*, 127) cs7
 127     format('    RAND         = ', A)
    
     write (*,130)
 130     format(//' Please send all errors/feedbacks to:'// &
              ' NPB Development Team'/ &
              ' npb@nas.nasa.gov'//)

  return
end
function randlc ( x, a )

!*****************************************************************************80
!
!! RANDLC returns a uniform pseudorandom double precision number.
!
!  Discussion:
!
!    The number returned is in the range (0, 1).  
!
!    The algorithm uses the linear congruential generator
!
!   x_{k+1} = a x_k  (mod 2^46)
!
!   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
!   before repeating.  The argument A is the same as 'a' in the above formula,
!   and X is the same as x_0.  A and X must be odd double precision integers
!   in the range (1, 2^46).  The returned value RANDLC is normalized to be
!   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
!   the new seed x_1, so that subsequent calls to RANDLC using the same
!   arguments will generate a continuous sequence.
!
  implicit none

  double precision x, a
  integer*8 i246m1, Lx, La
  double precision d2m46
  double precision randlc

  parameter(d2m46=0.5d0**46)

  save i246m1
  data i246m1/X'00003FFFFFFFFFFF'/

  Lx = X
  La = A

  Lx   = iand(Lx*La,i246m1)
  randlc = d2m46*dble(Lx)
  x    = dble(Lx)

  return
end
subroutine vranlc ( n, x, a, y )

!*****************************************************************************80
!
!! VRANLC returns a vector of uniform pseudorandom double precision numbers.
!
  implicit none

  integer n

  double precision a
  double precision, parameter :: d2m46 = 0.5D+00 ** 46
  integer i
  integer*8 i246m1
  integer*8 la
  integer*8 lx
  double precision x
  double precision y(n)

  save i246m1

  data i246m1 / X'00003FFFFFFFFFFF' /

  lx = x
  la = a
  do i = 1, n
    lx = iand ( lx * la, i246m1 )
    y(i) = d2m46 * dble ( lx )
  end do

  x = dble ( lx )

  return
end      
subroutine timer_clear ( n )

!*****************************************************************************80
!
!! TIMER_CLEAR clears the timer.
!
  implicit none

  integer n
  
  double precision start(64), elapsed(64)
  common /tt/ start, elapsed

  elapsed(n) = 0.0
  return
end
subroutine timer_start ( n )

!*****************************************************************************80
!
!! TIMER_START starts the timer.
!
  implicit none

  external         elapsed_time
  double precision elapsed_time
  integer n
  double precision start(64), elapsed(64)
  common /tt/ start, elapsed

  start(n) = elapsed_time()

  return
end
subroutine timer_stop ( n )

!*****************************************************************************80
!
!! TIMER_STOP stops the timer.
!
  implicit none

  external         elapsed_time
  double precision elapsed_time
  integer n
  double precision start(64), elapsed(64)
  common /tt/ start, elapsed
  double precision t, now
  now = elapsed_time()
  t = now - start(n)
  elapsed(n) = elapsed(n) + t

  return
end
function timer_read ( n )

!*****************************************************************************80
!
!! TIMER_READ reads the timer.
!
  implicit none

  integer n
  double precision start(64), elapsed(64)
  double precision timer_read

  common /tt/ start, elapsed
  
  timer_read = elapsed(n)
  return
end
function elapsed_time ( )

!*****************************************************************************80
!
!! ELAPSED_TIME measures wall clock time.
!
  implicit none

  double precision elapsed_time
  double precision t
!
!  This function must measure wall clock time, not CPU time. 
!  Since there is no portable timer in Fortran (77)
!  we call a routine compiled in C (though the C source may have
!  to be tweaked). 
!
  call wtime(t)
!
!  The following is not ok for "official" results because it reports
!  CPU time not wall clock time. It may be useful for developing/testing
!  on timeshared Crays, though. 
!     call second(t)
!
  elapsed_time = t

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
