program main

!*****************************************************************************80
!
!! MAIN is the main program for PCHIP_PRB.
!
!  Discussion:
!
!    PCHIP_PRB tests the PCHIP library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ipass
  integer ( kind = 4 ) kprint
  integer ( kind = 4 ) lun

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCHIP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PCHIP library.'

  call dpchev_test ( )
  call dpchev_test02 ( )
  call dpchev_test03 ( )
  call dpchqa_test ( )
!
!  Call the built-in quick checks.
!
  lun = 6
  kprint = 1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Calling built-in check routines:'
  write ( *, '(a)' ) ''

  call dpchq1 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  DPCHQ1 passed.'
  else
    write ( *, '(a)' ) '  DPCHQ1 failed.'
  end if

  call dpchq2 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  DPCHQ2 passed.'
  else
    write ( *, '(a)' ) '  DPCHQ2 failed.'
  end if

  call dpchq3 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  DPCHQ3 passed.'
  else
    write ( *, '(a)' ) '  DPCHQ3 failed.'
  end if

  call dpchq4 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  DPCHQ4 passed.'
  else
    write ( *, '(a)' ) '  DPCHQ4 failed.'
  end if

  call dpchq5 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  DPCHQ5 passed.'
  else
    write ( *, '(a)' ) '  DPCHQ5 failed.'
  end if

  call pchqk1 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  PCHQK1 passed.'
  else
    write ( *, '(a)' ) '  PCHQK1 failed.'
  end if

  call pchqk2 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  PCHQK2 passed.'
  else
    write ( *, '(a)' ) '  PCHQK2 failed.'
  end if

  call pchqk3 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  PCHQK3 passed.'
  else
    write ( *, '(a)' ) '  PCHQK3 failed.'
  end if

  call pchqk4 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  PCHQK4 passed.'
  else
    write ( *, '(a)' ) '  PCHQK4 failed.'
  end if

  call pchqk5 ( lun, kprint, ipass )

  if ( ipass == 1 ) then
    write ( *, '(a)' ) '  PCHQK5 passed.'
  else
    write ( *, '(a)' ) '  PCHQK5 failed.'
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCHIP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine dpchev_test ( )

!*****************************************************************************80
!
!! DPCHEV_TEST tests DPCHEV.
!
!  Modified:
!
!    05 April 2015
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: ne = 101
  integer ( kind = 4 ), parameter :: nwk = 2 * n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) edans(4)
  real ( kind = 8 ) erans(4)
  real ( kind = 8 ) error
  real ( kind = 8 ) errord
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fans(4)
  real ( kind = 8 ) fd(ne)
  real ( kind = 8 ) fe(ne) 
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) r
  real ( kind = 8 ) rp
  logical spline
  real ( kind = 8 ) u
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xans(4)
  real ( kind = 8 ) xe(ne)

  save edans
  save erans
  save fans
  save xans

  data edans / &
    0.2932883472D+00, 0.2677039506D+00, 0.1744906562D+00, 0.0D+00 /
  data erans/ &
    -0.6075110024D-02, -0.3219009901D-02, &
    -0.946234414D-03, 0.0D+00 /
  data fans / 0.971920D+00, 0.986880D+00, 0.996560D+00, 1.0D+00 /
  data xans / -0.3D-01, -0.2D-01, -0.1D-01, 0.0D+00 /
!
!  Arithmetic statement functions for Runge's function and derivative.
!
  r ( u ) = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * u * u )
  rp ( u ) = -50.0D+00 * u * r ( u ) ** 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DPCHEV_TEST'
  write ( *, '(a)' ) '  DPCHEV evaluates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = r ( x(i) )
  end do
!
!  Set SPLINE FALSE, to choose the cubic Hermite interpolant.
!
  spline = .false.
!
!  Compute the interpolant.
!
  call dpchez ( n, x, f, d, spline, wk, nwk, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEZ returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Evaluate interpolant and derivative at NE points from -1 to 0.
!
  do i = 1, ne
    xe(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / real ( ne - 1, kind = 8 )
  end do

  call dpchev ( n, x, f, d, ne, xe, fe, fd, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEV returned error flag IERR= ', ierr
    stop 1
  end if

  do i = 1, ne
    error = fe(i) - r(xe(i))
    errord = fd(i) - rp(xe(i))
    write ( *, '(1x,d17.10,3x,d17.10,3x,d17.10,3x,d17.10)' ) &
      xe(i), fe(i), error, errord
  end do
!
!  Print reference results.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Reference results:'

  do i = 1, 4
    write ( *, '(1x,d17.10,3x,d17.10,3x,d17.10,3x,d17.10)' ) &
      xans(i), fans(i), erans(i), edans(i)
  end do

  return
end
subroutine dpchev_test02 ( )

!*****************************************************************************80
!
!! DPCHEV_TEST02 tests DPCHEV.
!
!  Discussion:
!
!    This example shows how a piecewise cubic Hermite interpolant
!    can monotonically interpolate data.
!
!  Modified:
!
!    05 April 2015
!
!  Reference:
!
!    Mathworks documentation for the MATLAB PCHIP command.
!
  implicit none

  integer ( kind = 4 ), parameter :: nd = 7
  integer ( kind = 4 ), parameter :: ni = 101
  integer ( kind = 4 ), parameter :: nwk = 2 * nd

  character ( len = 255 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 255 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) dd(nd)
  real ( kind = 8 ) di(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  character ( len = 255 ) interp_filename
  integer ( kind = 4 ) interp_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) output_filename
  logical spline
  real ( kind = 8 ) u
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  save xd
  save yd
  
  data xd / &
    -3.0D+00, -2.0D+00, -1.0D+00, 0.0D+00, 1.0D+00,  &
     2.0D+00, 3.0D+00 /
  data yd / &
    -1.0D+00, -1.0D+00, -1.0D+00, 0.0D+00, 1.0D+00,  &
     1.0D+00, 1.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DPCHEV_TEST02'
  write ( *, '(a)' ) '  DPCHEV evaluates a piecewise cubic Hermite interpolant.'
  write ( *, '(a)' ) '  Use monotonic data.'
  write ( *, '(a)' ) '  See if interpolant is monotonic.'
  write ( *, '(a)' ) ' '
!
!  Set SPLINE FALSE, to choose the cubic Hermite interpolant.
!
  spline = .false.
!
!  Compute the interpolant.
!
  call dpchez ( nd, xd, yd, dd, spline, wk, nwk, ierr )

  if ( ierr .lt. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST01 - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEZ returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Evaluate interpolant and derivative at NE points from -3 to 3.
!
  do i = 1, ni
    xi(i) = ( ( ni - i     ) * ( -3.0D+00 ) &
            + (      i - 1 ) * ( +3.0D+00 ) )  &
            / ( ni     - 1 )
  end do

  call dpchev ( nd, xd, yd, dd, ni, xi, yi, di, ierr )

  if ( ierr .ne. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST02 - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEV returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Create GNUPLOT data file.
!
  data_filename = 'test02_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, nd
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) xd(j), yd(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "'  &
    // trim ( data_filename ) // '".'
!
!  Create GNUPLOT interp file.
!
  interp_filename = 'test02_interp.txt'
  call get_unit ( interp_unit )
  open ( unit = interp_unit, file = interp_filename, status = 'replace' )
  do j = 1, ni
    write ( interp_unit, '(2x,g14.6,2x,g14.6)' ) xi(j), yi(j)
  end do
  close ( unit = interp_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics interp file "'  &
    // trim ( interp_filename ) // '".'
!
!  Plot the data and the interpolant.
!
  command_filename = 'test02_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename,  &
    status = 'replace' )

  output_filename = 'test02.png'

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' )  &
    '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' )  &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' )  &
    'set title "Data and Piecewise Cubic Hermite Interpolant"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' )  &
    'plot "' // trim ( data_filename ) //  &
    '" using 1:2 with points pt 7 ps 2 lc rgb "blue",\'
  write ( command_unit, '(a)' )  &
    '     "' // trim ( interp_filename ) //  &
    '" using 1:2 lw 3 linecolor rgb "red"'

  close ( unit = command_unit )
  write ( *, '(a)' ) '  Created graphics command file "'  &
    // trim ( command_filename ) // '".'

  return
end
subroutine dpchev_test03 ( )

!*****************************************************************************80
!
!! DPCHEV_TEST03 tests DPCHEV.
!
!  Discussion:
!
!    This example shows how a spline interpolant
!    can interpolate data.
!
!  Modified:
!
!    05 April 2015
!
!  Reference:
!
!    Mathworks documentation for the MATLAB PCHIP command.
!
  implicit none

  integer ( kind = 4 ), parameter :: nd = 7
  integer ( kind = 4 ), parameter :: ni = 101
  integer ( kind = 4 ), parameter :: nwk = 2 * nd

  character ( len = 255 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 255 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) dd(nd)
  real ( kind = 8 ) di(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  character ( len = 255 ) interp_filename
  integer ( kind = 4 ) interp_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) output_filename
  logical spline
  real ( kind = 8 ) u
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  save xd
  save yd
  
  data xd / &
    -3.0D+00, -2.0D+00, -1.0D+00, 0.0D+00, 1.0D+00,  &
     2.0D+00, 3.0D+00 /
  data yd / &
    -1.0D+00, -1.0D+00, -1.0D+00, 0.0D+00, 1.0D+00,  &
     1.0D+00, 1.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DPCHEV_TEST03'
  write ( *, '(a)' ) '  DPCHEV evaluates a cubic spline interpolant.'
  write ( *, '(a)' ) '  Use monotonic data.'
  write ( *, '(a)' ) '  See if interpolant is monotonic.'
  write ( *, '(a)' ) ' '
!
!  Set SPLINE TRUE, to choose the cubic spline interpolant.
!
  spline = .true.
!
!  Compute the interpolant.
!
  call dpchez ( nd, xd, yd, dd, spline, wk, nwk, ierr )

  if ( ierr .lt. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST01 - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEZ returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Evaluate interpolant and derivative at NE points from -3 to 3.
!
  do i = 1, ni
    xi(i) = ( ( ni - i     ) * ( -3.0D+00 ) &
            + (      i - 1 ) * ( +3.0D+00 ) )  &
            / ( ni     - 1 )
  end do

  call dpchev ( nd, xd, yd, dd, ni, xi, yi, di, ierr )

  if ( ierr .ne. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV_TEST03 - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEV returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Create GNUPLOT data file.
!
  data_filename = 'test03_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, nd
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) xd(j), yd(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "'  &
    // trim ( data_filename ) // '".'
!
!  Create GNUPLOT interp file.
!
  interp_filename = 'test03_interp.txt'
  call get_unit ( interp_unit )
  open ( unit = interp_unit, file = interp_filename, status = 'replace' )
  do j = 1, ni
    write ( interp_unit, '(2x,g14.6,2x,g14.6)' ) xi(j), yi(j)
  end do
  close ( unit = interp_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics interp file "'  &
    // trim ( interp_filename ) // '".'
!
!  Plot the data and the interpolant.
!
  command_filename = 'test03_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  output_filename = 'test03.png'

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' )  &
    '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' )  &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' )  &
    'set title "Data and Cubic Spline interpolant"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' )  &
    'plot "' // trim ( data_filename ) //  &
    '" using 1:2 with points pt 7 ps 2 lc rgb "blue",\'
  write ( command_unit, '(a)' )  &
    '     "' // trim ( interp_filename ) //  &
    '" using 1:2 lw 3 linecolor rgb "red"'

  close ( unit = command_unit )
  write ( *, '(a)' ) '  Created graphics command file "'  &
    // trim ( command_filename ) // '".'

  return
end
subroutine dpchqa_test ( )

!*****************************************************************************80
!
!! DPCHQA_TEST tests DPCHQA.
!
!  Modified:
!
!    05 April 2015
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: nwk = 2 * n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dpchqa
  real ( kind = 8 ) error
  real ( kind = 8 ) errord
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ians
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) q
  real ( kind = 8 ) qans
  real ( kind = 8 ) r
  real ( kind = 8 ) rp
  logical spline
  real ( kind = 8 ) u
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)

  save ians
  save qans

  data ians / 0 /
  data qans / 0.274679262701527D+00 /
!
!  Arithmetic statement functions for Runge's function and derivative.
!
  r ( u ) = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * u * u )
  rp ( u ) = -50.0D+00 * u * r ( u ) ** 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DPCHQA_TEST'
  write ( *, '(a)' ) '  DPCHQA integrates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = r ( x(i) )
  end do
!
!  Set SPLINE FALSE, to choose the cubic Hermite interpolant.
!
  spline = .false.
!
!  Compute the interpolant.
!
  call dpchez ( n, x, f, d, spline, wk, nwk, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEZ_TEST - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEZ returned error flag IERR= ', ierr
    stop 1
  end if
!
!  Compute integral of the interpolant over the interval [0,1] 
!
  a = 0.0D+00
  b = 1.0D+00
  q = dpchqa ( n, x, f, d, a, b, ierr )

  write ( *, '(2x,a,d20.12,3x,a,i5)' ) &
    'Integral from 0 to 1 = ', q, ' IERR = ', ierr
!
!  Print reference results.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Reference results:'
  write ( *, '(a)' ) ''

  write ( *, '(2x,a,d20.12,3x,a,i5)' ) &
    'Integral from 0 to 1 = ', qans, ' IERR = ', ians   

  return
end
