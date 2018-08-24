program main

!*****************************************************************************80
!
!! MAIN is the main program for SPRING_ODE2.
!
!  Discussion:
!
!    This is a revision of the SPRING_ODE code.
!
!    In this revision of the program, we want to use vectors (C arrays) to 
!    store the data, and we want to write the data out to a file in a form 
!    that Gnuplot (or other plotting programs) can use.
!
!    Hooke's law for a spring observes that the restoring force is
!    proportional to the displacement: F = - k x
!
!    Newton's law relates the force to acceleration: F = m a
!
!    Putting these together, we have
!
!      m * d^2 x/dt^2 = - k * x
!
!    We can add a damping force with coefficient c:
!
!      m * d^2 x/dt^2 = - k * x - c * dx/dt
!
!    If we write this as a pair of first order equations for (x,v), we have
!
!          dx/dt = v
!      m * dv/dt = - k * x - c * v
!
!    and now we can approximate these values for small time steps.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2012
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

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) c
  character ( len = 255 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 255 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) dt
  integer ( kind = 4 ) i
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) t_final
  real ( kind = 8 ) v(0:n)
  real ( kind = 8 ) x(0:n)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPRING_ODE2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Approximate the solution of a spring equation.'
  write ( *, '(a)' ) '  Write data to a file for use by gnuplot.'
!
!  Data
!
  m = 1.0D+00
  k = 1.0D+00
  c = 0.3D+00
  t_final = 20.0D+00
  dt = t_final / real ( n, kind = 8 )
!
!  Initial conditions.
!
  t(0) = 0.0D+00
  x(0) = 1.0D+00
  v(0) = 0.0D+00
!
!  Compute the approximate solution at equally spaced times.
!
  do i = 1, n

    t(i) = real ( i, kind = 8 ) * t_final / real ( n, kind = 8 )
    x(i) = x(i-1) + dt * v(i-1)
    v(i) = v(i-1) + ( dt / m ) * ( - k * x(i-1) - c * v(i-1) )

  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )
  data_filename = 'spring_ode2_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 0, n
    write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) t(i), x(i), v(i)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )
  command_filename = 'spring_ode2_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "xv_time.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- T --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- X(T), V(T) --->"'
  write ( command_unit, '(a)' ) &
    'set title "Position and Velocity versus Time"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue",' // &
    ' "" using 1:3 lw 3 linecolor rgb "red"'
  write ( command_unit, '(a)' ) &
    'set output "xv_phase.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X(T) --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- V(T) --->"'
  write ( command_unit, '(a)' ) &
    'set title "Position versus Velocity"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 2:3 lw 3 linecolor rgb "green"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPRING_ODE2:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

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
!    06 August 2005
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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
