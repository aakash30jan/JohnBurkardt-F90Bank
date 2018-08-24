program main

!*****************************************************************************80
!
!! MAIN is the main program for ORBITAL_DECAY.
!
!  Discussion:
!
!    ORBITAL_DECAY uses RKF45 as an integrator for the a version
!    of the two-body problem.
!
!    Given two massive bodies subject to gravity, it is possible to write down
!    differential equations describing their motion.  These equations are
!    simpler to formulate in the frame of reference in which the center of 
!    mass of the two bodies does not move.  If one body is much more massive 
!    than the other, then our calculations in this new frame are essentially
!    the same as in the original geometry.  This is the case when one body
!    is the sun, and another a planet.  
!
!    This simulation would need to be modified if we wanted to consider
!    the behavior of two bodies of comparable mass, and expected to see
!    them both moving, or, even in the sun-planet case, if we wanted to
!    allow the sun to have a velocity while we stayed in a fixed frame
!    of reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) step_num
  real ( kind = 8 ), allocatable :: ts(:)
  real ( kind = 8 ), allocatable :: ys(:,:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ORBITAL_DECAY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  This simulation follows a small body for 20 orbits'
  write ( *, '(a)' ) '  around a relatively massive body - such as Mercury around'
  write ( *, '(a)' ) '  the sun.'
  write ( *, '(a)' ) '  Kepler''s equations for a two body system are used.'
  write ( *, '(a)' ) '  Initially, the orbit is NOT an ellipse, but as time passes,'
  write ( *, '(a)' ) '  the orbit decays into an elliptical shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use RKF45 for the ODE integrator.'

  neqn = 4
  step_num = 1000

  allocate ( ts(0:step_num) )
  allocate ( ys(neqn,0:step_num) )

  call rkf45_solve ( neqn, step_num, ts, ys )
!
!  Create graphics files for processing by gnuplot.
!
  call gnuplot_files ( neqn, step_num, ts, ys )
!
!  Free memory.
!
  deallocate ( ts )
  deallocate ( ys )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ORBITAL_DECAY'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine gnuplot_files ( neqn, step_num, ts, ys )

!*****************************************************************************80
!
!! GNUPLOT_FILES creates two files for processing by gnuplot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4) STEP_NUM, the number of steps to take.
!
!    Input, real ( kind = 8 ) TS(0:STEP_NUM), the time values.
!
!    Input, real ( kind = 8 ) YS(NEQN,0:STEP_NUM), the solution values.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) step_num 

  character ( len = 255 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 255 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) j
  real ( kind = 8 ) ts(0:step_num)
  real ( kind = 8 ) ys(neqn,0:step_num)

  call get_unit ( data_unit )
  data_filename = 'orbital_decay_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 0, step_num
    write ( data_unit, '(2x,i6,2x,5(2x,g14.6))' ) &
      j, ts(j), ys(1:neqn,j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = 'orbital_decay_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "orbital_decay.png"'
  write ( command_unit, '(a)' ) 'set xlabel "X"'
  write ( command_unit, '(a)' ) 'set ylabel "Y"'
  write ( command_unit, '(a)' ) 'set title "Orbital decay after 20 orbits"'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'set style fill solid'
  write ( command_unit, '(a)' ) 'set object 1 circle fc rgb "red"'
  write ( command_unit, '(a)' ) 'set object 1 circle at 0,0 size 0.05'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 3:5 lw 3 linecolor rgb "blue"'

  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
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
subroutine rkf45_solve ( neqn, step_num, ts, ys )

!*****************************************************************************80
!
!! RKF45_SOLVE runs the two body ODE system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4) STEP_NUM, the number of steps to take.
!
!    Output, real ( kind = 8 ) TS(0:STEP_NUM), the time values.
!
!    Output, real ( kind = 8 ) YS(NEQN,0:STEP_NUM), the solution values.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) step_num 

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag
  external kepler
  real ( kind = 8 ) relerr
  integer ( kind = 4 ) step
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) ts(0:step_num)
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)
  real ( kind = 8 ) ys(neqn,0:step_num)

  abserr = 1.0D-10
  relerr = 1.0D-10

  flag = 1

  t_start = 0.0D+00
  t_stop = 20.0D+00 * 3.895D+00

  t = 0.0D+00
  t_out = 0.0D+00

  y(1:neqn) = (/ 1.0D+00, 0.0D+00,  0.0D+00,  0.50D+00 /)

  call kepler ( t, y, yp )

  ys(1:neqn,0) = y(1:neqn)
  ts(0) = t

  do step = 1, step_num

    t = ( real ( step_num - step + 1, kind = 8 ) * t_start &
        + real (            step - 1, kind = 8 ) * t_stop ) &
        / real ( step_num,            kind = 8 )

    t_out = ( real ( step_num - step, kind = 8 ) * t_start &
            + real (            step, kind = 8 ) * t_stop ) &
            / real ( step_num,        kind = 8 )

    call r8_rkf45 ( kepler, neqn, y, yp, t, t_out, relerr, abserr, flag )

    if ( abs ( flag ) /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RKF45_SOLVE - Warning!'
      write ( *, '(a,i4,a,g14.6)' ) '  Output value of FLAG = ', flag, &
        ' at T_OUT = ', t_out
    end if

    ys(1:neqn,step) = y(1:neqn)
    ts(step) = t_out

  end do

  return
end
subroutine kepler ( t, u, up )

!*****************************************************************************80
!
!! KEPLER evaluates the right hand side of the Kepler ODE system.
!
!  Discussion:
!
!    The Kepler ODE system has the form
!
!      u' = kepler ( t, u )
!
!    where u is a vector of length 4 whose components are the position
!    and velocity of a small body orbiting a massive one.
!
!      u = [ x(t), x'(t), y(t), y'(t) ]
!      
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2013
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) U(4), the current state.
!
!    Output, real ( kind = 8 ) UP(4), the derivative of the current state.
!
  implicit none

  real ( kind = 8 ) r3
  real ( kind = 8 ) t
  real ( kind = 8 ) u(4)
  real ( kind = 8 ) up(4)

  r3 = sqrt ( ( u(1) ** 2 + u(3) ** 2 ) ** 3 )

  up = (/ u(2), -u(1) / r3, u(4), -u(3) / r3 /)

  return
end
