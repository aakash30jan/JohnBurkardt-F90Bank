program main

!*****************************************************************************80
!
!! FILE_COLUMN_REVERSE makes a copy of a file with each line reversed.
!
!  Example:
!
!    Input file:
!
!      This is the tale
!      of three little pigs
!      and their tails.
!
!    Output file:
!
!      elat eht si sihT
!      sgip elttil eerht fo
!      .sliat rieht dna
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Commandline argument, character ( len = * ) INPUT_FILENAME, the name of 
!    the file to be reversed.
!
!    Commandline argument, character ( len = * ) OUTPUT_FILENAME, the name of 
!    the file to be created, contained a copy of the input file with the 
!    columns reversed.
!
  implicit none

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) input_filename
  integer ( kind = 4 ) numarg
  character ( len = 255 ) output_filename

  if ( debug ) then

    call timestamp ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_REVERSE'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Reverse the columns of a file'

  end if

  numarg = iargc ( )
!
!  Argument 1 is presumed to be input filename.
!
  if ( numarg < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the input filename: '
    read ( *, '(a)' ) input_filename
  else
    iarg = 1
    call getarg ( iarg, input_filename )
  end if
!
!  Argument N is presumed to be output filename.
!
  if ( numarg < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the output filename: '
    read ( *, '(a)' ) output_filename
  else
    iarg = 2
    call getarg ( iarg, output_filename )
  end if
!
!  Create the reversed copy.
!
  call file_reverse_columns ( input_filename, output_filename )
!
!  Terminate.
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_REVERSE'
    write ( *, '(a)' ) '  Normal end of execution.'
    call timestamp ( )
  end if

  stop
end
subroutine file_reverse_columns ( input_filename, output_filename )

!*****************************************************************************80
!
!! FILE_REVERSE_COLUMNS makes a copy of a file with each lines reversed.
!
!  Example:
!
!    Input file:
!
!      This is the tale
!      of three little pigs
!      and their tails.
!
!    Output file:
!
!      elat eht si sihT
!      sgip elttil eerht fo
!      .sliat rieht dna
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file to 
!    be reversed.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file to be
!    created, contained a copy of the input file with the columns reversed.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
    return
  end if
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_filename ) // '".'
    return
  end if

  do 

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_reverse ( line )

    write ( output_unit, '(a)' ) trim ( line )

  end do

  close ( unit = input_unit )

  endfile ( unit = output_unit )
  close ( unit = output_unit )

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
subroutine s_reverse ( s )

!*****************************************************************************80
!
!! S_REVERSE reverses the characters in a string.
!
!  Example:
!
!    Input        Output
!
!    ' Cat'       'taC '
!    'Goo gol  '  'log ooG  '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to reverse.
!    Trailing blanks are ignored.
!
  implicit none

  character ch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )
 
  do i = 1, s_length / 2
    j = s_length + 1 - i
    ch     = s(i:i)
    s(i:i) = s(j:j)
    s(j:j) = ch
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
