program main

!*****************************************************************************80
!
!! MAIN is the main program for EXTRACT.
!
!  Discussion:
!
!    EXTRACT extracts a module from a FORTRAN file.
!
!    A "module" is a function, module, program, or subroutine 
!    specified by name.
!
!    The command
!
!      extract capchr extract.f90
!
!    searches the file "extract.f90" for a subroutine or module named
!    "capchr", and if successful, writes the text of that module to
!    the file "capchr.f90".
!
!    The module name may include an extension.  Thus, the command
!
!      extract capchr.txt extract.f90
!
!    is the same as the previous command, but the output file is
!    forced to be called "capchr.txt".
!
!    If the module name does not include an extension, then then
!    extension is taken from the input file.  The command
!
!      extract capchr extract.f
!
!    writes the text to the file "capchr.f".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    character ( len = * ) MODULE, the name of the module to be extracted.
!
!    character ( len = * ) INPUT_FILE, the file to be searched for the routine.
!
  implicit none

  logical found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) input_extension
  character ( len = 255 ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  character ( len = 255 ) module
  character ( len = 255 ) module_file_name
  character ( len = 255 ) module_name
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) output_unit
  logical, parameter ::  verbose = .false.

  if ( verbose ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTRACT'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Extract a named module from a FORTRAN file.'
  end if
!
!  Count the number of command line arguments.
!
  numarg = iargc ( )
!
!  Get the module name.
!
  if ( 1 <= numarg ) then

    iarg = 1
    call getarg ( iarg, module )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the module name?'
    read ( *, '(a)' ) module
    if ( module == ' ' ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTRACT - Fatal error!'
      write ( *, '(a)' ) '  The module name was blank.'
      stop
    end if

  end if
!
!  Get the input file name.
!
  if ( 2 <= numarg ) then

    iarg = 2
    call getarg ( iarg, input_file )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the name of the input file?'
    read ( *, '(a)' ) input_file
    if ( input_file == ' ' ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTRACT - Fatal error!'
      write ( *, '(a)' ) '  The input file name was blank.'
      stop
    end if

  end if

  call file_ext ( input_file, i, j )

  if ( i == 0 ) then
    input_extension = ' '
  else if ( 1 < i ) then
    input_extension = input_file(i-1:j)
  else
    input_extension = ' '
  end if
!
!  If ( MODULE = module_name // module_extension ) then
!
!    module_file_name := module_name // module_extension
!
!  else if ( MODULE = module_name ) then
!
!    module_file_name := module_name // input_extension
!
  call file_ext ( module, i, j )

  if ( i == 0 ) then
    module_name = module
    call s_cat ( module_name, input_extension, module_file_name )
  else if ( 1 < i ) then
    module_name = module(1:i-2)
    module_file_name = module
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTRACT found an illegal period in'
    write ( *, '(a)' ) 'the module name ' // trim ( module )
  end if
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) trim ( input_file )
    stop
  end if
!
!  Search for the line that begins the module.
!
  call module_find ( input_unit, found, line, module_name )
!
!  If the starting line was found, open the output file,
!  and copy the appropriate text to it.
!
  if ( found ) then

    call get_unit ( output_unit )

    open ( unit = output_unit, file = module_file_name, status = 'replace', &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTRACT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(4x,a)' ) trim ( module_file_name )
      close ( unit = input_unit )
      stop
    end if

    call module_write ( input_unit, output_unit, line, module_file_name )

    close ( unit = output_unit )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  Could not find the module ' // trim ( module_name )

  end if

  close ( unit = input_unit )
!
!  Terminate.
!
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTRACT:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_ext ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_EXT determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME  I  J
!
!    bob.for    5  7
!    N.B.C.D    7  7
!    Naomi.     0  0
!    Arthur     0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last 
!    characters in the file extension.
!
!    If at least one period occurs in the filename, and at least one
!    nonblank character follows that period, then I will be the index
!    of the first character after the period, and J the index of the
!    last nonblank character after the period.  The extension is
!    therefore equal to FILE_NAME(I:J).
!
!    Otherwise, I and J will be returned as 0, indicating that the file
!    has no extension.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

    if ( i == j ) then
      i = 0
      j = 0
    else
      i = i + 1
    end if

  else

    j = 0

  end if

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
!    18 September 2005
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
subroutine module_find ( input_unit, found, line, module )

!*****************************************************************************80
!
!! MODULE_FIND searches a file for the first line of a given module.
!
!  Discussion:
!
!    It is assumed that the input file has been opened, and is positioned
!    at line 1.  On return from this routine, if the module was found,
!    the file is positioned just after the declaration line of the module.  
!    Otherwise, the file is positioned at the end of the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the FORTRAN unit associated 
!    with the input file.
!
!    Output, logical FOUND, is TRUE if the first line was found.
! 
!    Output, character ( len = 80 ) LINE, the first line, if found, or ' '.
!
!    Input, character ( len = * ) MODULE, the name of the module.
!
  implicit none

  logical found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = 255 ) line2
  character ( len = * ) module
  logical s_eqi
  character ( len = 120 ) s1
  character ( len = 120 ) s2
  character ( len = 120 ) s3

  found = .false.
!
!  Read a line from the input file.
!
  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      line = ' '
      return
    end if
!
!  Make a "cleaned up" version of the line that is easier to process.
!
    line2 = line
    call s_cap ( line2 )
    call s_blanks_delete ( line2 )
!
!  Possible patterns:
!
!    DOUBLE PRECISION FUNCTION name
!
!    BLOCK DATA name
!    BLOCKDATA name
!
!    CHARACTER FUNCTION name
!    COMPLEX   FUNCTION name
!    INTEGER   FUNCTION name
!    LOGICAL   FUNCTION name
!    REAL      FUNCTION name
!
!    RECURSIVE FUNCTION name
!
!    FUNCTION name
!    PROGRAM name
!    SUBROUTINE name
!    MODULE name
!
!  We don't try to catch more elaborate declarations such as:
!
!    REAL*2 FUNCTION name
!    RECURSIVE CHARACTER FUNCTION name
!
    call s_split ( line2, 'DOUBLE PRECISION FUNCTION', s1, s2, s3 )

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'BLOCK DATA', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'BLOCKDATA', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'CHARACTER FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'COMPLEX FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'INTEGER FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'LOGICAL FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'MODULE', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'PROGRAM', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'REAL FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'RECURSIVE FUNCTION', s1, s2, s3 )
    end if

    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'SUBROUTINE', s1, s2, s3 )
    end if
!
!  This check has to occur AFTER the checks for typed functions.
!
    if ( s1 /= ' ' .or. s2 == ' ' ) then
      call s_split ( line2, 'FUNCTION', s1, s2, s3 )
    end if
!
!  If there is no prefix to the string,
!  and we matched the string,
!  and there's something afterwards, ...
!
    if ( s1 == ' ' .and. s2 /= ' ' .and. s3 /= ' ' ) then

      i = index ( s3, '(' )

      if ( i /= 0 ) then
        s3 = s3(1:i-1)
      end if

      if ( s_eqi ( s3, module ) ) then
        found = .true.
        return
      end if

    end if

  end do

  return
end
subroutine module_write ( input_unit, output_unit, line, module_file_name )

!*****************************************************************************80
!
!! MODULE_WRITE writes out the lines of a file until 'END' is reached.
!
!  Discussion:
!
!    On input, the input file is assumed to be positioned just after the
!    first line of the module.  On return, the input file is positioned
!    at the end of file, or just after the "END" statement.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the FORTRAN unit associated with 
!    the input file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the FORTRAN unit associated with 
!    the output file.
!
!    Input/workspace, character ( len = 80 ) LINE.  On input, LINE contains the
!    first line of a module.  Thereafter, LINE is used as temporary
!    storage to contain the successive lines of the module.  On return,
!    LINE contains the first line after the module, or nothing if the
!    end of file was reached.
!
!    Input, character ( len = * ) MODULE_FILE_NAME, the name of the output file.
!
  implicit none

  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) module_file_name
  integer ( kind = 4 ) nline
  integer ( kind = 4 ) output_unit
  logical s_eqi
  character ( len = 11 ) s_of_i4

  nline = 0
!
!  Write the input line to the output file.
!
  write ( output_unit, '(a)' ) trim ( line )
  nline = nline + 1
!
!  Now read the next line from the file, and if it's not "END", write
!  it to the output file.
!
  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( output_unit, '(a)' ) trim ( line )
    nline = nline + 1
!
!  Do you think this was the last line of the module?
!
    call s_blank_delete ( line )

    if ( s_eqi ( line,       'END' ) .or. &
         s_eqi ( line(1:9),  'ENDBLOCKDATA' ) .or. &
         s_eqi ( line(1:11), 'ENDFUNCTION' ) .or. &
         s_eqi ( line(1:9),  'ENDMODULE' ) .or. &
         s_eqi ( line(1:13), 'ENDSUBROUTINE' ) .or. &
         s_eqi ( line(1:4),  'END!' ) ) then

       exit

    end if

  end do

  write ( *, '(a)' ) trim ( s_of_i4 ( nline ) ) // &
    ' lines written to ' // trim ( module_file_name )

  return
end
subroutine s_before_ss_copy ( s, ss, s2 )

!*****************************************************************************80
!
!! S_BEFORE_SS_COPY copies a string up to a given substring.
!
!  Discussion:
!
!    S and S2 can be the same object, in which case the string is
!    overwritten by a copy of itself up to the substring, followed
!    by blanks.
!
!  Example:
!
!    Input:
!
!      S = 'ABCDEFGH'
!      SS = 'EF'
!
!    Output:
!
!      S2 = 'ABCD'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be copied.
!
!    Input, character ( len = * ) SS, the substring before which the copy stops.
!
!    Output, character ( len = * ) S2, the copied portion of S.
!
  implicit none

  integer ( kind = 4 ) last
  integer ( kind = 4 ) last_s2
  character ( len = * ) s
  character ( len = * ) s2
  character ( len = * ) ss
!
!  Find the first occurrence of the substring.
!
  last = index ( s, ss )
!
!  If the substring doesn't occur at all, behave as though it begins
!  just after the string terminates.
!
!  Now redefine LAST to point to the last character to copy before
!  the substring begins.
!
  if ( last == 0 ) then
    last = len ( s )
  else
    last = last - 1
  end if
!
!  Now adjust again in case the copy holder is "short".
!
  last_s2 = len ( s2 )

  last = min ( last, last_s2 )
!
!  Copy the beginning of the string.
!  Presumably, compilers now understand that if LAST is 0, we don't
!  copy anything.
!  Clear out the rest of the copy.
!
  s2(1:last) = s(1:last)
  s2(last+1:last_s2) = ' '

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nchar
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len ( s )
  j = 0
  newchr = ' '

  do i = 1, nchar

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len ( s )

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // trim ( s2 )
 
  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_index_last ( string, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in STRING at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If STRING is of length 80, and SUB is of
!    length 80, then if STRING = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either STRING or SUB is all blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    STRING.  Otherwise S_INDEX_LAST = I, where STRING(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  integer ( kind = 4 ) s_index_last
  character ( len = * ) string
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( string )
  llen2 = len_trim ( sub )
!
!  In case STRING or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( string )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( string(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
function s_indexi ( s, sub )

!*****************************************************************************80
!
!! S_INDEXI is a case-insensitive INDEX function.
!
!  Discussion:
!
!    The function returns the location in the string at which the
!    substring SUB is first found, or 0 if the substring does not
!    occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    Because of the suppression of trailing blanks, this routine cannot be
!    used to find, say, the first occurrence of the two-character
!    string 'A '.  However, this routine treats as a special case the
!    occurrence where S or SUB is entirely blank.  Thus you can
!    use this routine to search for occurrences of double or triple blanks
!    in a string, for example, although INDEX itself would be just as
!    suitable for that problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEXI.  0 if SUB does not occur in
!    the string.  Otherwise S(S_INDEXI:S_INDEXI+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the first place
!    this happens.  However, note that this routine ignores case,
!    unlike the standard FORTRAN INDEX function.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  logical s_eqi
  integer ( kind = 4 ) s_indexi
  character ( len = * ) sub

  s_indexi = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do i = 1, llen1 + 1 - llen2

    if ( s_eqi ( s(i:i+llen2-1), sub ) ) then
      s_indexi = i
      return
    end if

  end do

  return
end
function s_of_i4 ( i )

!*****************************************************************************80
!
!! S_OF_I4 converts an I4 to a left-justified string.
!
!  Example:
!
!         I  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an I4 to be converted.
!
!    Output, character ( len = 11 ) S_OF_I4, the representation of the
!    I4.  The representation will be left-justified.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  character ( len = 11 ) s
  character ( len = 11 ) s_of_i4

  s = ' '

  ilo = 1
  ihi = 11
!
!  Make a copy of the I4.
!
  ival = i
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer ( kind = 4 ) goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do j = 1, ihi
        s(j:j) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  s_of_i4 = s

  return
end
subroutine s_split ( s, sub, s1, s2, s3 )

!*****************************************************************************80
!
!! S_SPLIT divides a string into three parts, given the middle.
!
!  Discussion:
!
!    This version of the routine is case-insensitive.
!
!  Example:
!
!    Input:
!
!      S = 'aBCdEfgh'
!      S2 = 'eF'
!
!    Output:
!
!      S1 = 'aBCd'
!      S2 =  'gh'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be analyzed.
!
!    Input, character ( len = * ) SUB, the substring used to "split" S.
!    Trailing blanks in SUB are ignored.
!
!    Output, character ( len = * ) S1, the entries in the string, up
!    to, but not including, the first occurrence, if any,
!    of SUB.  If SUB occurs immediately, then S1 = ' '.
!    If SUB is not long enough, trailing entries will be lost.
!
!    Output, character ( len = * ) S2, the part of the string that matched SUB.
!    If S2 is ' ', then there wasn't a match.
!
!    Output, character ( len = * ) S3, the part of the string after the match.
!    If there was no match, then S3 is blank.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenm
  integer ( kind = 4 ) lens
  character ( len = * ) s
  integer ( kind = 4 ) s_indexi
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3
  character ( len = * ) sub

  lens = len_trim ( s )

  lenm = len_trim ( sub )
  if ( lenm == 0 ) then
    lenm = 1
  end if

  i = s_indexi ( s, sub )
!
!  The substring did not occur.
!
  if ( i == 0 ) then
    s1 = s
    s2 = ' '
    s3 = ' '
!
!  The substring begins immediately.
!
  else if ( i == 1 ) then
    s1 = ' '
    s2 = s(1:lenm)
    s3 = s(lenm+1:)
!
!  What am I checking here?
!
  else if ( lens < i + lenm ) then
    s1 = s
    s2 = ' '
    s3 = ' '
!
!  The substring occurs in the middle.
!
  else
    s1 = s(1:i-1)
    s2 = s(i:i+lenm-1)
    s3 = s(i+lenm: )
  end if
!
!  Drop leading blanks.
!
  s1 = adjustl ( s1 )
  s2 = adjustl ( s2 )
  s3 = adjustl ( s3 )

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
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Special cases:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( lenc < next )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do

  word = s(ilo:next-1)

  return
end
