subroutine assemble ( p, t, u0, nt, np, Pr, Ir, Jc, L, nu, ndof )

!*****************************************************************************80
!
!! ASSEMBLE assembles the local stiffness and residual into global arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P(2,NP), the coordinates of the mesh nodes.
!
!    Input, integer ( kind = 4 ) T(6,NT), the indices of the nodes in each
!    element.
!
!    Input, real ( kind = 8 ) U0(NDOF), the current solution vector.
!
!    Input, integer ( kind = 4 ) NT, the number of elements.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input/output, real ( kind = 8 ) Pr(NNZ), the values of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Ir(NNZ), the row indices of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Jc(N+1), points to the first element of each
!    column.
!
!    Input/output, real ( kind = 8 ) L(NDOF), the residual vector.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Input, integer ( kind = 4 ) NDOF, the number of degrees of freedom.
!    NDOF = 2*NP + NP0 + NE.
!
  implicit none

  integer ( kind = 4 ) ndof
  integer ( kind = 4 ), parameter :: ngp = 7
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nt

  real ( kind = 8 ) gp(3,ngp)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) Jc(*)
  integer ( kind = 4 ) j2
  real ( kind = 8 ) L(ndof)
  real ( kind = 8 ) lK(15,15)
  real ( kind = 8 ) lL(15)
  real ( kind = 8 ) nu
  real ( kind = 8 ) p(2,np)
  real ( kind = 8 ) pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real ( kind = 8 ) sh1(3,ngp)
  real ( kind = 8 ) sh1r(3,ngp)
  real ( kind = 8 ) sh1s(3,ngp)
  real ( kind = 8 ) sh2(6,ngp)
  real ( kind = 8 ) sh2r(6,ngp)
  real ( kind = 8 ) sh2s(6,ngp)
  integer ( kind = 4 ) t(6,nt)
  real ( kind = 8 ) u0(ndof)
  real ( kind = 8 ) w(ngp)
!
!  Get the quadrature rule.
!
  call quad_rule ( ngp, w, gp )
!
!  Initialize the shape functions.
!
  call init_shape ( ngp, gp, sh1, sh1r, sh1s, sh2, sh2r, sh2s )
!
!  For each element, determine the local matrix and right hand side.
!
  do it = 1, nt

    call localKL ( ngp, w, sh1, sh1r, sh1s, sh2, sh2r, sh2s, p, t(1,it), &
      u0, np, nu, lK, lL )
!
!  Add the local right hand side and matrix to the global data.
!
    do i = 1, 15

      i2 = mod ( i - 1, 6 ) + 1
      ri = t(i2,it) + np * ( ( i - 1 ) / 6 )
      L(ri) = L(ri) + lL(i)

      do j = 1, 15
        j2 = mod ( j - 1, 6 ) + 1
        rj = t(j2,it) + np * ( ( j - 1 ) / 6 )
        call sparse_set ( Pr, Ir, Jc, ri, rj, lK(i,j) )
      end do

    end do

  end do

  return
end
subroutine assemble_constraints ( e, u0, np, ne, nv, Pr, Ir, Jc, L, ndof )

!*****************************************************************************80
!
!! ASSEMBLE_CONSTRAINTS assembles the constraints.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2013
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) E(3,NE), the node, variable type (0,1,2) and
!    value of each constraint.
!
!    Input, real ( kind = 8 ) U0(NDOF), the current solution vector.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, integer ( kind = 4 ) NE, the number of constraints.
!
!    Input, integer ( kind = 4 ) NV, the number of variables.
!
!    Input/output, real ( kind = 8 ) Pr(NNZ), the values of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Ir(NNZ), the row indices of the nonzero entries
!    of the sparse matrix.
!
!    Input, integer ( kind = 4 ) Jc(N+1), points to the first element of each
!    column.
!
!    Input/output, real ( kind = 8 ) L(NDOF), the residual vector.
!    On output, the residual vector has been updated to include the constraints.
!
!    Input, integer ( kind = 4 ) NDOF, the number of degrees of freedom.
!    NDOF = 2*NP+NP0+NE.
!
  implicit none

  integer ( kind = 4 ) ndof
  integer ( kind = 4 ) ne

  real ( kind = 8 ) e(3,ne)
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) Ir(*)
  integer ( kind = 4 ) Jc(*)
  real ( kind = 8 ) l(ndof)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nv
  real ( kind = 8 ) Pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real ( kind = 8 ) u0(ndof)
  real ( kind = 8 ) value

  value = 1.0D+00

  do ie = 1, ne

    ri = nv + ie
    rj = int ( e(2,ie) ) * np + int ( e(1,ie) )

    call sparse_set ( Pr, Ir, Jc, ri, rj, value )
    call sparse_set ( Pr, Ir, Jc, rj, ri, value )

    L(rj) = L(rj) + u0(ri)
    L(ri) = u0(rj) - e(3,ie)

  end do

  return
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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
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
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

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
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine i4mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! I4MAT_HEADER_READ reads the header from an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  t    = a(1)
  a(1) = a(n)
  a(n) = t
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    t    = a(1)
    a(1) = a(n1)
    a(n1) = t

  end do

  return
end
subroutine init_shape ( ngp, gp, sh1, sh1r, sh1s, sh2, sh2r, sh2s )

!*****************************************************************************80
!
!! INIT_SHAPE evaluates the shape functions at the quadrature points.
!
!  Discussion:
!
!    The shape functions and their spatial derivatives are evaluated at the
!    Gauss points in the reference triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of Gauss quadrature points.
!
!    Input, real ( kind = 8 ) GP(3,NGP), the Gauss quadrature points.
!
!    Output, real ( kind = 8 ) SH1(3,NGP), SH1R(3,NGP), SH1S(3,NGP), the
!    P1 (pressure) shape functions, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Output, real ( kind = 8 ) SH2(6,NGP), SH2R(6,NGP), SH2S(6,NGP), the
!    P2 (velocity) shape functions, and R and S derivatives, evaluated at the
!    Gauss points.
!
  implicit none

  integer ( kind = 4 ) ngp

  real ( kind = 8 ) gp(3,ngp)
  integer ( kind = 4 ) j
  real ( kind = 8 ) sh1(3,ngp)
  real ( kind = 8 ) sh1r(3,ngp)
  real ( kind = 8 ) sh1s(3,ngp)
  real ( kind = 8 ) sh2(6,ngp)
  real ( kind = 8 ) sh2r(6,ngp)
  real ( kind = 8 ) sh2s(6,ngp)

  do j = 1, ngp

    sh1(1,j) = gp(3,j)
    sh1(2,j) = gp(1,j)
    sh1(3,j) = gp(2,j)

    sh1r(1,j) = -1.0D+00
    sh1r(2,j) =  1.0D+00
    sh1r(3,j) =  0.0D+00

    sh1s(1,j) = -1.0D+00
    sh1s(2,j) =  0.0D+00
    sh1s(3,j) =  1.0D+00

    sh2(1,j) = 1.0D+00 - 3.0D+00 * gp(1,j) - 3.0D+00 * gp(2,j) &
      + 2.0D+00 * gp(1,j) * gp(1,j) + 4.0D+00 * gp(1,j) * gp(2,j) &
      + 2.0D+00 * gp(2,j) * gp(2,j)
    sh2(2,j) = - gp(1,j) + 2.0D+00 * gp(1,j) * gp(1,j)
    sh2(3,j) = - gp(2,j) + 2.0D+00 * gp(2,j) * gp(2,j)
    sh2(4,j) = 4.0D+00 * gp(1,j) * gp(2,j)
    sh2(5,j) = 4.0D+00 * gp(2,j) - 4.0D+00 * gp(1,j) * gp(2,j) &
      - 4.0D+00 * gp(2,j) * gp(2,j)
    sh2(6,j) = 4.0D+00 * gp(1,j) - 4.0D+00 * gp(1,j) * gp(2,j) &
      - 4.0D+00 * gp(1,j) * gp(1,j)

    sh2r(1,j) = - 3.0D+00 + 4.0D+00 * gp(1,j) + 4.0D+00 * gp(2,j)
    sh2r(2,j) = - 1.0D+00 + 4.0D+00 * gp(1,j)
    sh2r(3,j) = 0.0D+00
    sh2r(4,j) = 4.0D+00 * gp(2,j)
    sh2r(5,j) = - 4.0D+00 * gp(2,j)
    sh2r(6,j) = 4.0D+00 - 8.0D+00 * gp(1,j) - 4.0D+00 * gp(2,j)

    sh2s(1,j) = - 3.0D+00 + 4.0D+00 * gp(1,j) + 4.0D+00 * gp(2,j)
    sh2s(2,j) =   0.0D+00
    sh2s(3,j) = - 1.0D+00 + 4.0D+00 * gp(2,j)
    sh2s(4,j) =   4.0D+00 * gp(1,j)
    sh2s(5,j) =   4.0D+00 - 8.0D+00 * gp(2,j) - 4.0D+00 * gp(1,j)
    sh2s(6,j) = - 4.0D+00 * gp(1,j)

  end do

  return
end
subroutine localKL ( ngp, w, sh1, sh1r, sh1s, sh2, sh2r, sh2s, p, tt, &
  u0, np, nu, lK, lL )

!*****************************************************************************80
!
!! LOCALKL assembles the local stiffness matrix and residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of quadrature points.
!
!    Input, real ( kind = 8 ) W(NGP), the quadrature weights.
!
!    Input, real ( kind = 8 ) SH1(3,NGP), SH1R(3,NGP), SH1S(3,NGP), the
!    P1 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Input, real ( kind = 8 ) SH2(6,NGP), SH2R(6,NGP), SH2S(6,NGP), the
!    P2 shape function, and R and S derivatives, evaluated at the
!    Gauss points.
!
!    Input, real ( kind = 8 ) P(2,NP), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TT(6), the nodes for the current element.
!
!    Input, real ( kind = 8 ) U0(2*NP+NP0+NE), called "U" elsewhere in the
!    program, but renamed "U0" here because we want to use "U" for the
!    horizontal velocity.  U0 contains the current finite element solution.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Output, real ( kind = 8 ) lK(15,15), the local stiffness matrix.
!
!    Output, real ( kind = 8 ) lL(15), the local residual vector.
!
  implicit none

  integer ( kind = 4 ) ngp
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) igp
  integer ( kind = 4 ) j
  real ( kind = 8 ) Jdet
  real ( kind = 8 ) Jinv(2,2)
  real ( kind = 8 ) lK(15,15)
  real ( kind = 8 ) lL(15)
  real ( kind = 8 ) mul
  real ( kind = 8 ) nu
  real ( kind = 8 ) p(3,np)
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) sh1(3,ngp)
  real ( kind = 8 ) sh1r(3,ngp)
  real ( kind = 8 ) sh1s(3,ngp)
  real ( kind = 8 ) sh2(6,ngp)
  real ( kind = 8 ) sh2r(6,ngp)
  real ( kind = 8 ) sh2s(6,ngp)
  real ( kind = 8 ) sh1x(3)
  real ( kind = 8 ) sh1y(3)
  real ( kind = 8 ) sh2x(6)
  real ( kind = 8 ) sh2y(6)
  integer ( kind = 4 ) tt(6)
  real ( kind = 8 ) u
  real ( kind = 8 ) u0(*)
  real ( kind = 8 ) ux
  real ( kind = 8 ) uy
  real ( kind = 8 ) v
  real ( kind = 8 ) vx
  real ( kind = 8 ) vy
  real ( kind = 8 ) w(ngp)
  real ( kind = 8 ) xr
  real ( kind = 8 ) xs
  real ( kind = 8 ) yr
  real ( kind = 8 ) ys
!
!  Zero out lK and lL.
!
  lK(1:15,1:15) = 0.0D+00
  lL(1:15) = 0.0D+00

  do igp = 1, ngp
!
!  Compute the Jacobian.
!
    xr = 0.0D+00
    xs = 0.0D+00
    yr = 0.0D+00
    ys = 0.0D+00

    do i = 1, 6
      xr = xr + sh2r(i,igp) * p(1,tt(i))
      xs = xs + sh2s(i,igp) * p(1,tt(i))
      yr = yr + sh2r(i,igp) * p(2,tt(i))
      ys = ys + sh2s(i,igp) * p(2,tt(i))
    end do

    Jdet = xr * ys - xs * yr

    Jinv(1,1) =  ys / Jdet
    Jinv(1,2) = -xs / Jdet
    Jinv(2,1) = -yr / Jdet
    Jinv(2,2) =  xr / Jdet
!
!  Set the X and Y derivatives of the shape functions.
!
    do i = 1, 3
      sh1x(i) = sh1r(i,igp) * Jinv(1,1) + sh1s(i,igp) * Jinv(2,1)
      sh1y(i) = sh1r(i,igp) * Jinv(1,2) + sh1s(i,igp) * Jinv(2,2)
    end do

    do i = 1, 6
      sh2x(i) = sh2r(i,igp) * Jinv(1,1) + sh2s(i,igp) * Jinv(2,1)
      sh2y(i) = sh2r(i,igp) * Jinv(1,2) + sh2s(i,igp) * Jinv(2,2)
    end do
!
!  Solution and derivatives.
!
    u  = 0.0D+00
    ux = 0.0D+00
    uy = 0.0D+00

    v  = 0.0D+00
    vx = 0.0D+00
    vy = 0.0D+00

    do i = 1, 6
      u  = u  + sh2(i,igp) * u0(tt(i))
      ux = ux + sh2x(i)    * u0(tt(i))
      uy = uy + sh2y(i)    * u0(tt(i))

      v  = v  + sh2(i,igp) * u0(np+tt(i))
      vx = vx + sh2x(i)    * u0(np+tt(i))
      vy = vy + sh2y(i)    * u0(np+tt(i))
    end do

    px = 0.0D+00
    py = 0.0D+00
    do i = 1, 3
      px = px + sh1x(i) * u0(2*np+tt(i))
      py = py + sh1y(i) * u0(2*np+tt(i))
    end do
!
!  Local K.
!
    mul = w(igp) * Jdet / 2.0D+00

    do i = 1, 6

      lL(i)   = lL(i)   + mul * ( &
        ( u * ux + v * uy + px ) * sh2(i,igp) &
        + nu * ( ux * sh2x(i) + uy * sh2y(i) ) )

      lL(6+i) = lL(6+i) + mul * ( &
        ( u * vx + v * vy + py ) * sh2(i,igp) &
        + nu * ( vx * sh2x(i) + vy * sh2y(i) ) )

      do j = 1, 6

        lK(i,j) = lK(i,j) &
          + nu * ( sh2x(i) * sh2x(j) + sh2y(i) * sh2y(j) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + nu * ( sh2x(i) * sh2x(j) + sh2y(i) * sh2y(j) ) * mul

        lK(i,j) = lK(i,j) &
          + ( u  * sh2(i,igp) * sh2x(j) &
          +   v  * sh2(i,igp) * sh2y(j) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + ( u  * sh2(i,igp) * sh2x(j) &
          +   v  * sh2(i,igp) * sh2y(j) ) * mul

        lK(i,j) = lK(i,j) &
          +  ( ux * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(i,6+j) = lK(i,6+j) &
          + ( uy * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(6+i,j) = lK(6+i,j) &
          + ( vx * sh2(i,igp) * sh2(j,igp) ) * mul

        lK(6+i,6+j) = lK(6+i,6+j) &
          + ( vy * sh2(i,igp) * sh2(j,igp) ) * mul

      end do

      do j = 1, 3
        lK(i,12+j)   = lK(i,12+j)   + ( sh2(i,igp) * sh1x(j) ) * mul
        lK(6+i,12+j) = lK(6+i,12+j) + ( sh2(i,igp) * sh1y(j) ) * mul
      end do

    end do
!
!  Local L.
!
    do i = 1, 3
      lL(12+i) = lL(12+i) + ( ux + vy ) * sh1(i,igp) * mul
      do j = 1, 6
        lK(12+i,j) =   lK(12+i,j)   + ( sh1(i,igp) * sh2x(j) ) * mul
        lK(12+i,6+j) = lK(12+i,6+j) + ( sh1(i,igp) * sh2y(j) ) * mul
      end do
    end do

  end do

  return
end
subroutine nsasm_interface ( p_file, t_file, e_file, np0, nu, &
  ndof, k_nz, k_row, k_col, k_val, l, u0  )

!*****************************************************************************80
!
!! NSASM_INTERFACE is an interface for NSASM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002
!
!  Parameters:
!
!    Input, character ( len = * ) P_FILE, contains 2 rows and NP columns of 
!    (X,Y) node coordinates.
!
!    Input, character ( len = * ) T_FILE contains 6 rows and NT columns, with 
!    each column containing the (1-based) indices of nodes forming 
!    the triangles, in a particular order.
!
!    Input, character ( len = * ) E_FILE contains 3 rows and NE columns, 
!    defining constraints on the data, including Dirichlet boundary values in 
!    particular.
!    Item #1 is a node index, 
!    Item #2 is a variable index 
!      0 = horizontal velocity,
!      1 = vertical velocity, 
!      2 = pressure;
!    Item #3 is the prescribed value.
!
!    Input, integer ( kind = 4 ) NP0, the number of pressure nodes.
!
!    Input, real ( kind = 8 ) NU, the fluid viscosity.
!
!    Output, integer ( kind = 4 ) NDOF, the number of degrees of freedom.
!
!    Output, integer ( kind = 4 ) K_NZ, the number of nonzeros in K.
!
!    Output, integer ( kind = 4 ) K_ROW(K_NZ), the row indexes of nonzeros.
!
!    Output, integer ( kind = 4 ) K_COL(K_NZ), the column indexes of nonzeros.
!
!    Output, real ( kind = 8 ) K_VAL(K_NZ), the nonzero matrix values.
!
!    Output, real ( kind = 8 ) L(NDOF), the right hand side.
!
!    Output, real ( kind = 8 ) U0(NDOF), the estimated solution.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) NE, the number of constraints.
!
!    Local, integer ( kind = 4 ) NGP, the order of Gauss quadrature.
!
!    Local, integer ( kind = 4 ) NP, the number of nodes.
!
!    Local, integer ( kind = 4 ) NT, the number of triangles.
!
  implicit none

  real ( kind = 8 ), allocatable :: e(:,:)
  character ( len = * ) e_file
  integer ( kind = 4 ) e_order
  real ( kind = 8 ), allocatable :: gp(:,:)
  integer ( kind = 4 ) k_nz
  integer ( kind = 4 ), allocatable :: k_col(:)
  integer ( kind = 4 ), allocatable :: k_row(:)
  real ( kind = 8 ), allocatable :: k_val(:)
  real ( kind = 8 ), allocatable :: l(:)
  integer ( kind = 4 ) ndof
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) ngp
  integer ( kind = 4 ) np
  integer ( kind = 4 ) np0
  integer ( kind = 4 ) nt
  real ( kind = 8 ) nu
  integer ( kind = 4 ) nv
  real ( kind = 8 ), allocatable :: p(:,:)
  character ( len = * ) p_file
  integer ( kind = 4 ) p_order
  integer ( kind = 4 ), allocatable :: t(:,:)
  character ( len = * ) t_file
  integer ( kind = 4 ) t_order
  real ( kind = 8 ), allocatable :: u0(:)
  real ( kind = 8 ), allocatable :: w(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NSASM_INTERFACE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  An interface to the NSASM function.'
!
!  Load problem data from the files.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Loading user node data from "' // trim ( p_file ) // '".'
  call r8mat_header_read ( p_file, p_order, np )
  write ( *, '(a,i6)' ) '  Data dimension = ', p_order
  write ( *, '(a,i6)' ) '  Number of nodes = ', np
  allocate ( p(1:p_order,1:np) )
  call r8mat_data_read ( p_file, p_order, np, p )

  write ( *, '(a)' ) &
    '  Loading user element data from "' // trim ( t_file ) // '".'
  call i4mat_header_read ( t_file, t_order, nt )
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Element order = ', t_order
  write ( *, '(a,i6)' ) '  Number of elements NT = ', nt
  allocate ( t(1:t_order,1:nt) )
  call i4mat_data_read ( t_file, t_order, nt, t )

  write ( *, '(a)' ) &
    '  Loading user constraint data from "' // trim ( e_file ) // '".'
  call r8mat_header_read ( e_file, e_order, ne )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Constraint array order = ', e_order
  write ( *, '(a,i6)' ) '  Number of constraints NE = ', ne
  allocate ( e(1:e_order,1:ne) )
  call r8mat_data_read ( e_file, e_order, ne, e )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of pressure nodes NP0 = ', np0
  write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
!
!  U0 is the set of finite element coefficients.
!  We set this to zero.
!  In typical usage, NSASM is used inside a Newton iteration,
!  and we are trying to improve the estimate of U0.
!
  ndof = 2 * np + np0 + ne
  write ( *, '(a,i6)' ) '  Degrees of freedom NDOF = ', ndof

  allocate ( l(1:ndof) )
  l(1:ndof) = 0.0D+00

  allocate ( u0(1:ndof) )
  u0(1:ndof) = 0.0D+00
!
!  Need to determine the number of nonzeros in K!
!
! k_nz = ?
! allocate ( k_col(1:k_nz) )
! allocate ( k_row(1:k_nz) )
! allocate ( k_val(1:k_nz) )
!
!  Construct K and L.
!  
! call assemble ( p, t, u0, nt, np, k_nz, k_col, k_row, k_val, l, nu, ndof )

  nv = 2 * np + np0

! call assemble_constraints ( e, u0, np, ne, nv, k_nz, k_col, k_row, ...
!   k_val, l, ndof )
!
!  Free memory.
!
  deallocate ( e )
  deallocate ( gp )
  deallocate ( p )
  deallocate ( t )
  deallocate ( w )

  return
end
subroutine quad_rule ( ngp, w, gp )

!*****************************************************************************80
!
!! QUAD_RULE returns the points and weights of a quadrature rule.
!
!  Discussion:
!
!    At the moment, only a 7-point rule is available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGP, the number of quadrature points.
!
!    Output, real ( kind = 8 ) W(NGP), the quadrature weights.
!
!    Output, real ( kind = 8 ) GP(3,NGP), the quadrature points.
!
  implicit none

  integer ( kind = 4 ) ngp

  real ( kind = 8 ) gp(3,ngp)
  real ( kind = 8 ), save :: gp_7(3,7) = reshape ( (/ &
    0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333334D+00, &
    0.059715871789770D+00, 0.470142064105115D+00, 0.470142064105115D+00, &
    0.470142064105115D+00, 0.059715871789770D+00, 0.470142064105115D+00, &
    0.470142064105115D+00, 0.470142064105115D+00, 0.059715871789770D+00, &
    0.797426985353087D+00, 0.101286507323457D+00, 0.101286507323457D+00, &
    0.101286507323457D+00, 0.797426985353087D+00, 0.101286507323457D+00, &
    0.101286507323457D+00, 0.101286507323457D+00, 0.797426985353087D+00 /), &
  (/ 3, 7 /) )
  real ( kind = 8 ) w(ngp)
  real ( kind = 8 ), save :: w_7(7) = (/ &
    0.225000000000000D+00, &
    0.132394152788506D+00, &
    0.132394152788506D+00, &
    0.132394152788506D+00, &
    0.125939180544827D+00, &
    0.125939180544827D+00, &
    0.125939180544827D+00 /)

  if ( ngp == 7 ) then
    w(1:ngp)      = w_7(1:ngp)
    gp(1:3,1:ngp) = gp_7(1:3,1:ngp)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUAD_RULE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of NGP.'
    stop
  end if

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8sp_print_some ( m, n, nz_num, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SP_PRINT_SOME prints some of an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT_SOME has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij(incx)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '
    write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      aij(1:inc) = 0.0D+00
!
!  Is matrix entry K actually the value of A(I,J), with J2LO <= J <= J2HI?
!  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
!  to sum up what we find.
!
      do k = 1, nz_num

        if ( i == row(k) .and. &
             j2lo <= col(k) .and. &
             col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)

        end if

      end do

      if ( any ( aij(1:inc) /= 0.0D+00 ) ) then
        write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
      end if

    end do

  end do

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length + 1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine sparse_set ( pr, ir, jc, ri, rj, val )

!*****************************************************************************80
!
!! SPARSE_SET increments an entry of the sparse matrix.
!
!  Discussion:
!
!    We know RJ, the column in which the entry occurs.  We know that entries
!    in column RJ occur in positions K1 = Jc(RJ) through K2 = JC(RJ+1)-1.
!    We may assume that the corresponding entries Ir(K1) through Ir(K2)
!    are ascending sorted and that one of these entries is equal to RI.
!
!    We now simply use binary search on Ir to locate the index K for which
!    Ir(K) = RI, and then increment Pr(K).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2011
!
!  Author:
!
!    Original C version by Per-Olof Persson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) PR(NNZ), the values of the nonzero
!    entries of the sparse matrix.
!
!    Input, integer ( kind = 4 ) IR(NNZ), the row indices of the nonzero
!    entries of the sparse matrix.
!
!    Input, integer ( kind = 4 ) JC(N+1), points to the location in IR and
!    PR of the first element of each column.
!
!    Input, integer ( kind = 4 ) RI, RJ, the row and column index of the
!    matrix entry that is to be incremented.
!
!    Input, real ( kind = 8 ) VAL, the amount to be added to the matrix entry.
!
  implicit none

  integer ( kind = 4 ) cr
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) pr(*)
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real ( kind = 8 ) val

  cr = ri
!
!  Set the range K1 <= K2 of entries in Ir to be searched.
!
  k1 = jc(rj)
  k2 = jc(rj+1) - 1
!
!  We seek the index K so that IR(K) = RI.
!
  do

    k = ( k1 + k2 ) / 2

    if ( cr < ir(k) ) then
      k2 = k - 1
    else
      k1 = k + 1
    end if

    if ( cr == ir(k) ) then
      exit
    end if

    if ( k2 < k1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_SET - Fatal error!'
      write ( *, '(a)' ) '  Could not locate the sparse matrix storage'
      write ( *, '(a)' ) '  index K for logical matrix entry (RI,RJ).'
      stop
    end if

  end do
!
!  Increment the matrix entry.
!
  pr(k) = pr(k) + val

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
