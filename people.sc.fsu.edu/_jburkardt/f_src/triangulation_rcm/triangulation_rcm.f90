program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGULATION_RCM.
!
!  Discussion:
!
!    TRIANGULATION_RCM applies the RCM reordering to a triangulation.
!
!    The user supplies a node file and a triangle file, containing
!    the coordinates of the nodes, and the indices of the nodes that
!    make up each triangle.  Either 3-node or 6-node triangles may
!    be used.
!
!    The program reads the data, computes the adjacency information,
!    carries out the RCM algorithm to get the permutation, applies
!    the permutation to the nodes and triangles, and writes out
!    new node and triangle files that correspond to the RCM permutation.
!
!    Note that node data is normally two dimensional, that is,
!    each node has an X and Y coordinate.  In some applications, it
!    may be desirable to specify more information.  This program
!    will accept node data that includes DIM_NUM entries on each line,
!    as long as DIM_NUM is the same for each entry.
!
!  Usage:
!
!    triangulation_rcm prefix
!
!    where 'prefix' is the common filename prefix:
!
!    * prefix_nodes.txt contains the node coordinates,
!    * prefix_elements.txt contains the element definitions.
!    * prefix_rcm_nodes.txt will contain the RCM node coordinates,
!    * prefix_rcm_elements.txt will contain the RCM element definitions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_perm_bandwidth
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) :: node_filename = ' '
  character ( len = 255 ) :: element_filename = ' '
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  character ( len = 255 ) :: node_rcm_filename = ' '
  character ( len = 255 ) :: element_rcm_filename = ' '
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm_inv
  character ( len = 255 ) prefix
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_RCM:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a node dataset of NODE_NUM points in 2 dimensions.'
  write ( *, '(a)' ) &
    '  Read an associated triangulation dataset of TRIANGLE_NUM '
  write ( *, '(a)' ) '  triangles using 3 or 6 nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Apply the RCM reordering (Reverse Cuthill-McKee).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reorder the data and write it out to files.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Argument 1 is the common filename prefix.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, prefix )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_RCM:'
    write ( *, '(a)' ) '  Please enter the filename prefix.'

    read ( *, '(a)' ) prefix

  end if
!
!  Create the filenames.
!
  node_filename = trim ( prefix ) // '_nodes.txt'
  element_filename = trim ( prefix ) // '_elements.txt'
  node_rcm_filename = trim ( prefix ) // '_rcm_nodes.txt'
  element_rcm_filename = trim ( prefix ) // '_rcm_elements.txt'
!
!  Read the node data.
!
  call r8mat_header_read (  node_filename, dim_num, node_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( node_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  Number of nodes NODE_NUM  = ', node_num

  allocate ( node_xy(1:dim_num,1:node_num) )

  call r8mat_data_read ( node_filename, dim_num, node_num, node_xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( node_filename ) //'".'

  call r8mat_transpose_print_some ( dim_num, node_num, node_xy, 1, 1, &
    dim_num, 5, '  Coordinates of first 5 nodes:' )
!
!  Read the triangulation data.
!
  call i4mat_header_read (  element_filename, triangle_order, &
    triangle_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( element_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Triangle order TRIANGLE_ORDER = ', triangle_order
  write ( *, '(a,i8)' ) '  Number of triangles TRIANGLE_NUM  = ', triangle_num

  if ( triangle_order /= 3 .and. triangle_order /= 6 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_RCM - Fatal error!'
    write ( *, '(a)' ) '  This program can only handle triangulations'
    write ( *, '(a)' ) '  of orders 3 and 6.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input triangulation seems to have'
    write ( *, '(a,i8)' ) '  order = ', triangle_order
    stop
  end if

  allocate ( triangle_node(1:triangle_order,1:triangle_num) )

  call i4mat_data_read ( element_filename, triangle_order, &
    triangle_num, triangle_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( element_filename ) //'".'

  call i4mat_transpose_print_some ( triangle_order, triangle_num, &
    triangle_node, 1, 1, triangle_order, 5, '  First 5 triangles:' )
!
!  Detect and correct 0-based indexing.
!
  call mesh_base_one ( node_num, triangle_order, triangle_num, triangle_node )
!
!  Create the triangle neighbor array.
!
  allocate ( triangle_neighbor(1:triangle_order,1:triangle_num) )

  call triangulation_neighbor_triangles ( triangle_order, triangle_num, &
      triangle_node, triangle_neighbor )
!
!  Next code depends on triangle order.
!
  if ( triangle_order == 3 ) then
!
!  Count the number of adjacencies, and set up the ADJ_ROW
!  adjacency pointer array.
!
    allocate ( adj_row(1:node_num+1) )

    call triangulation_order3_adj_count ( node_num, triangle_num, &
      triangle_node, triangle_neighbor, adj_num, adj_row )
!
!  Set up the ADJ adjacency array.
!
    allocate ( adj(1:adj_num) )

    call triangulation_order3_adj_set ( node_num, triangle_num, triangle_node, &
      triangle_neighbor, adj_num, adj_row, adj )

  else if ( triangle_order == 6 ) then
!
!  Count the number of adjacencies, and set up the ADJ_ROW
!  adjacency pointer array.
!
    allocate ( adj_row(1:node_num+1) )

    call triangulation_order6_adj_count ( node_num, triangle_num, &
      triangle_node, triangle_neighbor, adj_num, adj_row )
!
!  Set up the ADJ adjacency array.
!
    allocate ( adj(1:adj_num) )

    call triangulation_order6_adj_set ( node_num, triangle_num, triangle_node, &
      triangle_neighbor, adj_num, adj_row, adj )

  end if

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth
!
!  Compute the RCM permutation.
!
  allocate ( perm(1:node_num) )

  call genrcm ( node_num, adj_num, adj_row, adj, perm )

  allocate ( perm_inv(1:node_num) )

  call perm_inverse3 ( node_num, perm, perm_inv )

  bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, &
    perm, perm_inv )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Permuted ADJ bandwidth = ', bandwidth
!
!  Permute the nodes according to the permutation vector.
!
  call r8col_permute ( dim_num, node_num, perm, node_xy )
!
!  Permute the node indices in the triangle array.
!
  do j = 1, triangle_num
    do i = 1, triangle_order
      node = triangle_node(i,j)
      triangle_node(i,j) = perm_inv ( node )
    end do
  end do
!
!  Write the nodes.
!
  call r8mat_write ( node_rcm_filename, dim_num, node_num, node_xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the node file "' &
    // trim ( node_rcm_filename ) //'".'

  call i4mat_write ( element_rcm_filename, triangle_order, &
    triangle_num, triangle_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the triangulation file "' &
    // trim ( element_rcm_filename ) //'".'
!
!  Free memory.
!
  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( perm )
  deallocate ( perm_inv )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_RCM'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function adj_bandwidth ( node_num, adj_num, adj_row, adj )

!*****************************************************************************80
!
!! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) ADJ_BANDWIDTH, the bandwidth of the adjacency
!    matrix.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) band_hi
  integer ( kind = 4 ) band_lo
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(i), adj_row(i+1)-1
      col = adj(j)
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_bandwidth = band_lo + 1 + band_hi

  return
end
function adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, perm, perm_inv )

!*****************************************************************************80
!
!! ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information and a permutation.
!
!    The routine also computes the bandwidth and the size of the envelope.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) PERM(NODE_NUM), integer PERM_INV(NODE_NUM), 
!    the permutation and inverse permutation.
!
!    Output, integer ( kind = 4 ) ADJ_PERM_BANDWIDTH, the bandwidth of the 
!    permuted adjacency matrix.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_perm_bandwidth
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) band_hi
  integer ( kind = 4 ) band_lo
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) perm(node_num)
  integer ( kind = 4 ) perm_inv(node_num)

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(perm(i)), adj_row(perm(i)+1)-1
      col = perm_inv(adj(j))
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_perm_bandwidth = band_lo + 1 + band_hi

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
!! CH_TO_DIGIT returns the value of a base 10 digit.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
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
subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
  node_num )

!*****************************************************************************80
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   05 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the 
!    connected component.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I 
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes 
!    which are to be considered.
!
!    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in 
!    the connected component, its degree.
!
!    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the connected
!    component.
!
!    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through 
!    ICCSIZE the nodes in the connected component, starting with ROOT, and 
!    proceeding by levels.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) deg(node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) ideg
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ls(node_num)
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) lvsize
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            iccsze = iccsze + 1
            ls(iccsze) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

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
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
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
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
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
subroutine genrcm ( node_num, adj_num, adj_row, adj, perm )

!*****************************************************************************80
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row
!    I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) level_num
  integer ( kind = 4 ) level_row(node_num+1)
  integer ( kind = 4 ) num
  integer ( kind = 4 ) perm(node_num)
  integer ( kind = 4 ) root

  mask(1:node_num) = 1

  num = 1

  do i = 1, node_num
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!
      call root_find ( root, adj_num, adj_row, adj, mask, level_num, &
        level_row, perm(num), node_num )
!
!  RCM orders the component using ROOT as the starting node.
!
      call rcm ( root, adj_num, adj_row, adj, mask, perm(num), iccsze, &
        node_num )

      num = num + iccsze
!
!  We can stop once every node is in one of the connected components.
!
      if ( node_num < num ) then
        return
      end if

    end if

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
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
!    Output, integer ( kind = 4 ) TABLE(M,N), the table data.
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

  character ( len = * )  input_filename
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
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*****************************************************************************80
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!    The root node chosen by the user is assigned level 1, and masked.
!    All (unmasked) nodes reachable from a node in level 1 are
!    assigned level 2 and masked.  The process continues until there
!    are no unmasked nodes adjacent to any node in the current level.
!    The number of levels may vary between 2 and NODE_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) MASK(NODE_NUM).  On input, only nodes 
!    with nonzero MASK are to be processed.  On output, those nodes which were
!    included in the level set have MASK set to 1.
!
!    Output, integer ( kind = 4 ) LEVEL_NUM, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
!    rooted level structure.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) level_num
  integer ( kind = 4 ) level_row(node_num+1)
  integer ( kind = 4 ) level(node_num)
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) lvsize
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root

  mask(root) = 0
  level(1) = root
  level_num = 0
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
    level_num = level_num + 1
    level_row(level_num) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
    do i = lbegin, lvlend

      node = level(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1)-1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          iccsze = iccsze + 1
          level(iccsze) = nbr
          mask(nbr) = 0
        end if

      end do

    end do
!
!  Compute the current level width (the number of nodes encountered.)
!  If it is positive, generate the next level.
!
    lvsize = iccsze - lvlend

    if ( lvsize <= 0 ) then
      exit
    end if

  end do

  level_row(level_num+1) = lvlend + 1
!
!  Reset MASK to 1 for the nodes in the level structure.
!
  mask(level(1:iccsze)) = 1

  return
end
subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

!*****************************************************************************80
!
!! MESH_BASE_ONE ensures that the element definition is one-based.
!
!  Discussion:
!
!    The ELEMENT_NODE array contains nodes indices that form elements.
!    The convention for node indexing might start at 0 or at 1.
!    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
!    necessary to check a given element definition and, if it is actually
!    0-based, to convert it.
!
!    This function attempts to detect 9-based node indexing and correct it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int NODE_NUM, the number of nodes.
!
!    Input, int ELEMENT_ORDER, the order of the elements.
!
!    Input, int ELEMENT_NUM, the number of elements.
!
!    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
!    definitions.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_max
  integer ( kind = 4 ) node_min
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) order

  node_min = node_num + 1
  node_max = -1

  node_min = minval ( element_node(1:element_order,1:element_num) )
  node_max = maxval ( element_node(1:element_order,1:element_num) )

  if ( node_min == 0 .and. node_max == node_num - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 0-based!'
    write ( *, '(a)' )'  This will be converted to 1-based.'
    element_node(1:element_order,1:element_num) = &
      element_node(1:element_order,1:element_num) + 1
  else if ( node_min == 1 .and. node_max == node_num  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 1-based!'
    write ( *, '(a)' )'  No conversion is necessary.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
    write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
    write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
    write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  end if

  return
end
subroutine perm_check ( n, p, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end
subroutine perm_inverse3 ( n, perm, perm_inv )

!*****************************************************************************80
!
!! PERM_INVERSE3 produces the inverse of a given permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation.
!
!    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end
subroutine r8col_permute ( m, n, p, a )

!*****************************************************************************80
!
!! R8COL_PERMUTE permutes an R8COL in place.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    The same logic can be used to permute an array of objects of any
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      M = 2
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of objects.
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
!    Input/output, real ( kind = 8 ) A(M,N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_temp(m)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp(1:m) = a(1:m,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:m,iput) = a_temp(1:m)
          exit
        end if

        a(1:m,iput) = a(1:m,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
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
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
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
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

!*****************************************************************************80
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROOT, the node that defines the connected 
!    component.  It is used as the starting point for the RCM ordering.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer ( kind = 4 ) MASK(NODE_NUM), a mask for the nodes. 
!    Only those nodes with nonzero input mask values are considered by the
!    routine.  The nodes numbered by RCM will have their mask values
!    set to zero.
!
!    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ( kind = 4 ) ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!  Local parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
!    the degree of the nodes in the section graph specified by mask and root.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) deg(node_num)
  integer ( kind = 4 ) fnbr
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) lnbr
  integer ( kind = 4 ) lperm
  integer ( kind = 4 ) lvlend
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) node
  integer ( kind = 4 ) perm(node_num)
  integer ( kind = 4 ) root
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

  mask(root) = 0

  if ( iccsze <= 1 ) then
    return
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  If no neighbors, skip to next node in this level.
!
      if ( lnbr <= fnbr ) then
        cycle
      end if
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      k = fnbr

      do while ( k < lnbr )

        l = k
        k = k + 1
        nbr = perm(k)

        do while ( fnbr < l )

          lperm = perm(l)

          if ( deg(lperm) <= deg(nbr) ) then
            exit
          end if

          perm(l+1) = lperm
          l = l-1

        end do

        perm(l+1) = nbr

      end do

    end do

  end do
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call i4vec_reverse ( iccsze, perm )

  return
end
subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*****************************************************************************80
!
!! ROOT_FIND finds a pseudo-peripheral node.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    NODE_NUM-1 from the pseudo-peripheral node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George and Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!    Norman Gibbs, William Poole, Paul Stockmeyer,
!    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, pages 236-250, 1976.
!
!    Norman Gibbs,
!    Algorithm 509: A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, pages 378-387, 1976.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer ( kind = 4 ) MASK(NODE_NUM), specifies a section subgraph.  
!    Nodes for which MASK is zero are ignored by FNROOT.
!
!    Output, integer ( kind = 4 ) LEVEL_NUM, is the number of levels in the 
!    level structure rooted at the node ROOT.
!
!    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), 
!    integer LEVEL(NODE_NUM), the level structure array pair containing the 
!    level structure found.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) iccsze
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstrt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstop
  integer ( kind = 4 ) kstrt
  integer ( kind = 4 ) level(node_num)
  integer ( kind = 4 ) level_num
  integer ( kind = 4 ) level_num2
  integer ( kind = 4 ) level_row(node_num+1)
  integer ( kind = 4 ) mask(node_num)
  integer ( kind = 4 ) mindeg
  integer ( kind = 4 ) nabor
  integer ( kind = 4 ) ndeg
  integer ( kind = 4 ) node
  integer ( kind = 4 ) root
!
!  Determine the level structure rooted at ROOT.
!
  call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
    level_row, level, node_num )
!
!  Count the number of nodes in this level structure.
!
  iccsze = level_row(level_num+1) - 1
!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
  if ( level_num == 1 ) then
    return
  end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
  if ( level_num == iccsze ) then
    return
  end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
  do

    mindeg = iccsze

    jstrt = level_row(level_num)
    root = level(jstrt)

    if ( jstrt < iccsze ) then

      do j = jstrt, iccsze

        node = level(j)
        ndeg = 0
        kstrt = adj_row(node)
        kstop = adj_row(node+1)-1

        do k = kstrt, kstop
          nabor = adj(k)
          if ( 0 < mask(nabor) ) then
            ndeg = ndeg+1
          end if
        end do

        if ( ndeg < mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do

    end if
!
!  Generate the rooted level structure associated with this node.
!
    call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
      level_row, level, node_num )
!
!  If the number of levels did not increase, accept the new ROOT.
!
    if ( level_num2 <= level_num ) then
      exit
    end if

    level_num = level_num2
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
    if ( iccsze <= level_num ) then
      exit
    end if

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
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
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
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
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

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

  logical ch_eqi
  character c
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
  if ( iterm /= 1 .and. length+1 == nchar ) then
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
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
!    I and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

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
subroutine triangulation_neighbor_triangles ( triangle_order, triangle_num, &
  triangle_node, triangle_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_NEIGHBOR_TRIANGLES determines triangle neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
!    data items.
!
!    Note that ROW is a work array allocated dynamically inside this
!    routine.  It is possible, for very large values of TRIANGLE_NUM,
!    that the necessary amount of memory will not be accessible, and the
!    routine will fail.  This is a limitation of the implementation of
!    dynamic arrays in FORTRAN90.  One way to get around this would be
!    to require the user to declare ROW in the calling routine
!    as an allocatable array, get the necessary memory explicitly with
!    an ALLOCATE statement, and then pass ROW into this routine.
!
!    Of course, the point of dynamic arrays was to make it easy to
!    hide these sorts of temporary work arrays from the poor user!
!
!    This routine was revised to store the edge data in a column
!    array rather than a row array.
!
!  Example:
!
!    The input information from TRIANGLE_NODE:
!
!    Triangle   Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in TRIANGLE_NEIGHBOR:
!
!    Triangle  Neighboring Triangles
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
!    the nodes that make up each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three
!    triangles that are direct neighbors of a given triangle.
!    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1,
!    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative
!    if there is no neighbor on that side.  In this case, that side of the
!    triangle lies on the boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  integer ( kind = 4 ) col(4,3*triangle_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  integer ( kind = 4 ) tri1
  integer ( kind = 4 ) tri2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,3,T) or (J,I,3,T),
!    (J,K,1,T) or (K,J,1,T),
!    (K,I,2,T) or (I,K,2,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, triangle_num

    i = triangle_node(1,tri)
    j = triangle_node(2,tri)
    k = triangle_node(3,tri)

    if ( i < j ) then
      col(1:4,3*(tri-1)+1) = (/ i, j, 3, tri /)
    else
      col(1:4,3*(tri-1)+1) = (/ j, i, 3, tri /)
    end if

    if ( j < k ) then
      col(1:4,3*(tri-1)+2) = (/ j, k, 1, tri /)
    else
      col(1:4,3*(tri-1)+2) = (/ k, j, 1, tri /)
    end if

    if ( k < i ) then
      col(1:4,3*(tri-1)+3) = (/ k, i, 2, tri /)
    else
      col(1:4,3*(tri-1)+3) = (/ i, k, 2, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call i4col_sort_a ( 4, 3*triangle_num, col )
!
!  Step 3. Neighboring triangles show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRIANGLE_NEIGHBOR.
!
  triangle_neighbor(1:3,1:triangle_num) = -1

  icol = 1

  do

    if ( 3 * triangle_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    tri1 = col(4,icol)
    side2 = col(3,icol+1)
    tri2 = col(4,icol+1)

    triangle_neighbor(side1,tri1) = tri2
    triangle_neighbor(side2,tri2) = tri1

    icol = icol + 2

  end do

  return
end
subroutine triangulation_order3_adj_count ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2 6
!    2:  1  *  3  6  7
!    3:  2  *  4  7  8
!    4:  3  *  5  8  9
!    5:  4  *  9 10
!    6:  1  2  *  7 11
!    7:  2  3  6  *  8 11 12
!    8:  3  4  7  *  9 12 13
!    9:  4  5  8  * 10 13 14
!   10:  5  9  * 14 15
!   11:  6  7  * 12 16
!   12:  7  8 11  * 13 16 17
!   13:  8  9 12  * 14 17 18
!   14:  9 10 13  * 15 18 19
!   15: 10 14  * 19 20
!   16: 11 12  * 17 21
!   17: 12 13 16  * 18 21 22
!   18: 13 14 17  * 19 22 23
!   19: 14 15 18  * 20 23 24
!   20: 15 19  * 24 25
!   21: 16 17  * 22
!   22: 17 18 21  * 23
!   23: 18 19 22  * 24
!   24: 19 20 23  * 25
!   25: 20 24  *
!
!    Counting the adjacencies, we get:
!
!      2, 4, 4, 4, 3,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      3, 4, 4, 4, 2.
!
!    yielding ADJ_ROW:
!
!      1,   3,   7,  11,  15,
!     18,  22,  28,  34,  40,
!     44,  48,  54,  60,  66,
!     70,  74,  80,  86,  92,
!     96,  99, 103, 107, 111,
!    113
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
!    make up each triangle, in counterclockwise order.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) triangle_neighbor(3,tri_num)
  integer ( kind = 4 ) triangle_node(3,tri_num)

  adj_num = 0
  adj_row(1:node_num) = 0
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n2) = adj_row(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n3) = adj_row(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n3) = adj_row(n3) + 1
    end if

  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  adj_num = adj_row(node_num+1) - 1

  return
end
subroutine triangulation_order3_adj_set ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row, adj )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2 6
!    2:  1  *  3  6
!    3:  2  *  4  7  8
!    4:  3  *  5  8  9
!    5:  4  *  9 10
!    6:  1  2  *  7 11
!    7:  2  3  6  *  8 11 12
!    8:  3  4  7  *  9 12 13
!    9:  4  5  8  * 10 13 14
!   10:  5  9  * 14 15
!   11:  6  7  * 12 16
!   12:  7  8 11  * 13 16 17
!   13:  8  9 12  * 14 17 18
!   14:  9 10 13  * 15 18 19
!   15: 10 14  * 19 20
!   16: 11 12  * 17 21
!   17: 12 13 16  * 18 21 22
!   18: 13 14 17  * 19 22 23
!   19: 14 15 18  * 20 23 24
!   20: 15 19  * 24 25
!   21: 16 17  * 22
!   22: 17 18 21  * 23
!   23: 18 19 22  * 24
!   24: 19 20 23  * 25
!   25: 20 24  *
!
!    Counting the adjacencies, we get:
!
!      2, 4, 4, 4, 3,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      3, 4, 4, 4, 2.
!
!    yielding ADJ_ROW:
!
!      1,   3,   7,  11,  15,
!     18,  22,  28,  34,  40,
!     44,  48,  54,  60,  66,
!     70,  74,  80,  86,  92,
!     96,  99, 103, 107, 111,
!    113
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
!    make up each triangle in counterclockwise order.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_copy(node_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) triangle_neighbor(3,tri_num)
  integer ( kind = 4 ) triangle_node(3,tri_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_row(1:node_num)
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
    end if

  end do

  return
end
subroutine triangulation_order6_adj_count ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Counting only the strict lower triangle of adjacency, we get
!
!      0, 1, 2, 1, 2,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5
!
!    Yielding ADJ_ROW:
!
!      1,  1,  2,  4,  5,
!      7, 10, 14, 18, 22,
!     24, 30, 34, 43, 47,
!     52, 55, 59, 63, 67,
!     69, 75, 79, 88, 92,
!     97.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
!    make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) triangle_neighbor(3,tri_num)
  integer ( kind = 4 ) triangle_node(6,tri_num)

  adj_num = 0
  adj_row(1:node_num) = 0
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
    n4 = triangle_node(4,triangle)
    n5 = triangle_node(5,triangle)
    n6 = triangle_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj_row(n3) = adj_row(n3) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n1) = adj_row(n1) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n2) = adj_row(n2) + 1
    adj_row(n6) = adj_row(n6) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n6) = adj_row(n6) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n6) = adj_row(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n4) = adj_row(n4) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n4) = adj_row(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n5) = adj_row(n5) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n5) = adj_row(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n6) = adj_row(n6) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n6) = adj_row(n6) + 1
    end if

  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  adj_num = adj_row(node_num+1) - 1

  return
end
subroutine triangulation_order6_adj_set ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row, adj )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Counting only the strict lower triangle of adjacency, we get
!
!      0, 1, 2, 1, 2,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5
!
!    Yielding ADJ_ROW:
!
!      1,  1,  2,  4,  5,
!      7, 10, 14, 18, 22,
!     24, 30, 34, 43, 47,
!     52, 55, 59, 63, 67,
!     69, 75, 79, 88, 92,
!     97.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
!    make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
!    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_copy(node_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) triangle_neighbor(3,tri_num)
  integer ( kind = 4 ) triangle_node(6,tri_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_row(1:node_num)
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
    n4 = triangle_node(4,triangle)
    n5 = triangle_node(5,triangle)
    n6 = triangle_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj(adj_copy(n3)) = n4
    adj_copy(n3) = adj_copy(n3) + 1
    adj(adj_copy(n4)) = n3
    adj_copy(n4) = adj_copy(n4) + 1

    adj(adj_copy(n1)) = n5
    adj_copy(n1) = adj_copy(n1) + 1
    adj(adj_copy(n5)) = n1
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n4)) = n5
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n5)) = n4
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n2)) = n6
    adj_copy(n2) = adj_copy(n2) + 1
    adj(adj_copy(n6)) = n2
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n4)) = n6
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n6)) = n4
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n5)) = n6
    adj_copy(n5) = adj_copy(n5) + 1
    adj(adj_copy(n6)) = n5
    adj_copy(n6) = adj_copy(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n1)) = n4
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n4)) = n1
      adj_copy(n4) = adj_copy(n4) + 1
      adj(adj_copy(n2)) = n4
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n4)) = n2
      adj_copy(n4) = adj_copy(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n2)) = n5
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n5)) = n2
      adj_copy(n5) = adj_copy(n5) + 1
      adj(adj_copy(n3)) = n5
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n5)) = n3
      adj_copy(n5) = adj_copy(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n1)) = n6
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n6)) = n1
      adj_copy(n6) = adj_copy(n6) + 1
      adj(adj_copy(n3)) = n6
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n6)) = n3
      adj_copy(n6) = adj_copy(n6) + 1
    end if

  end do

  return
end
