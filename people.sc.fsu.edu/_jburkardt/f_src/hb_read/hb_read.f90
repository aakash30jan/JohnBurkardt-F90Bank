subroutine c4_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
  valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
  rhsval )

!*****************************************************************************80
!
!! C4_HB_DATA_READ reads the data for a C4 Harwell-Boeing file.
!
!  Discussion:
!
!    It is assumed that C4_HB_HEADER_READ has already been called, and that
!    the file is still open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Input, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Input, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Input, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Input, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Input, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
!    Output, complex ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
!    Output, complex ( kind = 4 ) RHSVAL(NROW*NRHS), HB stuff.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  complex ( kind = 4 ) rhsval(nrow*nrhs)
  integer ( kind = 4 ) rowind(nnzero)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 4 ) values(nnzero)
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if
!
!  Read right hand side values.
!
  if ( 0 < rhscrd ) then

    if ( rhstyp(1:1) == 'F' .or. rhstyp(1:1) == 'f' ) then
      read ( input, rhsfmt ) rhsval(1:nrow*nrhs)
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'C4_HB_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unsupported RHS input type.'
      stop 1
    end if

  end if

  return
end
subroutine c4_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! C4_HB_HEADER_PRINT prints the header for a C4 Harwell-Boeing file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, the number of lines of VALUE data.
!
!    Input, integer ( kind = 4 ) RHSCRD, the number of lines of RHS data.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading column
!    pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading row indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading VALUE data.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading RHS data.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero
  write ( *, '(a,i8)' ) '  Right hand sides           NRHS =   ', nrhs
  write ( *, '(a,i8)' ) '                             NRHSIX = ', nrhsix
  write ( *, '(a,i8)' ) '  Number of VALUE lines      VALCRD = ', valcrd
  write ( *, '(a,i8)' ) '  Number of RHS lines        RHSCRD = ', rhscrd
  write ( *, '(a)'    ) &
    '  RHS descriptor             RHSTYP = "' // rhstyp // '"'
  write ( *, '(a)'    ) &
    '  Format for column pointers PTRFMT = "' // trim ( ptrfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for row indices     INDFMT = "' // trim ( indfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for VALUE data      VALFMT = "' // trim ( valfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for RHS data        RHSFMT = "' // trim ( rhsfmt ) // '"'

  return
end
subroutine c4_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! C4_HB_HEADER_READ reads the header for a C4 Harwell-Boeing file.
!
!  Discussion:
!
!    This function only works for a particular kind of HB file.
!
!    It is assumed that the matrix is complex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Output, integer ( kind = 4 ) NRHS, HB stuff
!
!    Output, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Output, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Output, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Output, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Output, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Output, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Output, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
  implicit none

  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  read ( input, '(a72, a8 )' ) title, key

  read ( input, '(5i14)' ) totcrd, ptrcrd, indcrd, valcrd, rhscrd

  read ( input, '(a3, 11x, 4i14)' ) mxtype, nrow, ncol, nnzero, neltvl

  read ( input, '(a16, a16, a20, a20)' ) ptrfmt, indfmt, valfmt, rhsfmt

  if ( 0 < rhscrd ) then
    read ( input, '(a3,11x,i14,i14)' ) rhstyp, nrhs, nrhsix
  else
    rhstyp = '***'
    nrhs = 0
    nrhsix = 0
  end if

  return
end
subroutine c4_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! C4_HB_QUICK_PRINT prints a sparse matrix in C4 Harwell-Boeing format.
!
!  Discussion:
!
!    This function is designed to print the information that was read by
!    the simplified Harwell-Boeing reader C4_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, complex ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  character ( len = * ) filename
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero)
  complex ( kind = 4 ) values(nnzero)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero

  call c4_hb_structure_print ( ncol, nnzero, colptr, rowind )

  call c4_hb_values_print ( ncol, colptr, nnzero, values )

  return
end
subroutine c4_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! C4_HB_QUICK_READ reads a sparse matrix in C4 Harwell-Boeing format.
!
!  Discussion:
!
!    The Harwell-Boeing file format for sparse matrix storage is complicated,
!    and has many options.
!
!    This routine can only handle a particular simple case of the format,
!    in which the file contains a "CUA" matrix, complex, unsymmetric, assembled.
!
!    While the file may contain additional information, this routine cannot
!    retrieve it.
!
!    The user must allocate at least enough space for VALUES, ROWIND, and 
!    COLPTR before calling this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, complex ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) colptr(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  integer ( kind = 4 ) rowind(*)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 4 ) values(*)
!
!  Read header.
!
  read ( input, '(a72, a8 / 5i14 / a3, 11x, 4i14 / a16, a16, a20, a20 )' ) &
    title, key, &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
    mxtype, nrow, ncol, nnzero, neltvl, &
    ptrfmt, indfmt, valfmt, rhsfmt
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if

  return
end
subroutine c4_hb_structure_print ( ncol, nnzero, colptr, rowind )

!*****************************************************************************80
!
!! C4_HB_STRUCTURE_PRINT prints the structure of a C4 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row index of each item.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) rowind(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '------------------------------------------------------------'
  write ( *, '(a)' ) ' '
  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8,3x,a)' ) j, 'EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(3i8,3x,5i6)' ) &
            j, colptr(j), colptr(j+1)-1, rowind(klo:khi)
        else
          write ( *, '(27x,5i6)' ) rowind(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '------------------------------------------------------------'

  return
end
subroutine c4_hb_values_print ( ncol, colptr, nnzero, values )

!*****************************************************************************80
!
!! C4_HB_VALUES_PRINT prints the values of a C4 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, complex ( kind = 4 ) VALUES(NNZERO), the nonzero values of
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  complex ( kind = 4 ) values(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '--------------------------------------------------------'
  write ( *, '(a)' ) ' '

  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8)' ) j, '   EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 2
        khi = min ( klo + 1, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(2x,i6,2x,i6,2x,i6,3x)', advance = 'no' ) &
            j, colptr(j), colptr(j+1)-1
        else
          write ( *, '(27x)', advance = 'no' )
        end if
        do k = klo, khi
          write ( *, '(''('',g12.4,'','',g12.4,'')'')', advance = 'no' ) &
            values(k)
        end do
      write ( *, '(a)' ) ''
      end do


    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '--------------------------------------------------------'

  return
end
subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C4MAT_PRINT_SOME prints some of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  character ( len = * ) title
  complex ( kind = 4 ) zero

  zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( imag ( a(i,j) ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 4 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8_hb_read ( nrow, ncol, nnzero, values, rowind, colptr )

!*****************************************************************************80
!
!! C8_HB_READ reads a sparse matrix in C8 Harwell-Boeing format.
!
!  Discussion:
!
!    The Harwell-Boeing file format for sparse matrix storage is complicated,
!    and has many options.
!
!    This routine can only handle a particular simple case of the format,
!    in which the file contains a "CUA" matrix, complex, unsymmetric, assembled.
!
!    While the file may contain additional information, this routine cannot
!    retrieve it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, complex ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) colptr(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  integer ( kind = 4 ) rowind(*)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 8 ) values(*)
!
!  Read header.
!
  read ( *, '(a72, a8 / 5i14 / a3, 11x, 4i14 / a16, a16, a20, a20 )' ) &
    title, key, &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
    mxtype, nrow, ncol, nnzero, neltvl, &
    ptrfmt, indfmt, valfmt, rhsfmt
!
!  Read matrix structure.
!
  read ( *, ptrfmt ) colptr(1:ncol+1)
  read ( *, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( *, valfmt ) values(1:nnzero)
  end if

  return
end
subroutine c8_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
  valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
  rhsval )

!*****************************************************************************80
!
!! C8_HB_DATA_READ reads the data for a C8 Harwell-Boeing file.
!
!  Discussion:
!
!    It is assumed that C8_HB_HEADER_READ has already been called, and that
!    the file is still open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Input, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Input, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Input, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Input, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Input, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
!    Output, complex ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
!    Output, complex ( kind = 8 ) RHSVAL(NROW*NRHS), HB stuff.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  complex ( kind = 8 ) rhsval(nrow*nrhs)
  integer ( kind = 4 ) rowind(nnzero)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 8 ) values(nnzero)
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if
!
!  Read right hand side values.
!
  if ( 0 < rhscrd ) then

    if ( rhstyp(1:1) == 'F' .or. rhstyp(1:1) == 'f' ) then
      read ( input, rhsfmt ) rhsval(1:nrow*nrhs)
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'C8_HB_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unsupported RHS input type.'
      stop 1
    end if

  end if

  return
end
subroutine c8_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! C8_HB_HEADER_PRINT prints the header for a C8 Harwell-Boeing file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, number of lines of VALUE data.
!
!    Input, integer ( kind = 4 ) RHSCRD, number of lines of RHS data.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading column
!    pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading row indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading VALUE data.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading RHS data.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero
  write ( *, '(a,i8)' ) '  Right hand sides           NRHS =   ', nrhs
  write ( *, '(a,i8)' ) '                             NRHSIX = ', nrhsix
  write ( *, '(a,i8)' ) '  Number of VALUE lines      VALCRD = ', valcrd
  write ( *, '(a,i8)' ) '  Number of RHS lines        RHSCRD = ', rhscrd
  write ( *, '(a)'    ) &
    '  RHS descriptor             RHSTYP = "' // rhstyp // '"'
  write ( *, '(a)'    ) &
    '  Format for column pointers PTRFMT = "' // trim ( ptrfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for row indices     INDFMT = "' // trim ( indfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for VALUE data      VALFMT = "' // trim ( valfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for RHS data        RHSFMT = "' // trim ( rhsfmt ) // '"'

  return
end
subroutine c8_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! C8_HB_HEADER_READ reads the header for a C8 Harwell-Boeing file.
!
!  Discussion:
!
!    This function only works for a particular kind of HB file.
!
!    It is assumed that the matrix is complex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Output, integer ( kind = 4 ) NRHS, HB stuff
!
!    Output, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Output, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Output, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Output, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Output, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Output, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Output, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
  implicit none

  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  read ( input, '(a72, a8 )' ) title, key

  read ( input, '(5i14)' ) totcrd, ptrcrd, indcrd, valcrd, rhscrd

  read ( input, '(a3, 11x, 4i14)' ) mxtype, nrow, ncol, nnzero, neltvl

  read ( input, '(a16, a16, a20, a20)' ) ptrfmt, indfmt, valfmt, rhsfmt

  if ( 0 < rhscrd ) then
    read ( input, '(a3,11x,i14,i14)' ) rhstyp, nrhs, nrhsix
  else
    rhstyp = '***'
    nrhs = 0
    nrhsix = 0
  end if

  return
end
subroutine c8_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! C8_HB_QUICK_PRINT prints a sparse matrix in C8 Harwell-Boeing format.
!
!  Discussion:
!
!    This function is designed to print the information that was read by
!    the simplified Harwell-Boeing reader C8_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, complex ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  character ( len = * ) filename
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero)
  complex ( kind = 8 ) values(nnzero)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero

  call c8_hb_structure_print ( ncol, nnzero, colptr, rowind )

  call c8_hb_values_print ( ncol, colptr, nnzero, values )

  return
end
subroutine c8_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! C8_HB_QUICK_READ reads a sparse matrix in C8 Harwell-Boeing format.
!
!  Discussion:
!
!    The Harwell-Boeing file format for sparse matrix storage is complicated,
!    and has many options.
!
!    This routine can only handle a particular simple case of the format,
!    in which the file contains a "CUA" matrix, complex, unsymmetric, assembled.
!
!    While the file may contain additional information, this routine cannot
!    retrieve it.
!
!    The user must allocate at least enough space for VALUES, ROWIND, and 
!    COLPTR before calling this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, complex ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) colptr(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  integer ( kind = 4 ) rowind(*)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 8 ) values(*)
!
!  Read header.
!
  read ( input, '(a72, a8 / 5i14 / a3, 11x, 4i14 / a16, a16, a20, a20 )' ) &
    title, key, &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
    mxtype, nrow, ncol, nnzero, neltvl, &
    ptrfmt, indfmt, valfmt, rhsfmt
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if

  return
end
subroutine c8_hb_structure_print ( ncol, nnzero, colptr, rowind )

!*****************************************************************************80
!
!! C8_HB_STRUCTURE_PRINT prints the structure of a C8 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row index of each item.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) rowind(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '------------------------------------------------------------'
  write ( *, '(a)' ) ' '
  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8,3x,a)' ) j, 'EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(3i8,3x,5i6)' ) &
            j, colptr(j), colptr(j+1)-1, rowind(klo:khi)
        else
          write ( *, '(27x,5i6)' ) rowind(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '------------------------------------------------------------'

  return
end
subroutine c8_hb_values_print ( ncol, colptr, nnzero, values )

!*****************************************************************************80
!
!! C8_HB_VALUES_PRINT prints the values of a C8 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, complex ( kind = 8 ) VALUES(NNZERO), the nonzero values of
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  complex ( kind = 8 ) values(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '--------------------------------------------------------'
  write ( *, '(a)' ) ' '

  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8)' ) j, '   EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 2
        khi = min ( klo + 1, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(2x,i6,2x,i6,2x,i6,3x)', advance = 'no' ) &
            j, colptr(j), colptr(j+1)-1
        else
          write ( *, '(27x)', advance = 'no' )
        end if
        do k = klo, khi
          write ( *, '(''('',g12.4,'','',g12.4,'')'')', advance = 'no' ) &
            values(k)
        end do
      write ( *, '(a)' ) ''
      end do


    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '--------------------------------------------------------'

  return
end
subroutine c8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! C8MAT_PRINT prints a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  character ( len = * ) title
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

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
subroutine r4_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
  valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
  rhsval )

!*****************************************************************************80
!
!! R4_HB_DATA_READ reads the data for an R4 Harwell-Boeing file.
!
!  Discussion:
!
!    It is assumed that R4_HB_HEADER_READ has already been called, and that
!    the file is still open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Input, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Input, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Input, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Input, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Input, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
!    Output, real ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
!    Output, real ( kind = 4 ) RHSVAL(NROW*NRHS), HB stuff.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  real ( kind = 4 ) rhsval(nrow*nrhs)
  integer ( kind = 4 ) rowind(nnzero)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 4 ) values(nnzero)
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if
!
!  Read right hand side values.
!
  if ( 0 < rhscrd ) then

    if ( rhstyp(1:1) == 'F' .or. rhstyp(1:1) == 'f' ) then
      read ( input, rhsfmt ) rhsval(1:nrow*nrhs)
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'R4_HB_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unsupported RHS input type.'
      stop 1
    end if

  end if

  return
end
subroutine r4_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! R4_HB_HEADER_PRINT prints the header for an R4 Harwell-Boeing file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, number of lines of VALUE data.
!
!    Input, integer ( kind = 4 ) RHSCRD, number of lines of RHS data.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading column
!    pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading row indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading VALUE data.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading RHS data.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero
  write ( *, '(a,i8)' ) '  Right hand sides           NRHS =   ', nrhs
  write ( *, '(a,i8)' ) '                             NRHSIX = ', nrhsix
  write ( *, '(a,i8)' ) '  Number of VALUE lines      VALCRD = ', valcrd
  write ( *, '(a,i8)' ) '  Number of RHS lines        RHSCRD = ', rhscrd
  write ( *, '(a)'    ) &
    '  RHS descriptor             RHSTYP = "' // rhstyp // '"'
  write ( *, '(a)'    ) &
    '  Format for column pointers PTRFMT = "' // trim ( ptrfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for row indices     INDFMT = "' // trim ( indfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for VALUE data      VALFMT = "' // trim ( valfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for RHS data        RHSFMT = "' // trim ( rhsfmt ) // '"'

  return
end
subroutine r4_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! R4_HB_HEADER_READ reads the header for an R4 Harwell-Boeing file.
!
!  Discussion:
!
!    This function only works for a particular kind of HB file.
!
!    It is assumed that the matrix is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Output, integer ( kind = 4 ) NRHS, HB stuff
!
!    Output, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Output, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Output, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Output, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Output, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Output, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Output, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
  implicit none

  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  read ( input, '(a72, a8 )' ) title, key

  read ( input, '(5i14)' ) totcrd, ptrcrd, indcrd, valcrd, rhscrd

  read ( input, '(a3, 11x, 4i14)' ) mxtype, nrow, ncol, nnzero, neltvl

  read ( input, '(a16, a16, a20, a20)' ) ptrfmt, indfmt, valfmt, rhsfmt

  if ( 0 < rhscrd ) then
    read ( input, '(a3,11x,i14,i14)' ) rhstyp, nrhs, nrhsix
  else
    rhstyp = '***'
    nrhs = 0
    nrhsix = 0
  end if

  return
end
subroutine r4_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! R4_HB_QUICK_PRINT prints a sparse matrix in R4 Harwell-Boeing format.
!
!  Discussion:
!
!    This function is designed to print the information that was read by
!    the simplified Harwell-Boeing reader R4_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, real ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  character ( len = * ) filename
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero)
  real ( kind = 4 ) values(nnzero)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero

  call r4_hb_structure_print ( ncol, nnzero, colptr, rowind )

  call r4_hb_values_print ( ncol, colptr, nnzero, values )

  return
end
subroutine r4_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! R4_HB_QUICK_READ reads a sparse matrix in R4 Harwell-Boeing format.
!
!  Discussion:
!
!    The Harwell-Boeing file format for sparse matrix storage is complicated,
!    and has many options.
!
!    This routine can only handle a particular simple case of the format,
!    in which the file contains a "RUA" matrix, real, unsymmetric, assembled.
!
!    While the file may contain additional information, this routine cannot
!    retrieve it.
!
!    The user must allocate at least enough space for VALUES, ROWIND, and 
!    COLPTR before calling this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, real ( kind = 4 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) colptr(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  integer ( kind = 4 ) rowind(*)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 4 ) values(*)
!
!  Read header.
!
  read ( input, '(a72, a8 / 5i14 / a3, 11x, 4i14 / a16, a16, a20, a20 )' ) &
    title, key, &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
    mxtype, nrow, ncol, nnzero, neltvl, &
    ptrfmt, indfmt, valfmt, rhsfmt
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if

  return
end
subroutine r4_hb_structure_print ( ncol, nnzero, colptr, rowind )

!*****************************************************************************80
!
!! R4_HB_STRUCTURE_PRINT prints the structure of an R4 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row index of each item.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) rowind(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '------------------------------------------------------------'
  write ( *, '(a)' ) ' '
  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8,3x,a)' ) j, 'EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(3i8,3x,5i6)' ) &
            j, colptr(j), colptr(j+1)-1, rowind(klo:khi)
        else
          write ( *, '(27x,5i6)' ) rowind(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '------------------------------------------------------------'

  return
end
subroutine r4_hb_values_print ( ncol, colptr, nnzero, values )

!*****************************************************************************80
!
!! R4_HB_VALUES_PRINT prints the values of an R4 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, real ( kind = 4 ) VALUES(NNZERO), the nonzero values of
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  real ( kind = 4 ) values(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '--------------------------------------------------------'
  write ( *, '(a)' ) ' '

  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8)' ) j, '   EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(2x,i6,2x,i6,2x,i6,3x,5g12.4)' ) &
            j, colptr(j), colptr(j+1)-1, values(klo:khi)
        else
          write ( *, '(27x,5g12.4)' ) values(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '--------------------------------------------------------'

  return
end
subroutine r4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R4MAT_PRINT prints an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an MxN array of R4's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R4MAT_PRINT_SOME prints some of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an MxN array of R4's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 4 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
  valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
  rhsval )

!*****************************************************************************80
!
!! R8_HB_DATA_READ reads the data for an R8 Harwell-Boeing file.
!
!  Discussion:
!
!    It is assumed that R8_HB_HEADER_READ has already been called, and that
!    the file is still open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Input, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Input, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Input, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Input, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Input, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
!    Output, real ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
!    Output, real ( kind = 8 ) RHSVAL(NROW*NRHS), HB stuff.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  real ( kind = 8 ) rhsval(nrow*nrhs)
  integer ( kind = 4 ) rowind(nnzero)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 8 ) values(nnzero)
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if
!
!  Read right hand side values.
!
  if ( 0 < rhscrd ) then

    if ( rhstyp(1:1) == 'F' .or. rhstyp(1:1) == 'f' ) then
      read ( input, rhsfmt ) rhsval(1:nrow*nrhs)
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'R8_HB_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unsupported RHS input type.'
      stop 1
    end if

  end if

  return
end
subroutine r8_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! R8_HB_HEADER_PRINT prints the header for an R8 Harwell-Boeing file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Input, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Input, integer ( kind = 4 ) VALCRD, number of lines of VALUE data.
!
!    Input, integer ( kind = 4 ) RHSCRD, number of lines of RHS data.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading column
!    pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading row indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading VALUE data.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading RHS data.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero
  write ( *, '(a,i8)' ) '  Right hand sides           NRHS =   ', nrhs
  write ( *, '(a,i8)' ) '                             NRHSIX = ', nrhsix
  write ( *, '(a,i8)' ) '  Number of VALUE lines      VALCRD = ', valcrd
  write ( *, '(a,i8)' ) '  Number of RHS lines        RHSCRD = ', rhscrd
  write ( *, '(a)'    ) &
    '  RHS descriptor             RHSTYP = "' // rhstyp // '"'
  write ( *, '(a)'    ) &
    '  Format for column pointers PTRFMT = "' // trim ( ptrfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for row indices     INDFMT = "' // trim ( indfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for VALUE data      VALFMT = "' // trim ( valfmt ) // '"'
  write ( *, '(a)'    ) &
    '  Format for RHS data        RHSFMT = "' // trim ( rhsfmt ) // '"'

  return
end
subroutine r8_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, &
  nrhsix, valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )

!*****************************************************************************80
!
!! R8_HB_HEADER_READ reads the header for an R8 Harwell-Boeing file.
!
!  Discussion:
!
!    This function only works for a particular kind of HB file.
!
!    It is assumed that the matrix is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, character ( len = 3 ) RHSTYP, HB stuff.
!
!    Output, integer ( kind = 4 ) NRHS, the number of right hand sides.
!
!    Output, integer ( kind = 4 ) NRHSIX, HB stuff.
!
!    Output, integer ( kind = 4 ) VALCRD, the number of input lines of values.
!
!    Output, integer ( kind = 4 ) RHSCRD, the number of input lines of
!    right hand side data.
!
!    Output, character ( len = 16 ) PTRFMT, the format to read column pointer data.
!
!    Output, character ( len = 16 ) INDFMT, the format to read row index data.
!
!    Output, character ( len = 20 ) VALFMT, the format to read matrix value data.
!
!    Output, character ( len = 20 ) RHSFMT, the format to read right hand side data.
!
  implicit none

  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt

  read ( input, '(a72, a8 )' ) title, key

  read ( input, '(5i14)' ) totcrd, ptrcrd, indcrd, valcrd, rhscrd

  read ( input, '(a3, 11x, 4i14)' ) mxtype, nrow, ncol, nnzero, neltvl

  read ( input, '(a16, a16, a20, a20)' ) ptrfmt, indfmt, valfmt, rhsfmt

  if ( 0 < rhscrd ) then
    read ( input, '(a3,11x,i14,i14)' ) rhstyp, nrhs, nrhsix
  else
    rhstyp = '***'
    nrhs = 0
    nrhsix = 0
  end if

  return
end
subroutine r8_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! R8_HB_QUICK_PRINT prints a sparse matrix in R8 Harwell-Boeing format.
!
!  Discussion:
!
!    This function is designed to print the information that was read by
!    the simplified Harwell-Boeing reader R8_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Input, real ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  character ( len = * ) filename
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero)
  real ( kind = 8 ) values(nnzero)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Harwell Boeing data read from "' // trim ( filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Matrix rows                NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Matrix columns             NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  Matrix nonzeros            NNZERO = ', nnzero

  call r8_hb_structure_print ( ncol, nnzero, colptr, rowind )

  call r8_hb_values_print ( ncol, colptr, nnzero, values )

  return
end
subroutine r8_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, &
  colptr )

!*****************************************************************************80
!
!! R8_HB_QUICK_READ reads a sparse matrix in R8 Harwell-Boeing format.
!
!  Discussion:
!
!    The Harwell-Boeing file format for sparse matrix storage is complicated,
!    and has many options.
!
!    This routine can only handle a particular simple case of the format,
!    in which the file contains a "RUA" matrix, real, unsymmetric, assembled.
!
!    While the file may contain additional information, this routine cannot
!    retrieve it.
!
!    The user must allocate at least enough space for VALUES, ROWIND, and 
!    COLPTR before calling this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, the input unit associated with the file.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Output, integer ( kind = 4 ) NNZERO, the number of nonzeros.
!
!    Output, real ( kind = 8 ) VALUES(NNZERO), the nonzero values.
!
!    Output, integer ( kind = 4 ) ROWIND(NNZERO), the row indices of nonzeros.
!
!    Output, integer ( kind = 4 ) COLPTR(NCOL+1), compressed column indices.
!
  implicit none

  integer ( kind = 4 ) colptr(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neltvl
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  integer ( kind = 4 ) rowind(*)
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 8 ) values(*)
!
!  Read header.
!
  read ( input, '(a72, a8 / 5i14 / a3, 11x, 4i14 / a16, a16, a20, a20 )' ) &
    title, key, &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
    mxtype, nrow, ncol, nnzero, neltvl, &
    ptrfmt, indfmt, valfmt, rhsfmt
!
!  Read matrix structure.
!
  read ( input, ptrfmt ) colptr(1:ncol+1)
  read ( input, indfmt ) rowind(1:nnzero)
!
!  Read matrix values.
!
  if ( 0 < valcrd ) then
    read ( input, valfmt ) values(1:nnzero)
  end if

  return
end
subroutine r8_hb_structure_print ( ncol, nnzero, colptr, rowind )

!*****************************************************************************80
!
!! R8_HB_STRUCTURE_PRINT prints the structure of an R8 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) ROWIND(NNZERO), the row index of each item.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) rowind(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '------------------------------------------------------------'
  write ( *, '(a)' ) ' '
  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8,3x,a)' ) j, 'EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(3i8,3x,5i6)' ) &
            j, colptr(j), colptr(j+1)-1, rowind(klo:khi)
        else
          write ( *, '(27x,5i6)' ) rowind(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '------------------------------------------------------------'

  return
end
subroutine r8_hb_values_print ( ncol, colptr, nnzero, values )

!*****************************************************************************80
!
!! R8_HB_VALUES_PRINT prints the values of an R8 Harwell-Boeing matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns.
!
!    Input, integer ( kind = 4 ) COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeroes.
!
!    Input, real ( kind = 8 ) VALUES(NNZERO), the nonzero values of
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) colptr(ncol+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  real ( kind = 8 ) values(nnzero)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Column   Begin     End   ' // &
    '--------------------------------------------------------'
  write ( *, '(a)' ) ' '

  do j = 1, ncol

    if ( 5 < j .and. j < ncol ) then
      cycle
    end if

    if ( j == ncol .and. 6 < ncol ) then
      write ( *, '(a)' ) '  (Skipping intermediate columns...)'
    end if

    if ( colptr(j+1)-1 < colptr(j) ) then

      write ( *, '(i8)' ) j, '   EMPTY'

    else

      do klo = colptr(j), colptr(j+1)-1, 5
        khi = min ( klo + 4, colptr(j+1)-1 )
        if ( klo == colptr(j) ) then
          write ( *, '(2x,i6,2x,i6,2x,i6,3x,5g12.4)' ) &
            j, colptr(j), colptr(j+1)-1, values(klo:khi)
        else
          write ( *, '(27x,5g12.4)' ) values(klo:khi)
        end if
      end do

    end if

  end do

  write ( *, '(a)' ) &
    '                           ' // &
    '--------------------------------------------------------'

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

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
