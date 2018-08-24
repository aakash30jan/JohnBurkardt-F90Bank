subroutine file_delete ( file_name )

!*****************************************************************************80
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical, parameter :: verbose = .false.
!
!  Does the file exist?
!
  if ( .not. file_exist ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  There is no file of the given name.'
    return
  end if
!
!  Is the file open?
!
  if ( file_is_open ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  The file "' // trim ( file_name ) &
      // '" is currently open.'
    write ( *, '(a)' ) '  It must be closed before it can be deleted.'
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE:'
    write ( *, '(a)' ) '  Deleting "' // trim ( file_name ) // '".'
  end if

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
function file_is_open ( file_name )

!*****************************************************************************80
!
!! FILE_IS_OPEN reports whether a file (specified by filename) is open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_IS_OPEN, is TRUE if the file is open.
!
  implicit none

  character ( len = * ) file_name
  logical file_is_open

  inquire ( file = file_name, opened = file_is_open )

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
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10 
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the 
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8GE_PRINT prints an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME prints some of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
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
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
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
      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8s3_diagonal ( m, n, nz_num, sym, row, col, a )

!*****************************************************************************80
!
!! R8S3_DIAGONAL reorders an R8S3 matrix so diagonal entries are first.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine reorders the entries of A so that the first N entries
!    are exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input/output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) t

  found = 0

  do k = 1, nz_num

    do while ( row(k) == col(k) )

      if ( row(k) == k ) then
        found = found + 1
        exit
      end if

      i = row(k)

      j = row(i)
      row(i) = row(k)
      row(k) = j

      j = col(i)
      col(i) = col(k)
      col(k) = j

      t    = a(i)
      a(i) = a(k)
      a(k) = t
 
      found = found + 1

      if ( min ( m, n ) <= found ) then
        exit
      end if
     
    end do

    if ( min ( m, n ) <= found ) then
      exit
    end if

  end do

  if ( found < min ( m, n ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_DIAGONAL - Warning!'
    write ( *, '(a,i8)' ) &
      '  Number of diagonal entries expected: ', min ( m, n )
    write ( *, '(a,i8)' ) '  Number found was ', found
  end if

  return
end
subroutine r8s3_dif2 ( m, n, nz_num, sym, row, col, a )

!*****************************************************************************80
!
!! R8S3_DIF2 sets up an R8S3 second difference matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  k = 0
!
!  Diagonal entries.
!
  do j = 1, n
    i = j
    k = k + 1
    row(k) = i
    col(k) = j
    a(k) = 2.0D+00
  end do
!
!  Offdiagonal nonzeros, by column.
!
  do j = 1, n

    if ( sym .ne. 1 ) then
      if ( 1 < j ) then
        i = j - 1
        k = k + 1
        row(k) = i
        col(k) = j
        a(k) = -1.0D+00
      end if
    end if

    if ( j + 1 <= m ) then
      i = j + 1
      k = k + 1
      row(k) = i
      col(k) = j
      a(k) = -1.0D+00
    end if

  end do

  return
end
subroutine r8s3_indicator ( m, n, nz_num, sym, row, col, a )

!*****************************************************************************80
!
!! R8S3_INDICATOR sets up an R8S3 indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    a(k) = real ( fac * i + j, kind = 8 )
  end do

  return
end
subroutine r8s3_jac_sl ( n, nz_num, sym, row, col, a, b, x, it_max )

!*****************************************************************************80
!
!! R8S3_JAC_SL solves an R8S3 system using Jacobi iteration.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine REQUIRES that the matrix be square, that the matrix
!    have nonzero diagonal entries, and that the first N entries of
!    the array A be exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution 
!    to the system.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_new(n)

  do it_num = 1, it_max
!
!  Initialize to right hand side.
!
    x_new(1:n) = b(1:n)
!
!  Subtract off-diagonal terms.
!
    do k = n + 1, nz_num
      i = row(k)
      j = col(k)
      x_new(i) = x_new(i) - a(k) * x(j)
      if ( sym == 1 ) then
        x_new(j) = x_new(j) - a(k) * x(i)
      end if
    end do
!
!  Divide by diagonal terms.
!
    x_new(1:n) = x_new(1:n) / a(1:n)
!
!  Update.
!
    x(1:n) = x_new(1:n)

  end do

  return
end
subroutine r8s3_mtv ( m, n, nz_num, sym, row, col, a, x, b )

!*****************************************************************************80
!
!! R8S3_MTV multiplies an R8VEC times an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(m)

  b(1:n) = 0.0D+00

  do k = 1, nz_num
    i = col(k)
    j = row(k)
    b(i) = b(i) + a(k) * x(j)
  end do
!
!  Handle the symmetric option.
!
  if ( sym == 1 .and. m == n ) then
    do k = 1, nz_num
      i = row(k)
      j = col(k)
      if ( i /= j ) then
        b(i) = b(i) + a(k) * x(j)
      end if
    end do
  end if

  return
end
subroutine r8s3_mv ( m, n, nz_num, sym, row, col, a, x, b )

!*****************************************************************************80
!
!! R8S3_MV multiplies an R8S3 matrix by an R8VEC.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    b(i) = b(i) + a(k) * x(j)
  end do
!
!  Handle the symmetric option.
!
  if ( sym == 1 .and. m == n ) then
    do k = 1, nz_num
      i = col(k)
      j = row(k)
      if ( i /= j ) then
        b(i) = b(i) + a(k) * x(j)
      end if
    end do
  end if

  return
end
subroutine r8s3_print ( m, n, nz_num, sym, row, col, a, title )

!*****************************************************************************80
!
!! R8S3_PRINT prints an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  call r8s3_print_some ( m, n, nz_num, sym, row, col, a, 1, 1, m, &
    n, title )

  return
end
subroutine r8s3_print_some ( m, n, nz_num, sym, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8S3_PRINT_SOME prints some of an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
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
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nz_num)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nonzero
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

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
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
      nonzero = .false.

      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do

      do k = 1, nz_num

        if ( i == row(k) .and. j2lo <= col(k) .and. col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          write ( ctemp(j2), '(g14.6)' ) aij

        else if ( sym == 1 .and. m == n .and. &
          i == col(k) .and. j2lo <= row(k) .and. row(k) <= j2hi ) then

          j2 = row(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          write ( ctemp(j2), '(g14.6)' ) aij

        end if

      end do

      if ( nonzero ) then
        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
      end if

    end do

  end do

  return
end
subroutine r8s3_random ( m, n, nz_num, sym, row, col, seed, a )

!*****************************************************************************80
!
!! R8S3_RANDOM randomizes an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sym
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) row(nz_num)

  do k = 1,  nz_num
    a(k) = r8_uniform_01 ( seed )
  end do

  return
end
subroutine r8s3_read ( input_file, m, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8S3_READ reads a square R8S3 matrix from a file.
!
!  Discussion:
!
!    This routine needs the value of NZ_NUM, which can be determined
!    by a call to R8S3_READ_SIZE.
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
!
!    Unused, integer M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) row(k), col(k), a(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  I/O error while reading record ', k
      stop 1
    end if 

  end do

  close ( unit = input_unit )

  return
end
subroutine r8s3_read_size ( input_file, m, n, nz_num )

!*****************************************************************************80
!
!! R8S3_READ_SIZE reads the size of a square R8S3 matrix from a file.
!
!  Discussion:
!
!    The value of NZ_NUM is simply the number of records in the input file.
!
!    The value of N is determined as the maximum entry in the row and column
!    vectors.
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to 
!    be read.
!
!    Output, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
  implicit none

  real ( kind = 8 ) a_k
  integer ( kind = 4 ) col_k
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row_k

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop 1
  end if

  nz_num = 0
  m = 0
  n = 0

  do

    read ( input_unit, *, iostat = ios ) row_k, col_k, a_k

    if ( ios /= 0 ) then
      exit
    end if

    nz_num = nz_num + 1
    m = max ( m, row_k )
    n = max ( n, col_k )

  end do

  close ( unit = input_unit )

  return
end
subroutine r8s3_res ( m, n, nz_num, sym, row, col, a, x, b, r )

!*****************************************************************************80
!
!! R8S3_RES computes the residual R = B-A*X for R8S3 matrices.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Input, real ( kind = 8 ) B(M), the desired result A * x.
!
!    Output, real ( kind = 8 ) R(M), the residual R = B - A * X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) col(nz_num)
  real ( kind = 8 ) r(m)
  integer ( kind = 4 ) row(nz_num)
  integer ( kind = 4 ) sym
  real ( kind = 8 ) x(n)

  call r8s3_mv ( m, n, nz_num, sym, row, col, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r8s3_to_r8ge ( m, n, nz_num, sym, row, col, a, b )

!*****************************************************************************80
!
!! R8S3_TO_R8GE copies an R8S3 matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  b(1:m,1:n) = 0.0D+00

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    b(i,j) = b(i,j) + a(k)
    if ( sym == 1 .and. m == n .and. i /= j ) then
      b(j,i) = b(j,i) + a(k)
    end if
  end do

  return
end
subroutine r8s3_write ( m, n, nz_num, sym, row, col, a, output_file )

!*****************************************************************************80
!
!! R8S3_WRITE writes a square R8S3 matrix to a file.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( output_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num
    write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)' ) row(k), col(k), a(k)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8s3_zeros ( m, n, nz_num, sym, row, col, a )

!*****************************************************************************80
!
!! R8S3_ZEROS zeroes an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) SYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) sym
  integer ( kind = 4 ) row(nz_num)

  a(1:nz_num) = 0.0D+00

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator1 vector.
!
!  Discussion:
!
!    A(1:N) = (/ 1 : N /)
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
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
