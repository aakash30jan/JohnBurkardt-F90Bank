function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
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
subroutine r8sr_dif2 ( n, nz, row, col, diag, off )

!*****************************************************************************80
!
!! R8SR_DIF2 sets up an R8SR second difference matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in the matrix.
!
!    Output, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Output, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nz2
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  do i = 1, n
    diag(i) = - 2.0D+00
  end do

  row(1) = 1
  nz2 = 0

  do i = 1, n

    if ( i == 1 ) then

      nz2 = nz2 + 1
      col(nz2) = i + 1
      off(nz2) = 1.0D+00

      row(i+1) = row(i) + 1

    else if ( i < n ) then

      nz2 = nz2 + 1
      col(nz2) = i - 1
      off(nz2) = 1.0D+00

      nz2 = nz2 + 1
      col(nz2) = i + 1
      off(nz2) = 1.0D+00

      row(i+1) = row(i) + 2

    else

      nz2 = nz2 + 1
      col(nz2) = i - 1
      off(nz2) = 1.0D+00

      row(i+1) = row(i) + 1

    end if

  end do

  return
end
subroutine r8sr_indicator ( n, nz, row, col, diag, off )

!*****************************************************************************80
!
!! R8SR_INDICATOR sets up an R8SR indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n

    j = i
    diag(i) = real ( fac * i + j, kind = 8 )

    do k = row(i), row(i+1) - 1
      j = col(k)
      off(k) = real ( fac * i + j, kind = 8 )
    end do

  end do

  return
end
subroutine r8sr_mtv ( n, nz, row, col, diag, off, x, b )

!*****************************************************************************80
!
!! R8SR_MTV multiplies an R8VEC times an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplies by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  real ( kind = 8 ) x(n)

  b(1:n) = diag(1:n) * x(1:n)

  do i = 1, n
    do k = row(i), row(i+1) - 1
      j = col(k)
      b(j) = b(j) + off(k) * x(i)
    end do
  end do

  return
end
subroutine r8sr_mv ( n, nz, row, col, diag, off, x, b )

!*****************************************************************************80
!
!! R8SR_MV multiplies an R8SR matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by the matrix.
!
!    Output, real ( kind = 8 ) B(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  real ( kind = 8 ) x(n)

  b(1:n) = diag(1:n) * x(1:n)

  do i = 1, n
    do k = row(i), row(i+1) - 1
      j = col(k)
      b(i) = b(i) + off(k) * x(j)
    end do
  end do

  return
end
subroutine r8sr_print ( n, nz, row, col, diag, off, title )

!*****************************************************************************80
!
!! R8SR_PRINT prints an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of
!    the element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  character ( len = * ) title

  call r8sr_print_some ( n, nz, row, col, diag, off, 1, 1, n, n, title )

  return
end
subroutine r8sr_print_some ( n, nz, row, col, diag, off, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SR_PRINT_SOME prints some of an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nz)
  character ( len = 14 ) ctemp(incx)
  real ( kind = 8 ) diag(n)
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
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
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
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
!  1) Assume everything is zero.
!
      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do
!
!  2) Insert the diagonal entry, if appropriate.
!
      if ( j2lo <= i .and. i <= j2hi ) then
        j2 = i - j2lo + 1
        aij = diag(i)
        write ( ctemp(j2), '(g14.6)' ) aij
      end if
!
!  3) Now examine all the offdiagonal entries.
!
      do k = row(i), row(i+1) - 1
        if ( j2lo <= col(k) .and. col(k) <= j2hi ) then 
          j2 = col(k) - j2lo + 1
          aij = off(k)
          write ( ctemp(j2), '(g14.6)' ) aij
        end if
      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sr_random ( n, nz, row, col, diag, off, seed )

!*****************************************************************************80
!
!! R8SR_RANDOM randomizes an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements 
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  integer ( kind = 4 ) seed

  do i = 1, n
    diag(i) = r8_uniform_01 ( seed )
    do j = row(i), row(i+1) - 1
      off(j) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8sr_to_r8ge ( n, nz, row, col, diag, off, b )

!*****************************************************************************80
!
!! R8SR_TO_R8GE converts an R8SR matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
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
!    21 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SR_TO_R8GE - Fatal error!'
    write ( *, '(a,i8)' ) '  N is less than or equal to zero, N = ', n
    stop 1
  end if

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    b(i,i) = diag(i)
  end do

  do i = 1, n
    do j = row(i), row(i+1) - 1
      b(i,col(j)) = off(j)
    end do
  end do

  return
end
subroutine r8sr_zeros ( n, nz, row, col, diag, off )

!*****************************************************************************80
!
!! R8SR_ZEROS zeroes an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements 
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  diag(1:n) = 0.0D+00
  off(1:nz) = 0.0D+00

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
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
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
