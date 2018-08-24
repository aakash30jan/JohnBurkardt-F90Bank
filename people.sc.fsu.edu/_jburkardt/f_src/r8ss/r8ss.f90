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
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
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
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

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
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

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
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8ge_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MV multiplies an R8GE matrix by an R8VEC.
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
!    11 January 1999
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
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) )

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
subroutine r8ss_dif2 ( n, na, diag, a )

!*****************************************************************************80
!
!! R8SS_DIF2 sets up an R8SS second difference matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.   Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) NA, the dimension of the array A, which for
!    this special case will 2*N-1.
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A(2*N-1), the R8SS matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) na

  na = 0

  do j = 1, n

    if ( 1 < j ) then
      na = na + 1
      a(na) = -1.0D+00
    end if

    na = na + 1
    a(na) = 2.0D+00
    diag(j) = na

  end do

  return
end
subroutine r8ss_error ( n, na, diag, ierror )

!*****************************************************************************80
!
!! R8SS_ERROR checks dimensions for an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal 
!    elements.
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error.
!    1, N is less than 1.
!    2, NA is less than N.
!    3, DIAG(1) is not 1.
!    4, the elements of DIAG are not strictly increasing.
!    5, DIAG(N) is greater than NA.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) na

  ierror = 0

  if ( n < 1 ) then
    ierror = 1
    return
  end if

  if ( na < n ) then
    ierror = 2
    return
  end if

  if ( diag(1) /= 1 ) then
    ierror = 3
    return
  end if

  do i = 1, n - 1
    if ( diag(i+1) <= diag(i) ) then
      ierror = 4
      return
    end if
  end do

  if ( na < diag(n) ) then
    ierror = 5
    return
  end if

  return
end
subroutine r8ss_indicator ( n, na, diag, a )

!*****************************************************************************80
!
!! R8SS_INDICATOR sets up an R8SS indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.   Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) NA, the dimension of the array A, which for
!    this special case will be the maximum, ( N * ( N + 1 ) ) / 2
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8SS matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) na

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  na = 0

  do j = 1, n

    do i = 1, j
      na = na + 1
      a(na) = real ( fac * i + j, kind = 8 )
    end do

    diag(j) = na

  end do

  return
end
subroutine r8ss_mv ( n, na, diag, a, x, b )

!*****************************************************************************80
!
!! R8SS_MV multiplies an R8SS matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A*x.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 - ( diag(j) - diagold )

    do i = ilo, j - 1
      k = k + 1
      b(i) = b(i) + a(k) * x(j)
      b(j) = b(j) + a(k) * x(i)
    end do

    k = k + 1
    b(j) = b(j) + a(k) * x(j)

    diagold = diag(j)

  end do

  return
end
subroutine r8ss_print ( n, na, diag, a, title )

!*****************************************************************************80
!
!! R8SS_PRINT prints an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N 
!    diagonal elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(na)
  integer ( kind = 4 ) diag(n)
  character ( len = * ) title

  call r8ss_print_some ( n, na, diag, a, 1, 1, n, n, title )

  return
end
subroutine r8ss_print_some ( n, na, diag, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SS_PRINT_SOME prints some of an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) na
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ijm1
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

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
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
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00

        if ( j < i ) then
          if ( i == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(i-1)
          end if
          ij = diag(i)
          if ( ijm1 < ij + j - i ) then
            aij = a(ij+j-i)
          end if
        else if ( j == i ) then
          ij = diag(j)
          aij = a(ij)
        else if ( i < j ) then
          if ( j == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(j-1)
          end if
          ij = diag(j)
          if ( ijm1 < ij + i - j ) then
            aij = a(ij+i-j)
          end if
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ss_random ( n, seed, na, diag, a )

!*****************************************************************************80
!
!! R8SS_RANDOM randomizes an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
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
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8SS matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
!
!  Set the values of DIAG.
!
  diag(1) = 1
  na = 1
  do i = 2, n
    k = i4_uniform_ab ( 1, i, seed )
    diag(i) = diag(i-1) + k
    na = na + k
  end do
!
!  Now set the values of A.
!
  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j
      k = k + 1
      a(k) = r8_uniform_01 ( seed )
    end do

    diagold = diag(j)

  end do

  return
end
subroutine r8ss_to_r8ge ( n, na, diag, a, b  )

!*****************************************************************************80
!
!! R8SS_TO_R8GE copies an R8SS matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    11   0  13  0 15
!     0  22  23  0  0
!    31  32  33 34  0
!     0   0  43 44  0
!    51   0   0  0 55
!
!    A = ( 11 | 22 | 13, 23, 33 | 34, 44 | 15, 0, 0, 0, 55 )
!    NA = 12
!    DIAG = ( 1, 2, 5, 7, 12 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  b(1:n,1:n) = 0.0D+00

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j - 1
      k = k + 1
      b(i,j) = a(k)
      b(j,i) = a(k)
    end do

    k = k + 1
    b(j,j) = a(k)

    diagold = diag(j)

  end do

  return
end
subroutine r8ss_zeros ( n, na, diag, a )

!*****************************************************************************80
!
!! R8SS_ZEROS zeroes an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8SS matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na

  na = ( n * ( n + 1 ) ) / 2

  a(1:na) = 0.0D+00

  k = 0
  do i = 1, n
    k = k + i
    diag(i) = k
  end do

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
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
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
