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
subroutine r8bto_dif2 ( m, l, a )

!*****************************************************************************80
!
!! R8BTO_DIF2 sets up an R8BTO second difference matrix.
!
!  Discussion:
!
!    To get the second difference matrix, it is assumed that M is 1!
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 97, 129 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!    For the second different matrix, M should be 1.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Output, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  real ( kind = 8 ) value
!
!  Blocks 1 to L form the first row.
!
  j = 0
  do k = 1, l

    if ( k == 1 ) then
      value = 2.0D+00
    else if ( k == 2 ) then
      value = -1.0D+00
    else
      value = 0.0D+00
    end if

    do j2 = 1, m
      j = j + 1
      do i = 1, m
        a(i,j2,k) = value
      end do
    end do
  end do
!
!  Blocks L+1 through 2*L-1 form the remainder of the first column.
!
  i = m

  do k = l + 1, 2 * l - 1

    if ( k == l + 1 ) then
      value = -1.0D+00
    else
      value = 0.0D+00
    end if

    do i2 = 1, m
      i = i + 1
      do j = 1, m
        a(i2,j,k) = value
      end do
    end do

  end do

  return
end
subroutine r8bto_indicator ( m, l, a )

!*****************************************************************************80
!
!! R8BTO_INDICATOR sets up an R8BTO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 97, 129 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Output, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( m * l ) + 1 )
!
!  Blocks 1 to L form the first row.
!
  j = 0

  do k = 1, l

    do j2 = 1, m
      j = j + 1
      do i = 1, m
        a(i,j2,k) = real ( fac * i + j, kind = 8 )
      end do
    end do
  end do
!
!  Blocks L+1 through 2*L-1 form the remainder of the first column.
!
  i = m

  do k = l + 1, 2 * l - 1

    do i2 = 1, m
      i = i + 1
      do j = 1, m
        a(i2,j,k) = real ( fac * i + j, kind = 8 )
      end do
    end do

  end do

  return
end
subroutine r8bto_mtv ( m, l, a, x, b )

!*****************************************************************************80
!
!! R8BTO_MTV multiplies an R8VEC by an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 163, 122, 121, 130, 87, 96 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or 
!    column of A.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) X(M*L), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M*L), the product X * A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0D+00

    do j = 1, i
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a(1:m,1:m,i+1-j) ), x(1:m,j) )
    end do

    do j = i + 1, l
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a(1:m,1:m,l+j-i) ), x(1:m,j) )
    end do

  end do

  return
end
subroutine r8bto_mv ( m, l, a, x, b )

!*****************************************************************************80
!
!! R8BTO_MV multiplies an R8BTO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 79, 138 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) X(M*L), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M*L), the product A * X.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0D+00

    do j = 1, i - 1
      b(1:m,i) = b(1:m,i) + matmul ( a(1:m,1:m,l+i-j), x(1:m,j) )
    end do

    do j = i, l
      b(1:m,i) = b(1:m,i) + matmul ( a(1:m,1:m,j+1-i), x(1:m,j) )
    end do

  end do

  return
end
subroutine r8bto_print ( m, l, a, title )

!*****************************************************************************80
!
!! R8BTO_PRINT prints an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  character ( len = * ) title

  call r8bto_print_some ( m, l, a, 1, 1, m * l, m * l, title )

  return
end
subroutine r8bto_print_some ( m, l, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BTO_PRINT_SOME prints some of an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3hi
  integer ( kind = 4 ) i3lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) j3hi
  integer ( kind = 4 ) j3lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) n
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  n = m * l
!
!  Print the columns of the matrix, in strips of 5.
!
  do j3lo = jlo, jhi, incx

    j3hi = j3lo + incx - 1
    j3hi = min ( j3hi, n )
    j3hi = min ( j3hi, jhi )

    inc = j3hi + 1 - j3lo

    write ( *, '(a)' ) ' '

    do j = j3lo, j3hi
      j3 = j + 1 - j3lo
      write ( ctemp(j3), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j3), j3 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i3lo = max ( ilo, 1 )
    i3hi = min ( ihi, n )

    do i = i3lo, i3hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j3 = 1, inc

        j = j3lo - 1 + j3
!
!  i = M * ( i1 - 1 ) + i2
!  j = M * ( j1 - 1 ) + j2
!
        i1 = ( i - 1 ) / m + 1
        i2 = i - m * ( i1 - 1 )
        j1 = ( j - 1 ) / m + 1
        j2 = j - m * ( j1 - 1 )

        if ( i1 <= j1 ) then
          aij = a(i2,j2,j1+1-i1)
        else
          aij = a(i2,j2,l+i1-j1)
        end if

        write ( ctemp(j3), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j3), j3 = 1, inc )

    end do

  end do

  return
end
subroutine r8bto_random ( m, l, seed, a )

!*****************************************************************************80
!
!! R8BTO_RANDOM randomizes an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
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
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A1(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) l

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  do i = 1, m
    do j = 1, m
      do k = 1, 2 * l - 1 
        a(i,j,k) = r8_uniform_01 ( seed )
      end do
    end do
  end do

  return
end
subroutine r8bto_sl ( m, l, a, b, x )

!*****************************************************************************80
!
!! R8BTO_SL solves an R8BTO system.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M*M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) B(M*L), the right hand side vector.
!
!    Output, real ( kind = 8 ) X(M*L), the solution vector.  X and B
!    may share storage.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m*m,2*l-1)
  real ( kind = 8 ) b(m,l)
  real ( kind = 8 ) c1(m,m,l-1)
  real ( kind = 8 ) c2(m,m,l-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pivot(m)
  real ( kind = 8 ) r1(m,m)
  real ( kind = 8 ) r2(m,m)
  real ( kind = 8 ) r3(m,m)
  real ( kind = 8 ) r5(m,m)
  real ( kind = 8 ) r6(m,m)
  real ( kind = 8 ) x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1
  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a(i3,1)
      r1(i,j) = a(i3,1)
      i3 = i3 + 1
    end do
  end do

  r3(1:m,1:m) = r1(1:m,1:m)
  x(1:m,1) = b(1:m,1)

  call r8ge_fa ( m, r3, pivot, info )

  job = 0
  call r8ge_sl ( m, r3, pivot, x(1,1), job )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system
!  with the block Toeplitz matrix for N = 2 through L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    i3 = 1
    do j = 1, m
      do i = 1, m
        r5(i,j) = a(i3,l+n-1)
        r6(i,j) = a(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( 2 < n ) then

      c1(1:m,1:m,n-1) = r2(1:m,1:m)

      do i1 = 1, n - 2
        i2 = n - i1
        do j = 1, m
          i3 = 1
!
!  My apologies, but we have an I1+L, followed in the next line by I1+1.
!
          do i = 1, m
            r5(1:m,j) = r5(1:m,j) + c1(i,j,i2) * a(i3:i3+m-1,i1+l)
            r6(1:m,j) = r6(1:m,j) + c2(i,j,i1) * a(i3:i3+m-1,i1+1)
            i3 = i3 + m
          end do
        end do
      end do

    end if

    do j = 1, m
      r2(1:m,j) = - r5(1:m,j)
      job = 0
      call r8ge_sl ( m, r3, pivot, r2(1,j), job )
    end do

    r3(1:m,1:m) = r6(1:m,1:m)
    r6(1:m,1:m) = - c1(1:m,1:m,1)

    do j = 1, m
      do i = 1, m
        c1(1:m,j,1) = c1(1:m,j,1) + r2(i,j) * r3(1:m,i)
      end do
    end do

    call r8ge_fa ( m, r6, pivot, info )

    do j = 1, m
      call r8ge_sl ( m, r6, pivot, r3(1,j), job )
      do i = 1, m
        r1(1:m,j) = r1(1:m,j) + r3(i,j) * r5(1:m,i)
      end do
    end do

    if ( 2 < n ) then

      r6(1:m,1:m) = c2(1:m,1:m,1)

      do i1 = 2, n - 1

        if ( i1 /= n - 1 ) then
          r5(1:m,1:m) = c2(1:m,1:m,i1)
        end if

        do j = 1, m
          c2(1:m,j,i1) = r6(1:m,j)
          do i = 1, m
            c2(1:m,j,i1) = c2(1:m,j,i1) + r3(i,j) * c1(1:m,i,i1)
          end do
        end do

        do j = 1, m
          do i = 1, m
            c1(1:m,j,i1) = c1(1:m,j,i1) + r2(i,j) * r6(1:m,i)
          end do
        end do

        r6(1:m,1:m) = r5(1:m,1:m)

      end do

    end if

    c2(1:m,1:m,1) = r3(1:m,1:m)
!
!  Compute the solution of the system with the principal minor of order M*N.
!
    r3(1:m,1:m) = r1(1:m,1:m)
    x(1:m,n) = b(1:m,n)

    do i1 = 1, n - 1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        x(1:m,n) = x(1:m,n) - x(i,i2) * a(i3:i3+m-1,i1+l)
        i3 = i3 + m
      end do
    end do

    call r8ge_fa ( m, r3, pivot, info )

    call r8ge_sl ( m, r3, pivot, x(1,n), job )

    do i1 = 1, n - 1
      do i = 1, m
        x(1:m,i1) = x(1:m,i1) + x(i,n) * c2(1:m,i,i1)
      end do
    end do

  end do

  return
end
subroutine r8bto_to_r8ge ( m, l, a, b )

!*****************************************************************************80
!
!! R8BTO_TO_R8GE copies an R8BTO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
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
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the R8BTO matrix.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column
!    of the R8BTO matrix.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Output, real ( kind = 8 ) B(ML,ML), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m*l,m*l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n

  n = m * l

  b(1:n,1:n) = 0.0D+00

  do i = 1, n

    i1 = 1 + ( i - 1 ) / m
    i2 = i - m * ( i1 - 1 )

    do j = 1, n

      j1 = ( j - 1 ) / m + 1
      j2 = j - m * ( j1 - 1 )

      if ( i1 <= j1 ) then
        b(i,j) = a(i2,j2,j1+1-i1)
      else
        b(i,j) = a(i2,j2,l+i1-j1)
      end if

    end do

  end do

  return
end
subroutine r8bto_zeros ( m, l, a )

!*****************************************************************************80
!
!! R8BTO_ZEROS zeroes an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
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
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Output, real ( kind = 8 ) A1(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) l

  real ( kind = 8 ) a(m,m,2*l-1)

  a(1:m,1:m,1:2*l-1) = 0.0D+00

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
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
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8ge_indicator ( m, n, a )

!*****************************************************************************80
!
!! R8GE_INDICATOR sets up an R8GE indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
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
!    11 January 2004
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
!    Output, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      a(i,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

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
subroutine r8ge_sl ( n, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
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
!    R8GE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 1999
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
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

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
subroutine r8vec_zeros ( n, r )

!*****************************************************************************80
!
!! R8VEC_ZEROS zeros an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries 
!    in the vector.
!
!    Output, real ( kind = 8 ) R(N), the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)

  r(1:n) = 0.0D+00

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
