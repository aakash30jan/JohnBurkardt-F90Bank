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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2000
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
!    Input, character ( len = * ) TITLE, a title first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,i11)' ) i, a(i)
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
subroutine r8gd_dif2 ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8GD_DIF2 sets up an R8GD second difference matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 3.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.  The values -1, 0 and +1 should be included.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdiag
  integer ( kind = 4 ) offset(ndiag)

  if ( ndiag < 3 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8GD_DIF2 - Fatal error!'
    write ( *, '(a)' ) '  NDIAG must be at least 3.'
    stop 1
  end if

  a(1:n,1:ndiag) = 0.0D+00

  do i = 1, n
    do jdiag = 1, ndiag
      j = i + offset(jdiag)
      if ( 1 <= j .and. j <= n ) then
        if ( offset(jdiag) == 0 ) then
          a(i,jdiag) = 2.0D+00
        else if ( offset(jdiag) == -1 .or. offset(jdiag) == +1 ) then
          a(i,jdiag) = -1.0D+00
        end if
      end if
    end do
  end do

  return
end
subroutine r8gd_indicator ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8GD_INDICATOR sets up an R8GD indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)

  a(1:n,1:ndiag) = 0.0D+00

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = real ( fac * i + j, kind = 8 )
      end if
    end do
  end do

  return
end
subroutine r8gd_mtv ( n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8GD_MTV multiplies an R8VEC by an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(j) = b(j) + x(i) * a(i,diag)
      end if
    end do
  end do

  return
end
subroutine r8gd_mv ( n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8GD_MV multiplies an R8GD matrix by an R8VEC.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Example:
!
!    The "offset" value is printed near the first entry of each diagonal
!    in the original matrix, and above the columns in the new matrix.
!
!    Original matrix               New Matrix
!
!      0    1   2   3   4   5        -3  -2   0   1   3   5
!       \    \   \   \   \   \
!        11  12   0  14   0  16      --  --  11  12  14  16
!   -1 =  0  22  23   0  25   0      --  --  22  23  25  --
!   -2 = 31   0  33  34   0  36      --  31  33  34  36  --
!   -3 = 41  42   0  44  45   0      41  42  44  45  --  --
!   -4 =  0  52  53   0  55  56      52  53  55  56  --  --
!   -5 =  0   0  63  64  65  66      63  64  66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(i) = b(i) + a(i,diag) * x(j)
      end if
    end do
  end do

  return
end
subroutine r8gd_print ( n, ndiag, offset, a, title )

!*****************************************************************************80
!
!! R8GD_PRINT prints an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
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
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the 
!    diagonal storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)
  character ( len = * ) title

  call r8gd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

  return
end
subroutine r8gd_print_some ( n, ndiag, offset, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GD_PRINT_SOME prints some of an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) diag
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
  integer ( kind = 4 ) off
  integer ( kind = 4 ) offset(ndiag)
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
        off = j - i
        do diag = 1, ndiag
          if ( off == offset(diag) ) then
            aij = a(i,diag)
          end if
        end do

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8gd_random ( n, ndiag, offset, seed, a )

!*****************************************************************************80
!
!! R8GD_RANDOM randomizes an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) seed

  a(1:n,1:ndiag) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = r8_uniform_01 ( seed )
      end if
    end do
  end do

  return
end
subroutine r8gd_to_r8ge ( n, ndiag, offset, a, b )

!*****************************************************************************80
!
!! R8GD_TO_R8GE copies an R8GD matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
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
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(i,j) = a(i,diag)
      end if
    end do
  end do

  return
end
subroutine r8gd_zeros ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8GD_ZEROS zeroes an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  a(1:n,1:ndiag) = 0.0D+00

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
