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
subroutine r8ltt_det ( n, a, det )

!*****************************************************************************80
!
!! R8LTT_DET computes the determinant of a R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) det

  det = a(1) ** n

  return
end
subroutine r8ltt_indicator ( n, a )

!*****************************************************************************80
!
!! R8LTT_INDICATOR sets up a R8LTT indicator matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) j

  do j = 1, n
    a(j) = real ( j, kind = 8 )
  end do

  return
end
subroutine r8ltt_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8LTT_INVERSE computes the inverse of a R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pn(n)
!
!  Initialize B.
!
  d = 1.0 / a(1)
  b(1) = d
  b(2:n) = 0.0D+00
!
!  Set the strict upper triangle.
!
  p(1) = 0.0D+00
  p(2:n) = a(2:n)
!
!  PN will hold powers of P.
!
  pn(1) = 1.0D+00
  pn(2:n) = 0.0D+00
!
!  Add N-1 powers of strict upper triangle.
!
  do i = 2, n
    d = - d / a(1)
    call r8ltt_mm ( n, p, pn, pn )
    b = b + d  * pn
  end do

  return
end
subroutine r8ltt_mm ( n, a, b, c )

!*****************************************************************************80
!
!! R8LTT_MM computes C = A * B, where A and B are R8LTT matrices.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N), the first factor.
!
!    Input, real ( kind = 8 ) B(N), the second factor.
!
!    Output, real ( kind = 8 ) C(N), the product.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) k

  do k = 1, n
    d(n+1-k) = b(k)
  end do

  call r8ltt_mtv ( n, a, d, e )

  c(1:n) = e(n:1:-1)

  return
end
subroutine r8ltt_mtm ( n, a, b, c )

!*****************************************************************************80
!
!! R8LTT_MTM computes C = A' * B, where A and B are R8LTT matrices.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!    Note that the result C is a dense matrix, of type R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N), B(N), the factors.
!
!    Output, real ( kind = 8 ) C(N,N), the product.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) d(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  d(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n
      do k = max ( i, j ), n
        d(i,j) = d(i,j) + a(k-i+1) * b(k-j+1)
      end do
    end do
  end do

  c(1:n,1:n) = d(1:n,1:n)

  return
end
subroutine r8ltt_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! R8LTT_MTV computes b = A'*x, where A is an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do d = 1, n
    do i = d, n
      j = i + 1 - d
      b(j) = b(j) + a(i-j+1) * x(i)
    end do
  end do

  return
end
subroutine r8ltt_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R8LTT_MV computes b=A*x, where A is an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do d = 1, n
    do i = d, n
      j = i + 1 - d
      b(i) = b(i) + a(i-j+1) * x(j)
    end do
  end do

  return
end
subroutine r8ltt_print ( n, a, title )

!*****************************************************************************80
!
!! R8LTT_PRINT prints a R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  character ( len = * ) title

  call r8ltt_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8ltt_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8LTT_PRINT_SOME prints some of a R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!    1 <= ILO <= IHI <= N.
!    1 <= JLO <= JHI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
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

        if ( i < j ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+1)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ltt_random ( n, seed, a )

!*****************************************************************************80
!
!! R8LTT_RANDOM randomizes an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) A(N), the matrix.
!
!    Output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8ltt_sl ( n, a, b, x )

!*****************************************************************************80
!
!! R8LTT_SL solves a linear system A*x=b with an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!    No factorization of the lower triangular matrix is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  x(1:n) = b(1:n)

  do j = 1, n
    x(j) = x(j) / a(1)
    do i = j + 1, n
      x(i) = x(i) - a(i-j+1) * x(j)
    end do
  end do

  return
end
subroutine r8ltt_slt ( n, a, b, x )

!*****************************************************************************80
!
!! R8LTT_SLT solves a linear system A'*x=b with an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!    No factorization of the lower triangular matrix is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  x(1:n) = b(1:n)

  do i = n, 1, -1
    x(i) = x(i) / a(1)
    do j = 1, i - 1
      x(j) = x(j) - a(i-j+1) * x(i)
    end do
  end do

  return
end
subroutine r8ltt_to_r8ge ( n, a_ltt, a_ge )

!*****************************************************************************80
!
!! R8LTT_TO_R8GE copies an R8LTT matrix to an R8GE matrix.
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
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LTT(N), the matrix.
!
!    Output, real ( kind = 8 ) A_GE(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_ge(n,n)
  real ( kind = 8 ) a_ltt(n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  a_ge(1:n,1:n) = 0.0D+00

  do d = 1, n
    do i = d, n
      j = i + 1 - d
      a_ge(i,j) = a_ltt(i-j+1)
    end do
  end do

  return
end
subroutine r8ltt_zeros ( n, a )

!*****************************************************************************80
!
!! R8LTT_ZEROS zeroes an R8LTT matrix.
!
!  Discussion:
!
!    The R8LTT storage format is used for an N by N lower triangular Toeplitz
!    matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)

  a(1:n) = 0.0D+00

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator vector (1,2,3,...).
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
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array.
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
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
!    13 August 2014
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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
