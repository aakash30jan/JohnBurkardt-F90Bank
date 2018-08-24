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
subroutine r8to_dif2 ( n, a )

!*****************************************************************************80
!
!! R8TO_DIF2 sets the second difference as an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a real Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)

  a(1:2*n-1) = 0.0D+00

  a(1) = 2.0D+00
  a(2) = -1.0D+00
  a(n+1) = -1.0D+00

  return
end
subroutine r8to_indicator ( n, a )

!*****************************************************************************80
!
!! R8TO_INDICATOR sets up an R8TO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8TO storage format is used for a real Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1
  k = 0
  do j = 1, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = 8 )
  end do

  j = 1
  do i = 2, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = 8 )
  end do
  
  return
end
subroutine r8to_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! R8TO_MTV multiplies an R8VEC by an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = sum ( a(i:1:-1) * x(1:i) ) + &
           sum ( a(n+1:2*n-i) * x(i+1:n) )

  end do

  return
end
subroutine r8to_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R8TO_MV multiplies an R8TO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  b(1) = sum ( a(1:n) * x(1:n) )

  do i = 2, n
    b(i) = sum ( a(n+i-1:n+1:-1) * x(1:i-1) ) &
         + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8to_print ( n, a, title )

!*****************************************************************************80
!
!! R8TO_PRINT prints an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
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
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  character ( len = * ) title

  call r8to_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8to_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8TO_PRINT_SOME prints some of an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
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
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) aij
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

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8to_random ( n, seed, a )

!*****************************************************************************80
!
!! R8TO_RANDOM randomizes an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
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
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seed

  n2 = 2 * n - 1

  call r8vec_uniform_01 ( n2, seed, a )

  return
end
subroutine r8to_sl ( n, a, b, x )

!*****************************************************************************80
!
!! R8TO_SL solves an R8TO system.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) B(N) the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.  X and B may share the
!    same storage.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c1(n-1)
  real ( kind = 8 ) c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r5
  real ( kind = 8 ) r6
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    r5 = a(n+nsub-1)
    r6 = a(nsub)

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub - 2
        r5 = r5 + a(n+i) * c1(nsub-i)
        r6 = r6 + a(i+1) * c2(i)
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = 0.0D+00

      do i = 2, nsub - 1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = sum ( a(n+1:n+nsub-1) * x(nsub-1:1:-1) )

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine r8to_slt ( n, a, b, x )

!*****************************************************************************80
!
!! R8TO_SLT solves the system A'*x=b, where A is an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) B(N) the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.  X and B may share the
!    same storage.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c1(n-1)
  real ( kind = 8 ) c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r5
  real ( kind = 8 ) r6
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    r5 = a(nsub)
    r6 = a(n+nsub-1)

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub - 2
        r5 = r5 + a(i+1) * c1(nsub-i)
        r6 = r6 + a(n+i) * c2(i)
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = 0.0D+00

      do i = 2, nsub - 1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = sum ( a(2:nsub) * x(nsub-1:1:-1) )

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine r8to_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8TO_TO_R8GE copies an R8TO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
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
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i

  do i = 1, n
    b(i,1:i-1) = a(n+i-1:n+1:-1)
    b(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine r8to_zeros ( n, a )

!*****************************************************************************80
!
!! R8TO_ZEROS zeroes an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
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
!    N must be positive.
!
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) n2

  a(1:2*n-1) = 0.0D+00

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
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
!    Input, integer ( kind = 4 ) N, the number of entries 
!    in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
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
      seed = seed + 2147483647
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
