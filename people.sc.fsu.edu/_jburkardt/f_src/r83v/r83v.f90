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
subroutine i4_log_10_test ( )

!*****************************************************************************80
!
!! I4_LOG_10_TEST tests I4_LOG_10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 13

  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x
  integer ( kind = 4 ), dimension ( test_num ) :: x_test = (/ &
    0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_LOG_10_TEST'
  write ( *, '(a)' ) '  I4_LOG_10: whole part of log base 10,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, I4_LOG_10'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2x, i8, i12 )' ) x, i4_log_10 ( x )
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
subroutine r8_uniform_01_test ( )

!*****************************************************************************80
!
!! R8_UNIFORM_01_TEST tests R8_UNIFORM_01
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8_UNIFORM_01 samples a uniform random'
  write ( *, '(a)' ) '  distribution in [0,1].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Starting with seed = ', seed

  do i = 1, n
    x(i) = r8_uniform_01 ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few values:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  mean = 0.0D+00
  do i = 1, n
    mean = mean + x(i)
  end do
  mean = mean / real ( n, kind = 8 )
 
  variance = 0.0D+00
  do i = 1, n
    variance = variance + ( x(i) - mean ) ** 2
  end do
  variance = variance / real ( n, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,g14.6)' ) '  Minimum value was ', minval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Maximum value was ', maxval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  return
end
subroutine r83_indicator ( m, n, a )

!*****************************************************************************80
!
!! R83_INDICATOR sets up an R83 indicator matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) A(3,N), the R83 indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  a(1:3,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      a(i-j+2,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r83_indicator_test ( )

!*****************************************************************************80
!
!! R83_INDICATOR_TEST tests R83_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83_INDICATOR sets up an R83 indicator matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:3,1:n) )

    call r83_indicator ( m, n, a )

    call r83_print ( m, n, a, '  The R83 indicator matrix:' )

    deallocate ( a )

  end do

  return
end
subroutine r83_print ( m, n, a, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) m
  character ( len = * ) title

  call r83_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r83_print_test ( )

!*****************************************************************************80
!
!! R83_PRINT_TEST tests R83_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_PRINT_TEST'
  write ( *, '(a)' ) '  R83_PRINT prints an R83 matrix.'
!
!  Set the matrix.
!
  call r83_indicator ( m, n, a )

  call r83_print ( m, n, a, '  The R83 matrix:' )

  return
end
subroutine r83_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83_PRINT_SOME prints some of an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
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
  integer ( kind = 4 ) m
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
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
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i - j + 2 < 1 .or. 3 < i - j + 2 ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+2,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83_print_some_test ( )

!*****************************************************************************80
!
!! R83_PRINT_SOME_TEST tests R83_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83_PRINT_SOME prints some of an R83 matrix.'
!
!  Set the matrix.
!
  call r83_indicator ( m, n, a )

  call r83_print_some ( m, n, a, 2, 2, 5, 4, '  Rows 2-5, Cols 2-4:' )

  return
end
subroutine r83v_cg ( n, a, b, c, ax, x )

!*****************************************************************************80
!
!! R83V_CG uses the conjugate gradient method on an R83V system.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frank Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    in Mathematical Methods for Digital Computers,
!    edited by John Ralston, Herbert Wilf,
!    Wiley, 1967,
!    ISBN: 0471706892,
!    LC: QA76.5.R3.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N-1)), B(N), C(N-1), 
!    the R83V matrix.
!
!    Input, real ( kind = 8 ) AX(N), the right hand side vector.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n-1)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ap(n)
  real ( kind = 8 ) ax(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(n-1)
  integer ( kind = 4 ) it
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pap
  real ( kind = 8 ) pr
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rap
  real ( kind = 8 ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call r83v_mv ( n, n, a, b, c, x, ap )

  r(1:n) = ax(1:n) - ap(1:n)
  p(1:n) = ax(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r83v_mv ( n, n, a, b, c, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = dot_product ( p(1:n), ap(1:n) )
    pr = dot_product ( p(1:n), r(1:n) )

    if ( pap == 0.0D+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = dot_product ( r(1:n), ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine r83v_cg_test ( )

!*****************************************************************************80
!
!! R83V_CG_TEST tests R83V_CG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ) e_norm
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ) r_norm
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CG_TEST'
  write ( *, '(a)' ) '  R83V_CG applies CG to an R83V matrix.'
!
!  Set A to the second difference matrix.
!
  n = 10

  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )
  call r8vec_print ( n, x1, '  Desired solution X1:' )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00

  call r83v_cg ( n, a, b, c, ax, x2 )
  call r8vec_print ( n, x2, '  Computed solution X2' )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83v_res ( n, n, a, b, c, x2, ax, r )

  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_copy ( m, n, a1, a2, a3, b1, b2, b3 )

!*****************************************************************************80
!
!! R83V_COPY copies an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), 
!    the R83V matrix.
!
!    Output, real ( kind = 8 ) B1(min(M-1,N)), B2(min(M,N)), B3(min(M,N-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(min(m-1,n))
  real ( kind = 8 ) a2(min(m,n))
  real ( kind = 8 ) a3(min(m,n-1))
  real ( kind = 8 ) b1(min(m-1,n))
  real ( kind = 8 ) b2(min(m,n))
  real ( kind = 8 ) b3(min(m,n-1))

  b1 = a1(1:min(m-1,n  ))
  b2 = a2(1:min(m,  n  ))
  b3 = a3(1:min(m,  n-1))
 
  return
end
subroutine r83v_copy_test ( )

!*****************************************************************************80
!
!! R83V_COPY_TEST tests R83V_COPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  real ( kind = 8 ), allocatable :: b1(:)
  real ( kind = 8 ), allocatable :: b2(:)
  real ( kind = 8 ), allocatable :: b3(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_COPY_TEST'
  write ( *, '(a)' ) '  R83V_COPY copies an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  A:' )

    allocate ( b1(1:min(m-1,n)) )
    allocate ( b2(1:min(m,n)) )
    allocate ( b3(1:min(m,n-1)) )

    call r83v_copy ( m, n, a1, a2, a3, b1, b2, b3 )
    call r83v_print ( m, n, b1, b2, b3, '  B = copy of A:' )

    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )
    deallocate ( b1 )
    deallocate ( b2 )
    deallocate ( b3 )

  end do

  return
end
subroutine r83v_cr_fa ( n, a, b, c, a_cr )

!*****************************************************************************80
!
!! R83V_CR_FA decomposes an R83V matrix using cyclic reduction.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used to
!    solve linear systems A * x = b.
!
!    R83_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, R83_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause R83_CR_FA to fail.
!
!    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 R83 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             R83    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             R83   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by R83_CR_SL to solve linear
!    systems.
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2015
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)),
!    the R83V matrix.
!
!    Output, real ( kind = 8 ) A_CR(3,2*N+1), factorization information.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n-1)
  real ( kind = 8 ) a_cr(3,2*n+1)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n-1)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulp
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ilp
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp

  if ( n <= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R83V_CR_FA - Fatal error!'
    write ( *, '(a,i4)' ) '  Nonpositive N = ', n
    stop 1
  end if

  a_cr(1:3,1:2*n+1) = 0.0D+00

  if ( n == 1 ) then
    a_cr(2,2) = 1.0D+00 / b(1)
    return
  end if
!
!  Set the workspace entries.
!
  a_cr(1,2:n)   = c(1:n-1)
  a_cr(2,2:n+1) = b(1:n)
  a_cr(3,2:n)   = a(1:n-1)

  il = n
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod ( il, 2 ) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = ( inc / 2 )
    il = ( il / 2 )
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      a_cr(2,iful+1) = 1.0D+00 / a_cr(2,iful+1)
      a_cr(3,iful+1)  = a_cr(3,iful+1)  * a_cr(2,iful+1)
      a_cr(1,ifulp+1) = a_cr(1,ifulp+1) * a_cr(2,ifulp+2)
      a_cr(2,ihaf+1)  = a_cr(2,ifulp+1) - a_cr(1,iful+1)  * a_cr(3,iful+1) &
                                        - a_cr(1,ifulp+1) * a_cr(3,ifulp+1)
      a_cr(3,ihaf+1) = -a_cr(3,ifulp+1) * a_cr(3,ifulp+2)
      a_cr(1,ihaf+1) = -a_cr(1,ifulp+1) * a_cr(1,ifulp+2)
    end do

  end do

  a_cr(2,ipntp+2) = 1.0D+00 / a_cr(2,ipntp+2)

  return
end
subroutine r83v_cr_fa_test ( )

!*****************************************************************************80
!
!! R83V_CR_FA_TEST tests R83V_CR_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: a_cr(:,:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  logical ( kind = 4 ) debug
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x(:)

  n = 5
  debug = .false.

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_FA_TEST'
  write ( *, '(a)' ) '  R83V_CR_FA factors an R83V matrix'
  write ( *, '(a)' ) '  Once the matrix has been factored, we can call'
  write ( *, '(a)' ) '  R83V_CR_SL to solve a linear system.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  System matrix A:' )
!
!  Factor the matrix once.
!
  allocate ( a_cr(1:3,1:2*n+1) )
  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Print the factor information.
!
  if ( debug ) then
    call r83_print ( n, 2 * n + 1, a_cr, &
      '  Cyclic reduction factor information:' )
  end if

  allocate ( ax(1:n) )
  allocate ( x(1:n) )

  do j = 1, 2

    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Solve linear system number #', j

    if ( j == 1 ) then
      ax(1:n-1) = 0.0D+00
      ax(n) = real ( n + 1, kind = 8 )
    else
      ax(1) = 1.0D+00
      ax(2:n-1) = 0.0D+00
      ax(n) = 1.0D+00
    end if
!
!  Solve the linear system.
!
    call r83v_cr_sl ( n, a_cr, ax, x )

    call r8vec_print ( n, x, '  Solution:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_cr_sl ( n, a_cr, ax, x )

!*****************************************************************************80
!
!! R83V_CR_SL solves a real linear system factored by R83V_CR_FA.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or many
!    linear systems involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SL may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_CR(3,2*N+1), factorization information computed 
!    by R83V_CR_FA.
!
!    Input, real ( kind = 8 ) AX(N), the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_cr(3,2*n+1)
  real ( kind = 8 ) ax(n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  real ( kind = 8 ) rhs(2*n+1)
  real ( kind = 8 ) x(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R83V_CR_SL - Fatal error!'
    write ( *, '(a,i4)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = a_cr(2,2) * ax(1)
    return
  end if
!
!  Set up RHS.
!
  rhs(1) = 0.0D+00
  rhs(2:n+1) = ax(1:n)
  rhs(n+2:2*n+1) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = ( il / 2 )
    ndiv = ndiv * 2
    ihaf = ipntp

    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf+1) = rhs(iful+1) - a_cr(3,iful) * rhs(iful) &
        - a_cr(1,iful+1) * rhs(iful+2)
    end do

  end do

  rhs(ihaf+1) = rhs(ihaf+1) * a_cr(2,ihaf+1)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ( ndiv / 2 )
    il = ( n / ndiv )
    ipnt = ipnt - il
    ihaf = ipntp

    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful+1) = rhs(ihaf+1)
      rhs(ifulm+1) = a_cr(2,ifulm+1) * ( rhs(ifulm+1) & 
        - a_cr(3,ifulm) * rhs(ifulm) &
        - a_cr(1,ifulm+1) * rhs(iful+1) )
    end do

  end do

  x(1:n) = rhs(2:n+1)

  return
end
subroutine r83v_cr_sl_test ( )

!*****************************************************************************80
!
!! R83V_CR_SL_TEST tests R83V_CR_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: a_cr(:,:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x(:)

  n = 5

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_SL_TEST'
  write ( *, '(a)' ) '  R83V_CR_SL solves a factored system'
  write ( *, '(a)' ) '  after R83V_CR_FA has factored it.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  a(1:n-1) = 1.0D+00
  b(1:n) = 4.0D+00
  c(1:n-1) = 2.0D+00
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  Input matrix A:' )
!
!  Set the desired solution.
!
  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x,ax )

  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  allocate ( a_cr(1:3,1:2*n+1) )
  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Solve the linear system.
!
  call r83v_cr_sl ( n, a_cr, ax, x )

  call r8vec_print ( n, x, '  Solution:' )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_cr_sls ( n, a_cr, nb, ax, x )

!*****************************************************************************80
!
!! R83V_CR_SLS solves several real linear systems factored by R83V_CR_FA.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)). 
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
!    many linear systems involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SLS may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_CR(3,2*N+1), factorization information computed 
!    by R83V_CR_FA.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) AX(N,NB), the right hand side vectors.
!
!    Output, real ( kind = 8 ) X(N,NB), the solutions of the linear system.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a_cr(3,2*n+1)
  real ( kind = 8 ) ax(n,nb)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  real ( kind = 8 ) rhs(2*n+1,nb)
  real ( kind = 8 ) x(n,nb)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R83V_CR_SLS - Fatal error!'
    write ( *, '(a,i4)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    x(1,1:nb) = a_cr(2,2) * ax(1,1:nb)
    return
  end if
!
!  Set up RHS.
!
  rhs(1,1:nb) = 0.0D+00
  rhs(2:n+1,1:nb) = ax(1:n,1:nb)
  rhs(n+2:2*n+1,1:nb) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = ( il / 2 )
    ndiv = ndiv * 2
    ihaf = ipntp

    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf+1,1:nb) = rhs(iful+1,1:nb) &
        - a_cr(3,iful) * rhs(iful,1:nb) &
        - a_cr(1,iful+1) * rhs(iful+2,1:nb)
    end do

  end do

  rhs(ihaf+1,1:nb) = rhs(ihaf+1,1:nb) * a_cr(2,ihaf+1)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ( ndiv / 2 )
    il = ( n / ndiv )
    ipnt = ipnt - il
    ihaf = ipntp

    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful+1,1:nb) = rhs(ihaf+1,1:nb)
      rhs(ifulm+1,1:nb) = a_cr(2,ifulm+1) * ( rhs(ifulm+1,1:nb) &
        - a_cr(3,ifulm) * rhs(ifulm,1:nb) &
        - a_cr(1,ifulm+1) * rhs(iful+1,1:nb) )
    end do

  end do

  x(1:n,1:nb) = rhs(2:n+1,1:nb)

  return
end
subroutine r83v_cr_sls_test ( )

!*****************************************************************************80
!
!! R83V_CR_SLS_TEST tests R83V_CR_SLS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: a_cr(:,:)
  real ( kind = 8 ), allocatable :: ax(:,:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb
  real ( kind = 8 ), allocatable :: x(:,:)

  n = 5

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_SLS_TEST'
  write ( *, '(a)' ) '  R83V_CR_SLS solves multiple systems A*x1:xn=b1:bn'
  write ( *, '(a)' ) '  after R83V_CR_FA has factored it.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  Input matrix A:' )
!
!  Factor the matrix once.
!
  allocate ( a_cr(3,2*n+1) )

  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Set the number of right hand sides.
!
  nb = 2
  allocate ( ax(1:n,1:nb) )

  ax(1:n-1,1) = 0.0D+00
  ax(n,1) = n + 1

  ax(1,2) = 1.0D+00
  ax(2:n-1,2) = 0.0D+00
  ax(n,2) = 1.0D+00

  call r8ge_print ( n, 2, ax, '  Right hand sides b1:b2' )
!
!  Solve the linear systems.
!
  allocate ( x(1:n,1:nb) )

  call r83v_cr_sls ( n, a_cr, nb, ax, x )

  call r8ge_print ( n, nb, x, '  Solutions x1:x2' )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_dif2 ( m, n, a, b, c )

!*****************************************************************************80
!
!! R83V_DIF2 returns the DIF2 matrix in R83V format.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Properties:
!
!    A is banded, with bandwidth 3.
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is a special case of the TRIS or tridiagonal scalar matrix.
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is Toeplitz: constant along diagonals.
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is positive definite.
!
!    A is an M matrix.
!
!    A is weakly diagonally dominant, but not strictly diagonally dominant.
!
!    A has an LU factorization A = L * U, without pivoting.
!
!      The matrix L is lower bidiagonal with subdiagonal elements:
!
!        L(I+1,I) = -I/(I+1)
!
!      The matrix U is upper bidiagonal, with diagonal elements
!
!        U(I,I) = (I+1)/I
!
!      and superdiagonal elements which are all -1.
!
!    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
!
!      L(I,I) =    sqrt ( (I+1) / I )
!      L(I,I-1) = -sqrt ( (I-1) / I )
!
!    The eigenvalues are
!
!      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
!                = 4 SIN^2(I*PI/(2*N+2))
!
!    The corresponding eigenvector X(I) has entries
!
!       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
!
!    Simple linear systems:
!
!      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
!
!      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
!
!    det ( A ) = N + 1.
!
!    The value of the determinant can be seen by induction,
!    and expanding the determinant across the first row:
!
!      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
!                = 2 * N - (N-1)
!                = N + 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969,
!    ISBN: 0882756494,
!    LC: QA263.68
!
!    Morris Newman, John Todd,
!    Example A8,
!    The evaluation of matrix inversion programs,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 6, Number 4, pages 466-476, 1958.
!
!    John Todd,
!    Basic Numerical Mathematics,
!    Volume 2: Numerical Algebra,
!    Birkhauser, 1980,
!    ISBN: 0817608117,
!    LC: QA297.T58.
!
!    Joan Westlake,
!    A Handbook of Numerical Matrix Inversion and Solution of 
!    Linear Equations,
!    John Wiley, 1968,
!    ISBN13: 978-0471936756,
!    LC: QA263.W47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))

  a(1:min(m-1,n)  ) = -1.0D+00
  b(1:min(m,  n)  ) = +2.0D+00
  c(1:min(m,  n-1)) = -1.0D+00

  return
end
subroutine r83v_dif2_test ( )

!*****************************************************************************80
!
!! R83V_DIF2_TEST tests R83V_DIF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_DIF2_TEST'
  write ( *, '(a)' ) '  R83V_DIF2 sets an R83V matrix to the second difference.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_dif2 ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Second difference in R83V format:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_fs ( n, a1, a2, a3, b, x )

!*****************************************************************************80
!
!! R83V_FS solves a linear system with R83V matrix.
!
!  Discussion:
!
!    This function is based on the LINPACK SGTSL routine.
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt, based on the LINPACK SGTSL function.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the tridiagonal matrix.
!
!    Input, real ( kind = 8 ) A1(N-1), A2(N), A3(N-1), the R83V matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n-1)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n-1)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
!
!  Copy the input data.
!
  c(1) = 0.0D+00
  c(2:n) = a1(1:n-1)
  d(1:n) = a2(1:n)
  e(1:n-1) = a3(1:n-1)
  e(n) = 0.0D+00

  x(1:n) = b(1:n)
!
!  Factor.
!
  c(1) = a2(1)

  if ( 2 <= n ) then

    d(1) = e(1)
    e(1) = 0.0D+00
    e(n) = 0.0D+00

    do k = 2, n
!
!  Find the larger of the two rows and interchange if necessary.
!
      if ( abs ( c(k-1) ) <= abs ( c(k) ) ) then

        t = c(k)
        c(k) = c(k-1)
        c(k-1) = t

        t = d(k)
        d(k) = d(k-1)
        d(k-1) = t

        t = e(k)
        e(k) = e(k-1)
        e(k-1) = t

        t = x(k)
        x(k) = x(k-1)
        x(k-1) = t

      end if
!
!  Zero elements.
!
      if ( c(k-1) == 0.0D+00 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'R83V_FS - Fatal error!'
        write ( *, '(a,i4)' ) '  Zero pivot on step K = ', k
        stop 1
      end if

      t = - c(k) / c(k-1)
      c(k) = d(k) + t * d(k-1)
      d(k) = e(k) + t * e(k-1)
      e(k) = 0.0D+00
      x(k) = x(k) + t * x(k-1)

    end do

  end if

  if ( c(n) == 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R83V_FS - Fatal error!'
    write ( *, '(a,i4)' ) '  Zero pivot on step K = ', n
    stop 1
  end if
!
!  Back solve.
!
  x(n) = x(n) / c(n)

  if ( 1 < n ) then

    x(n-1) = ( x(n-1) - d(n-1) * x(n) ) / c(n-1)

    do k = n - 2, 1, -1
      x(k) = ( x(k) - d(k) * x(k+1) - e(k) * x(k+2) ) / c(k)
    end do

  end if

  return
end
subroutine r83v_fs_test ( )

!*****************************************************************************80
!
!! R83V_FS_TEST tests R83V_FS_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  n = 10

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_FS_TEST'
  write ( *, '(a)' ) '  R83V_FS factors and solves a linear system'
  write ( *, '(a)' ) '  for an R83V matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
!
!  Set the matrix values.
!
  allocate ( a1(1:n-1) )
  allocate ( a2(n) )
  allocate ( a3(n-1) )

  call r83v_dif2 ( n, n, a1, a2, a3 )
!
!  Set the desired solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83v_mv ( n, n, a1, a2, a3, x1, b )

  call r8vec_print  ( n, b, "  The right hand side:" )
!
!  Solve the linear system.
!
  allocate ( x2(1:n) )
  call r83v_fs ( n, a1, a2, a3, b, x2 )

  call r8vec_print ( n, x2, "  Solution:" )
!
!  Free memory.
!
  deallocate ( a1 )
  deallocate ( a2 )
  deallocate ( a3 )
  deallocate ( b )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_gs_sl ( n, a, b, c, ax, x, it_max )

!*****************************************************************************80
!
!! R83V_GS_SL solves a R83V system using Gauss-Seidel iteration.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
!    Input, real ( kind = 8 ) AX(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), on input,
!    an approximate solution to the system.  On output, an improved solution
!    estimate.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n-1)
  real ( kind = 8 ) ax(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_new(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( b(i) == 0.0D+00 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'R83V_GS_SL - Fatal error!'
      write ( *, '(a,i4)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  do it_num = 1, it_max

    x_new(1) = ( ax(1) - c(1) * x(2) ) / b(1)
    do i = 2, n - 1
      x_new(i) = ( ax(i) - a(i-1) * x_new(i-1) - c(i) * x(i+1) ) / b(i)
    end do
    x_new(n) = ( ax(n) - a(n-1) * x_new(n-1) ) / b(n)

    x(1:n) = x_new(1:n)

  end do

  return
end
subroutine r83v_gs_sl_test ( )

!*****************************************************************************80
!
!! R83V_GS_SL_TEST tests R83V_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_GS_SL_TEST'
  write ( *, '(a)' ) '  R83V_GS_SL applies Gauss-Seidel iteration with an'
  write ( *, '(a)' ) '  R83V matrix to solve a linear system A*x=b.'
!
!  Set A to the second difference matrix.
!
  n = 10
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Set the starting solution.
!
  allocate ( x2(1:n) )
  x2(1:n) = 0.0D+00
!
!  Call the Gauss-Seidel routine.
!
  it_max = 25

  do i = 1, 3

    call r83v_gs_sl ( n, a, b, c, ax, x2, it_max )

    call r8vec_print ( n, x2, '  Current solution estimate:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_indicator ( m, n, a, b, c )

!*****************************************************************************80
!
!! R83V_INDICATOR sets up an R83V indicator matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, min ( m - 1, n )
    a(i) = real ( fac * ( i + 1 ) + i, kind = 8 )
  end do

  do i = 1, min ( m, n )
    b(i) = real ( fac * i + i, kind = 8 )
  end do

  do i = 1, min ( m, n - 1 )
    c(i) = real ( fac * i + i + 1, kind = 8 )
  end do

  return
end
subroutine r83v_indicator_test ( )

!*****************************************************************************80
!
!! R83V_INDICATOR_TEST tests R83V_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83V_INDICATOR sets an R83V indicator matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  R83V indicator matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_jac_sl ( n, a, b, c, ax, x, it_max )

!*****************************************************************************80
!
!! R83V_JAC_SL solves a R83V system A*x=b using Jacobi iteration.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the 
!    R83V matrix.
!
!    Input, real ( kind = 8 ) AX(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution to the 
!    system, which is updated on output.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n-1)
  real ( kind = 8 ) ax(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_new(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( b(i) == 0.0D+00 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'R83V_JAC_SL - Fatal error!'
      write ( *, '(a,i4)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do
!
!  Iterate IT_MAX times.
!
  do it_num = 1, it_max

    x_new(1:n)   = ax(1:n)
    x_new(1:n-1) = x_new(1:n-1) - c(1:n-1) * x(2:n)
    x_new(2:n)   = x_new(2:n)   - a(1:n-1) * x(1:n-1)
!
!  Divide by diagonal terms.
!
    x_new(1:n) = x_new(1:n) / b(1:n)
!
!  Update.
!
    x(1:n) = x_new(1:n)

  end do

  return
end
subroutine r83v_jac_sl_test ( )

!*****************************************************************************80
!
!! R83V_JAC_SL_TEST tests R83V_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83V_JAC_SL applies Jacobi iteration with an R83V'
  write ( *, '(a)' ) '  matrix to solve a linear system A*x=b.'
!
!  Set A to the second difference matrix.
!
  n = 10
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Set the starting vector.
!
  allocate ( x2(1:n) )
  x2(1:n) = 0.0D+00
!
!  Call the Jacobi routine.
!
  it_max = 25

  do i = 1, 3

    call r83v_jac_sl ( n, a, b, c, ax, x2, it_max )

    call r8vec_print ( n, x2, '  Current solution:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_mtv ( m, n, a, b, c, x, ax )

!*****************************************************************************80
!
!! R83V_MTV multiplies a vector by an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)),
!    the R83V matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) AX(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  integer ( kind = 4 ) ahi
  real ( kind = 8 ) ax(n)
  real ( kind = 8 ) b(min(m,n))
  integer ( kind = 4 ) bhi
  real ( kind = 8 ) c(min(m,n-1))
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(m)
!
!  Find each nonzero A(I,J), multiply by X(I), add to B(J).
!
!  A(K) = A(K+1,K) = A'(K,K+1)
!
  ahi = min ( m - 1, n )
  bhi = min ( m, n )
  chi = min ( m, n - 1 )

  ax(1:n) = 0.0D+00

  ax(1:ahi)   = ax(1:ahi)   + a(1:ahi) * x(2:ahi+1)
  ax(1:bhi)   = ax(1:bhi)   + b(1:bhi) * x(1:bhi)
  ax(2:chi+1) = ax(2:chi+1) + c(1:chi) * x(1:chi)

  return
end
subroutine r83v_mtv_test ( )

!*****************************************************************************80
!
!! R83V_MTV_TEST tests R83V_MTV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a_83v(:)
  real ( kind = 8 ), allocatable :: a_ge(:,:)
  real ( kind = 8 ), allocatable :: ax_83v(:)
  real ( kind = 8 ), allocatable :: ax_ge(:)
  real ( kind = 8 ), allocatable :: b_83v(:)
  real ( kind = 8 ), allocatable :: c_83v(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_MTV_TEST'
  write ( *, '(a)' ) '  R83V_MV computes b=A''*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )

    allocate ( x(1:m) )
    call r8vec_indicator1 ( m, x )

    allocate ( ax_83v(1:n) )
    call r83v_mtv ( m, n, a_83v, b_83v, c_83v, x, ax_83v )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )

    allocate ( ax_ge(1:n) )
    call r8ge_mtv ( m, n, a_ge, x, ax_ge )

    call r8vec2_print ( n, ax_83v, ax_ge, '  Product comparison:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( ax_83v )
    deallocate ( ax_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )
    deallocate ( x )

  end do

  return
end
subroutine r83v_mv ( m, n, a, b, c, x, ax )

!*****************************************************************************80
!
!! R83V_MV multiplies a R83V matrix times a vector.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) AX(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  integer ( kind = 4 ) ahi
  real ( kind = 8 ) ax(m)
  real ( kind = 8 ) b(min(m,n))
  integer ( kind = 4 ) bhi
  real ( kind = 8 ) c(min(m,n-1))
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  ahi = min ( m - 1, n )
  bhi = min ( m, n )
  chi = min ( m, n - 1 )

  ax(1:m) = 0.0D+00
  ax(2:ahi+1) = ax(2:ahi+1) + a(1:ahi) * x(1:ahi)
  ax(1:bhi) = ax(1:bhi) + b(1:bhi) * x(1:bhi)
  ax(1:chi) = ax(1:chi) + c(1:chi) * x(2:chi+1)

  return
end
subroutine r83v_mv_test ( )

!*****************************************************************************80
!
!! R83V_MV_TEST tests R83V_MV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a_83v(:)
  real ( kind = 8 ), allocatable :: a_ge (:,:)
  real ( kind = 8 ), allocatable :: ax_83v(:)
  real ( kind = 8 ), allocatable :: ax_ge(:)
  real ( kind = 8 ), allocatable :: b_83v(:)
  real ( kind = 8 ), allocatable :: c_83v(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_MV_TEST'
  write ( *, '(a)' ) '  R83V_MV computes b=A*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )

    allocate ( x(1:n) )
    call r8vec_indicator1 ( n, x )

    allocate ( ax_83v(1:m) )
    call r83v_mv ( m, n, a_83v, b_83v, c_83v, x, ax_83v )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )

    allocate ( ax_ge(1:m) )
    call r8ge_mv ( m, n, a_ge, x, ax_ge )

    call r8vec2_print ( m, ax_83v, ax_ge, '  Product comparison:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( ax_83v )
    deallocate ( ax_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )
    deallocate ( x )

  end do

  return
end
subroutine r83v_print ( m, n, a, b, c, title )

!*****************************************************************************80
!
!! R83V_PRINT prints an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))
  character ( len = * ) title

  call r83v_print_some ( m, n, a, b, c, 1, 1, m, n, title )

  return
end
subroutine r83v_print_test ( )

!*****************************************************************************80
!
!! R83V_PRINT_TEST tests R83V_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRINT_TEST'
  write ( *, '(a)' ) '  R83V_PRINT prints an R83V matrix.'

  m = 5
  n = 5

  allocate ( a(1:min(m-1,n)) )
  allocate ( b(1:min(m,n)) )
  allocate ( c(1:min(m,n-1)) )

  call r83v_indicator ( m, n, a, b, c )
  call r83v_print ( m, n, a, b, c, '  R83V matrix:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine r83v_print_some ( m, n, a, b, c, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83V_PRINT_SOME prints some of a R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
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
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M+1,N)), 
!    the R83V matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) trim( title )
  end if

  incx = 5
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ''
    write ( *, '(a)', advance = 'no' ) '  Col: '
    do j = j2lo, j2hi
      write ( *, '(i7,7x)', advance = 'no' ) j
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      write ( *, '(i5,a)', advance = 'no' ) i, ':'

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j < i - 1 ) then
          write ( *, '(a)', advance = 'no' ) '              '
        else if ( j == i - 1 ) then
          write ( *, '(g14.6)', advance = 'no' ) a(i-1)
        else if ( j == i ) then
          write ( *, '(g14.6)', advance = 'no' ) b(i)
        else if ( j == i + 1  ) then
          write ( *, '(g14.6)', advance = 'no' ) c(i)
        else

        end if

      end do

      write ( *, '(a)' ) ''

    end do

  end do

  return
end
subroutine r83v_print_some_test ( )

!*****************************************************************************80
!
!! R83V_PRINT_SOME_TEST tests R83V_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83V_PRINT_SOME prints some of an R83V matrix.'

  m = 5
  n = 5

  allocate ( a(1:min(m-1,n)) )
  allocate ( b(1:min(m,n)) )
  allocate ( c(1:min(m,n-1)) )

  call r83v_indicator ( m, n, a, b, c )
  call r83v_print_some ( m, n, a, b, c, 2, 2, 5, 4, '  Rows 2-5, Cols 2-4:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine r83v_random ( m, n, seed, a, b, c )

!*****************************************************************************80
!
!! R83V_RANDOM randomizes an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)),
!    the R83V matrix.
!
!    Output, integer ( kind = 4 ) SEED, an updated seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( min ( m - 1, n     ), seed, a )
  call r8vec_uniform_01 ( min ( m,     n     ), seed, b )
  call r8vec_uniform_01 ( min ( m,     n - 1 ), seed, c )

  return
end
subroutine r83v_random_test ( )

!*****************************************************************************80
!
!! R83V_RANDOM_TEST tests R83V_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_RANDOM_TEST'
  write ( *, '(a)' ) '  R83V_RANDOM randomizes an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Random R83V matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_res ( m, n, a, b, c, x, ax, r )

!*****************************************************************************80
!
!! R83V_RES computes the residual R = b-A*x for R83V matrices.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2016
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
!    Input, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Input, real ( kind = 8 ) AX(M), the desired result A * x.
!
!    Output, real ( kind = 8 ) R(M), the residual R = b - A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) ax(m)
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) x(n)

  call r83v_mv ( m, n, a, b, c, x, r )

  r(1:m) = ax(1:m) - r(1:m)

  return
end
subroutine r83v_res_test ( )

!*****************************************************************************80
!
!! R83V_RES_TEST tests R83V_RES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: ax(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_RES_TEST'
  write ( *, '(a)' ) '  R83V_RES computes b-A*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    allocate ( ax(1:m) )
    allocate ( r(1:m) )
    allocate ( x(1:n) )

    seed = 123456789
    call r83v_random ( m, n, seed, a, b, c )
    call r8vec_indicator1 ( n, x )
    call r83v_mv ( m, n, a, b, c, x, ax )
    call r83v_res ( m, n, a, b, c, x, ax, r )
    call r8vec_print ( m, r, '  Residual A*x-b:' )

    deallocate ( a )
    deallocate ( ax )
    deallocate ( b )
    deallocate ( c )
    deallocate ( r )
    deallocate ( x )

  end do

  return
end
subroutine r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )

!*****************************************************************************80
!
!! R83V_TO_R8GE copies an R83V matrix to a R8GE matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_83V(min(M-1,N)), B_83V(min(M,N)), 
!    C_83V(min(M,N-1)), the R83V matrix.
!
!    Output, real ( kind = 8 ) A_GE(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_83v(min(m-1,n))
  real ( kind = 8 ) a_ge(m,n)
  integer ( kind = 4 ) ahi
  real ( kind = 8 ) b_83v(min(m,n))
  integer ( kind = 4 ) bhi
  real ( kind = 8 ) c_83v(min(m,n-1))
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) k

  a_ge(1:m,1:n) = 0.0D+00

  ahi = min ( m - 1, n )
  do k = 1, ahi
    a_ge(k+1,k) = a_83v(k)
  end do

  bhi = min ( m, n )
  do k = 1, bhi
    a_ge(k,k) = b_83v(k)
  end do

  chi = min ( m, n - 1 )
  do k = 1, chi
    a_ge(k,k+1) = c_83v(k)
  end do

  return
end
subroutine r83v_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83V_TO_R8GE_TEST tests R83V_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a_83v(:)
  real ( kind = 8 ), allocatable :: a_ge(:,:)
  real ( kind = 8 ), allocatable :: b_83v(:)
  real ( kind = 8 ), allocatable :: c_83v(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83V_TO_R8GE converts an R83V matrix to R8GE format.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )
    call r83v_print ( m, n, a_83v, b_83v, c_83v, '  R83V matrix:' )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )
    call r8ge_print ( m, n, a_ge, '  R8GE matrix:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )

  end do

  return
end
subroutine r83v_to_r8vec ( m, n, a1, a2, a3, a )

!*****************************************************************************80
!
!! R83V_TO_R8VEC copies an R83V matrix to an R8VEC.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), 
!    the matrix.
!
!    Output, real ( kind = 8 ) R83V_TO_R8VEC(min(M-1,N)+min(M,N)+min(M,N-1)),
!    the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n)+min(m,n)+min(m,n-1))
  real ( kind = 8 ) a1(min(m-1,n))
  real ( kind = 8 ) a2(min(m,n))
  real ( kind = 8 ) a3(min(m,n-1))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = 0
  do j = 1, n
    if ( j <= m + 1 .and. 2 <= j ) then
      k = k + 1
      a(k) = a3(j-1)
    end if
    if ( j <= m ) then
      k = k + 1
      a(k) = a2(j)
    end if
    if ( j <= m - 1 ) then
      k = k + 1
      a(k) = a1(j)
    end if
  end do

  return
end
subroutine r83v_to_r8vec_test ( )

!*****************************************************************************80
!
!! R83V_TO_R8VEC_TEST tests R83V_TO_R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TO_R8VEC_TEST'
  write ( *, '(a)' ) '  R83V_TO_R8VEC copies an R83V matrix to an R8VEC.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(min(m-1,n)) )
    allocate ( a2(min(m,  n)) )
    allocate ( a3(min(m,n-1)) )
    na = min ( m - 1, n ) + min ( m, n ) + min ( m, n - 1 )
    allocate ( a(na) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, "  R83V matrix A:" )

    call r83v_to_r8vec ( m, n, a1, a2, a3, a )
    call r8vec_print ( na, a, "  Vector version of A:" )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

  end do

  return
end
subroutine r83v_transpose ( m, n, a1, a2, a3, b1, b2, b3 )

!*****************************************************************************80
!
!! R83V_TRANSPOSE makes a transposed copy of an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), 
!    the R83V matrix.
!
!    Output, real ( kind = 8 ) B1(min(N-1,M)), B2(min(N,M)), B3(min(N,M-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(min(m-1,n))
  real ( kind = 8 ) a2(min(m,n))
  real ( kind = 8 ) a3(min(m,n-1))
  real ( kind = 8 ) b1(min(n-1,m))
  real ( kind = 8 ) b2(min(n,m))
  real ( kind = 8 ) b3(min(n,m-1))

  b1 = a3(1:min(m,  n-1))
  b2 = a2(1:min(m,  n  ))
  b3 = a1(1:min(m-1,n  ))
 
  return
end
subroutine r83v_transpose_test ( )

!*****************************************************************************80
!
!! R83V_TRANSPOSE_TEST tests R83V_TRANSPOSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  real ( kind = 8 ), allocatable :: b1(:)
  real ( kind = 8 ), allocatable :: b2(:)
  real ( kind = 8 ), allocatable :: b3(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TRANSPOSE_TEST'
  write ( *, '(a)' ) &
    '  R83V_TRANSPOSE makes a transposed copy of an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  A:' )

    allocate ( b1(1:min(n-1,m)) )
    allocate ( b2(1:min(n,m)) )
    allocate ( b3(1:min(n,m-1)) )

    call r83v_transpose ( m, n, a1, a2, a3, b1, b2, b3 )
    call r83v_print ( n, m, b1, b2, b3, '  B = tranposed copy of A:' )

    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )
    deallocate ( b1 )
    deallocate ( b2 )
    deallocate ( b3 )

  end do

  return
end
subroutine r83v_zeros ( m, n, a, b, c )

!*****************************************************************************80
!
!! R83V_ZEROS zeros an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the linear system.
!
!    Output, real ( kind = 8 ) A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n))
  real ( kind = 8 ) b(min(m,n))
  real ( kind = 8 ) c(min(m,n-1))

  a(1:min(m-1,n)) = 0.0D+00
  b(1:min(m,n)) = 0.0D+00
  c(1:min(m,n-1)) = 0.0D+00
 
  return
end
subroutine r83v_zeros_test ( )

!*****************************************************************************80
!
!! R83V_ZEROS_TEST tests R83V_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_ZEROS_TEST'
  write ( *, '(a)' ) '  R83V_ZEROS zeros an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_zeros ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Zeroed R83V matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

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
subroutine r8ge_indicator_test ( )

!*****************************************************************************80
!
!! R8GE_INDICATOR_TEST tests R8GE_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8GE_INDICATOR sets up the indicator matrix.'

  call r8ge_indicator ( m, n, a )

  call r8ge_print ( m, n, a, '  The R8GE indicator matrix:' )

  return
end
subroutine r8ge_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MTV multiplies an R8VEC by an R8GE matrix.
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
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8ge_mtv_test ( )

!*****************************************************************************80
!
!! R8GE_MTV_TEST tests R8GE_MTV
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8GE_MTV_TEST'
  write ( *, '(a)' ) '  R8GE_MTV computes a product b=A''*x for an R8GE matrix.'

  call r8ge_indicator ( m, n, a )
  call r8ge_print ( m, n, a, '  The R8GE matrix A:' )

  call r8vec_indicator1 ( m, x )
  call r8vec_print ( m, x, '  Vector x:' )

  call r8ge_mtv ( m, n, a, x, b )
  call r8vec_print ( n, b, '  Vector b = A''*x:' )

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
subroutine r8ge_mv_test ( )

!*****************************************************************************80
!
!! R8GE_MV_TEST tests R8GE_MV
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8GE_MV_TEST'
  write ( *, '(a)' ) '  R8GE_MV computes a product b=A*x for an R8GE matrix.'

  call r8ge_indicator ( m, n, a )
  call r8ge_print ( m, n, a, '  The R8GE matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )

  call r8ge_mv ( m, n, a, x, b )
  call r8vec_print ( m, b, '  Vector b = A*x:' )

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
subroutine r8ge_print_test ( )

!*****************************************************************************80
!
!! R8GE_PRINT_TEST tests R8GE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_PRINT_TEST'
  write ( *, '(a)' ) '  R8GET_PRINT prints an R8GE matrix.'

  do j = 1, n
    do i = 1, m
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  call r8ge_print ( m, n, a, '  The matrix:' )

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
subroutine r8ge_print_some_test ( )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME_TEST tests R8GE_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8GE_PRINT_SOME prints some of an R8GE matrix.'

  do j = 1, n
    do i = 1, m
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  call r8ge_print_some ( m, n, a, 2, 1, 4, 2, &
    '  The matrix, rows 2:4, cols 1:2:' )

  return
end
subroutine r8ge_to_r83v ( m, n, a, a1, a2, a3 )

!*****************************************************************************80
!
!! R8GE_TO_R83V copies (some of) an R8GE matrix to an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), 
!    the R83V matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a1(min(m-1,n))
  real ( kind = 8 ) a2(min(m,n))
  real ( kind = 8 ) a3(min(m,n-1))
  integer ( kind = 4 ) ahi
  integer ( kind = 4 ) bhi
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  ahi = min ( m - 1, n )
  bhi = min ( m,     n )
  chi = min ( m,     n - 1 )

  do k = 1, ahi
    a1(k) = a(k+1,k)
  end do

  do k = 1, bhi
    a2(k) = a(k,k)
  end do

  do k = 1, chi
    a3(k) = a(k,k+1)
  end do
   
  return
end
subroutine r8ge_to_r83v_test ( )

!*****************************************************************************80
!
!! R8GE_TO_R83V_TEST tests R8GE_TO_R83V.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8GE_TO_R83V_TEST'
  write ( *, '(a)' ) '  R8GE_TO_R83V copies an R8GE matrix to an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:m,1:n) )

    call r8ge_indicator ( m, n, a )
    call r8ge_print ( m, n, a, '  R8GE matrix A:' )

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r8ge_to_r83v ( m, n, a, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  R83V copy of (some of ) matrix A:' )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

  end do

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
subroutine r8vec_indicator1_test ( )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1_TEST tests R8VEC_INDICATOR1.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_INDICATOR1_TEST'
  write ( *, '(a)' ) '  R8VEC_INDICATOR1 returns an indicator1 vector.'

  call r8vec_indicator1 ( n, a )

  call r8vec_print ( n, a, '  The indicator1 vector:' )

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_norm_test ( )

!*****************************************************************************80
!
!! R8VEC_NORM_TEST tests R8VEC_NORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) r8vec_norm
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_NORM_TEST'
  write ( *, '(a)' ) '  R8VEC_NORM computes the L2 norm of an R8VEC.'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, x )

  call r8vec_print ( n, x, '  Vector X:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8VEC_NORM ( X ) = ', r8vec_norm ( n, x )

  return
end
function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the vectors.
!
!    Input, real ( kind = 8 ) V0(N), the base vector.
!
!    Input, real ( kind = 8 ) V1(N), the vector whose affine norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end
subroutine r8vec_norm_affine_test ( )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE_TEST tests R8VEC_NORM_AFFINE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_NORM_AFFINE_TEST'
  write ( *, '(a)' ) '  R8VEC_NORM_AFFINE computes the L2 norm of '
  write ( *, '(a)' ) '  the difference of two R8VECs.'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, x )
  call r8vec_uniform_01 ( n, seed, y )
  z(1:n) = x(1:n) - y(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8VEC_NORM_AFFINE(X,Y) = ', &
    r8vec_norm_affine ( n, x, y )
  write ( *, '(a,g14.6)' ) '  R8VEC_NORM (X-Y):        ', r8vec_norm ( n, z )

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
subroutine r8vec_print_test ( )

!*****************************************************************************80
!
!! R8VEC_PRINT_TEST tests R8VEC_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), dimension ( n ) :: a = (/ &
    123.456D+00, 0.000005D+00, -1.0D+06, 3.14159265D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_PRINT_TEST'
  write ( *, '(a)' ) '  R8VEC_PRINT prints an R8VEC.'

  call r8vec_print ( n, a, '  The R8VEC:' )

  return
end
subroutine r8vec_to_r83v ( m, n, a, a1, a2, a3 )

!*****************************************************************************80
!
!! R8VEC_TO_R83V copies an R8VEC to an R83V matrix.
!
!  Discussion:
!
!    The R83V storage format is used for a tridiagonal matrix.
!    The subdiagonal is in A(min(M-1,N)).
!    The diagonal is in B(min(M,N)).
!    The superdiagonal is in C(min(M,N-1)).
!
!  Example:
!
!    An R83V matrix of order 3x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!
!    An R83 matrix of order 5x5 would be stored:
!
!      B1  C1  **  **  **
!      A1  B2  C2  **  **
!      **  A2  B3  C3  **
!      **  **  A3  B4  C4
!      **  **  **  A4  B5
!
!    An R83 matrix of order 5x3 would be stored:
!
!      B1  C1  **
!      A1  B2  C2
!      **  A2  B3
!      **  **  A3
!      **  **  **
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(min(N-1,M)+min(N,M)+min(N,M-1)), the vector.
!
!    Output, real ( kind = 8 ) A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), 
!    the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(min(m-1,n)+min(m,n)+min(m,n-1))
  real ( kind = 8 ) a1(min(m-1,n))
  real ( kind = 8 ) a2(min(m,n))
  real ( kind = 8 ) a3(min(m,n-1))
  integer ( kind = 4 ) ahi
  integer ( kind = 4 ) bhi
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  ahi = min ( m - 1, n )
  bhi = min ( m,     n )
  chi = min ( m,     n - 1 )

  k = 0
  do j = 1, n
    if ( j <= m + 1 .and. 2 <= j ) then
      k = k + 1
      a3(j-1) = a(k)
    end if
    if ( j <= m ) then
      k = k + 1
      a2(j) = a(k)
    end if
    if ( j <= m - 1 ) then
      k = k + 1
      a1(j) = a(k)
    end if
  end do
   
  return
end
subroutine r8vec_to_r83v_test ( )

!*****************************************************************************80
!
!! R8VEC_TO_R83V_TEST tests R8VEC_TO_R83V.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: a1(:)
  real ( kind = 8 ), allocatable :: a2(:)
  real ( kind = 8 ), allocatable :: a3(:)
  integer ( kind = 4 ) ahi
  integer ( kind = 4 ) bhi
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_TO_R83V_TEST'
  write ( *, '(a)' ) '  R8VEC_TO_R83V copies an R8VEC to an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    ahi = min ( m - 1, n )
    bhi = min ( m,     n )
    chi = min ( m,     n - 1 )

    na = ahi + bhi + chi
    allocate ( a(1:na) )

    call r8vec_indicator1 ( na, a )
    call r8vec_print ( na, a, '  R8VEC:' )

    allocate ( a1(1:ahi) )
    allocate ( a2(1:bhi) )
    allocate ( a3(1:chi) )

    call r8vec_to_r83v ( m, n, a, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  R83V matrix:' )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

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
subroutine r8vec_uniform_01_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01_TEST tests R8VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01 returns a random R8VEC '
  write ( *, '(a)' ) '  with entries in [0,1].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Input SEED = ', seed

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print ( n, r, '  Random R8VEC:' )

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
subroutine r8vec2_print_test ( )

!*****************************************************************************80
!
!! R8VEC2_PRINT_TEST tests R8VEC2_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) :: a(n) = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC2_PRINT_TEST'
  write ( *, '(a)' ) '  R8VEC2_PRINT prints a pair of R8VEC''s.'

  do i = 1, n
    b(i) = a(i) ** 2
    c(i) = sqrt ( a(i) )
  end do

  call r8vec2_print ( n, b, c, '  Squares and square roots:' )

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
