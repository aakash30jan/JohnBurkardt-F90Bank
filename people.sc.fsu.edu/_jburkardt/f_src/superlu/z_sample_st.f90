program main

!*****************************************************************************80
!
!! Z_SAMPLE_ST calls SUPERLU to solve a system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: nst = 12

  complex ( kind = 8 ), allocatable :: acc(:)
  complex ( kind = 8 ), dimension ( nst ) :: ast = (/ &
    ( 19.0,  0.0 ), ( 12.0,  0.0 ), ( 12.0,  0.0 ), &
    ( 21.0,  0.0 ), ( 12.0,  0.0 ), ( 12.0,  0.0 ), &
    ( 21.0,  0.0 ), ( 16.0,  0.0 ), &
    ( 21.0,  0.0 ), (  5.0,  0.0 ), &
    ( 21.0,  0.0 ), ( 18.0,  0.0 ) /)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) b2(n)
  integer ( kind = 4 ), allocatable :: ccc(:)
  integer ( kind = 4 ) factors(8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: icc(:)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ), dimension ( nst ) :: ist = (/ &
    1, 2, 5, &
    2, 3, 5, &
    1, 3, &
    1, 4, &
    4, 5 /)
  integer ( kind = 4 ), dimension ( nst ) :: jst = (/ &
    1, 1, 1, &
    2, 2, 2, &
    3, 3, &
    4, 4, &
    5, 5 /)
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nrhs

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Z_SAMPLE_ST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  ZGSSV factors and solves a linear system'
  write ( *, '(a)' ) '  using single precision complex arithmetic.'

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Matrix order N = ', n
  write ( *, '(a,i6)' ) '  Matrix nonzeros NST = ', nst
!
!  Print the ST matrix.
!
  call st_print ( m, n, nst, ist, jst, ast, '  ST matrix:' )
!
!  Get the CC size.
!
  call st_to_cc_size ( nst, ist, jst, ncc )

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of CC values NCC = ', ncc
!
!  Create the CC indices.
!
  allocate ( icc(1:ncc) )
  allocate ( ccc(1:n+1) )

  call st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc )
!
!  Create the CC values.
!
  allocate ( acc(1:ncc) )

  call st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc, acc )

  nrhs = 1
  ldb = n
  do i = 1, n
    b(i) = ( 1.0E+00, 0.0E+00 )
  end do
!
!  Factor the matrix.
!
  iopt = 1
  call c_fortran_zgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Z_SAMPLE_ST - Fatal error!'
    write ( *, '(a)' ) '  Factorization failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
  end if

  write ( *, '(a)' ) '  Factorization succeeded.'
!
!  Solve the factored system.
!
  iopt = 2
  call c_fortran_zgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Z_SAMPLE_ST - Fatal error!'
    write ( *, '(a)' ) '  Backsolve failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) b(i)
  end do
!
!  B now contains the solution X.
!  Set B2 = A * X.
!
  call cc_mv ( m, n, ncc, icc, ccc, acc, b, b2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Product A*X:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) b2(i)
  end do
!
!  Free memory.
!
  iopt = 3
  call c_fortran_zgssv ( iopt, n, ncc, nrhs, acc, icc, &
    ccc, b, ldb, factors, info )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Z_SAMPLE_ST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine cc_mv ( m, n, ncc, icc, ccc, acc, x, b )

!*****************************************************************************80
!
!! CC_MV multiplies a CC matrix by a vector
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC values.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the compressed CC columns
!
!    Input, complex ( kind = 8 ) ACC(NCC), the CC values.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, complex ( kind = 8 ) B(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc

  complex ( kind = 8 ) acc(ncc)
  complex ( kind = 8 ) b(m)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  complex ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do k = ccc(j), ccc(j+1) - 1
      i = icc(k)
      b(i) = b(i) + acc(k) * x(j)
    end do
  end do

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares entries of an I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
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
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

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

      t     = a1(i)
      a1(i) = a1(j)
      a1(j) = t

      t     = a2(i)
      a2(i) = a2(j)
      a2(j) = t
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique_count ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE_COUNT counts unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), the items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  iu = 1
  unique_num = 1

  do i = 2, n

    if ( a1(i) /= a1(iu) .or. a2(i) /= a2(iu) ) then
      unique_num = unique_num + 1
      iu = i
    end if

  end do

  return
end
subroutine i4vec2_sorted_uniquely ( n1, a1, b1, n2, a2, b2 )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUELY copies unique elements from a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, the number of items.
!
!    Input, integer ( kind = 4 ) A1(N1), B1(N1), the array of items.
!
!    Input, integer ( kind = 4 ) N2, the number of unique items.
!
!    Output, integer ( kind = 4 ) A2(N2), B2(N2), the array of unique items.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  integer ( kind = 4 ) a1(n1)
  integer ( kind = 4 ) a2(n2)
  integer ( kind = 4 ) b1(n1)
  integer ( kind = 4 ) b2(n2)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2

  i1 = 1
  i2 = 1
  a2(i2) = a1(i1)
  b2(i2) = b1(i1)

  do i1 = 2, n1

    if ( a1(i1) /= a2(i2) .or. b1(i1) /= b2(i2) ) then

      i2 = i2 + 1

      a2(i2) = a1(i1)
      b2(i2) = b1(i1)

    end if

  end do

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
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
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)

  a2(1:n) = a1(1:n)

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
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!    On return, if INDX is
!    *greater than 0,
!    ...interchange items I and J;
!    ...call again.
!    *less than 0,
!    ...compare items I and J;
!    ...set ISGN = -1 if I < J, ISGN = +1 if J < I;
!    ...call again.
!    * equal to 0, 
!    ...the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
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
subroutine st_print ( m, n, nst, ist, jst, ast, title )

!*****************************************************************************80
!
!! ST_PRINT prints a sparse matrix in ST format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Input, complex ( kind = 8 ) AST(NST), the ST values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nst

  complex ( kind = 8 ) ast(nst)
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '     #     I     J       Ar              Ai'
  write ( *, '(a)' ) '  ----  ----  ----  --------------  --------------'
  write ( *, '(a)' ) ' '

  do k = 1, nst
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8,2x,g16.8)' ) &
      k, ist(k), jst(k), ast(k)
  end do

  return
end
subroutine st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc )

!*****************************************************************************80
!
!! ST_TO_CC_INDEX creates CC indices from ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Output, integer ( kind = 4 ) CCC(N+1), the compressed CC columns.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst

  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) ist2(nst)
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jcc(ncc)
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) jst2(nst)
!
!  Make copies so the sorting doesn't confuse the user.
!
  call i4vec_copy ( nst, ist, ist2 )
  call i4vec_copy ( nst, jst, jst2 )
!
!  Sort the elements.
!
  call i4vec2_sort_a ( nst, jst2, ist2 )
!
!  Get the unique elements.
!
  call i4vec2_sorted_uniquely ( nst, jst2, ist2, ncc, jcc, icc )
!
!  Compress the column index.
!
  ccc(1) = 1
  jlo = 1
  do i = 1, ncc
    jhi = jcc(i)
    if ( jhi /= jlo ) then
      ccc(jlo+1:jhi) = i
      jlo = jhi
    end if
  end do
  jhi = n + 1
  ccc(jlo+1:jhi) = ncc + 1

  return
end
subroutine st_to_cc_size ( nst, ist, jst, ncc )

!*****************************************************************************80
!
!! ST_TO_CC_SIZE sizes CC indexes based on ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Output, integer ( kind = 4 ) NCC, the number of CC elements.
!
  implicit none

  integer ( kind = 4 ) nst

  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) ist2(nst)
  integer ( kind = 4 ) jst2(nst)
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) ncc
!
!  Make copies so the sorting doesn't confuse the user.
!
  call i4vec_copy ( nst, ist, ist2 )
  call i4vec_copy ( nst, jst, jst2 )
!
!  Sort by column first, then row.
!
  call i4vec2_sort_a ( nst, jst2, ist2 )
!
!  Count the unique pairs.
!
  call i4vec2_sorted_unique_count ( nst, jst2, ist2, ncc )

  return
end
subroutine st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc, acc )

!*****************************************************************************80
!
!! ST_TO_CC_VALUES creates CC values from ST data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the number of ST elements.
!
!    Input, integer ( kind = 4 ) IST(NST), JST(NST), the ST rows and columns.
!
!    Input, complex ( kind = 8 ) AST(NST), the ST values.
!
!    Input, integer ( kind = 4 ) NCC, the number of CC elements.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, integer ( kind = 4 ) ICC(NCC), the CC rows.
!
!    Input, integer ( kind = 4 ) CCC(N+1), the CC compressed columns.
!
!    Output, complex ( kind = 8 ) ACC(NCC), the CC values.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst

  complex ( kind = 8 ) ast(nst)
  complex ( kind = 8 ) acc(ncc)
  integer ( kind = 4 ) ccc(n+1)
  integer ( kind = 4 ) chi
  integer ( kind = 4 ) clo
  logical fail
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icc(ncc)
  integer ( kind = 4 ) ist(nst)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jst(nst)
  integer ( kind = 4 ) kcc
  integer ( kind = 4 ) kst

  acc(1:ncc) = 0.0D+00

  do kst = 1, nst

    i = ist(kst)
    j = jst(kst)

    clo = ccc(j)
    chi = ccc(j+1)

    fail = .true.

    do kcc = clo, chi - 1
      if ( icc(kcc) == i ) then
        acc(kcc) = acc(kcc) + ast(kst)
        fail = .false.
        exit
      end if                  
    end do

    if ( fail ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'ST_TO_CC_VALUES - Fatal error!'
      write ( *, '(a)' ) '  ST entry cannot be located in CC array.'
      write ( *, '(a,i4)' ) '  ST index KST    = ', kst
      write ( *, '(a,i4)' ) '  ST row IST(KST) = ', ist(kst)
      write ( *, '(a,i4)' ) '  ST col JST(KST) = ', jst(kst)
      write ( *, '(a,2g14.6)' ) '  ST val AST(KST) = ', ast(kst)
      stop 1
    end if

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
