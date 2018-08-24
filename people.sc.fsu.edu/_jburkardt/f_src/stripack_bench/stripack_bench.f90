program main

!*****************************************************************************80
!
!! MAIN is the main program for STRIPACK_BENCH.
!
!  Discussion:
!
!    STRIPACK_BENCH times the Delaunay mesh calculation carried out by
!    STRIPACK, for random data sets of increasing size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) bisect
  real ( kind = 8 ), allocatable :: dist(:)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), allocatable :: lend(:)
  integer ( kind = 4 ), allocatable :: list(:)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ), allocatable :: lptr(:)
  integer ( kind = 4 ), allocatable :: ltri(:,:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable :: near(:)
  integer ( kind = 4 ), allocatable :: next(:)
  integer ( kind = 4 ), parameter :: nrow = 6
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) seed
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) time3
  real ( kind = 8 ), allocatable :: xyz(:,:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_BENCH'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Time STRIPACK''s Delaunay triangularion calculation'
  write ( *, '(a)' ) '  for sets of random nodes of increasing size.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        N        Time'
  write ( *, '(a)' ) ''
!
!  Allocate space for the nodes.
!
  n = 12
!
!  Originally, we wanted to go to level 10, but already level 8 takes
!  30 minutes.
!
  do bisect = 0, 8

    allocate ( xyz(1:3,1:n) )
!
!  Generate random node values on the sphere.
!
    seed = 123456789
    call uniform_on_sphere01_map ( 3, n, seed, xyz )

    call cpu_time ( time1 )
!
!  Compute the Delaunay triangulation.
!
    allocate ( dist(1:n) )
    allocate ( lend(1:n) )
    allocate ( list(1:6*(n-2)) )
    allocate ( lptr(1:6*(n-2)) )
    allocate ( ltri(1:nrow,1:2*(n-2)) )
    allocate ( near(1:n) )
    allocate ( next(1:n) )

    call trmesh ( n, xyz(1,:), xyz(2,:), xyz(3,:), list, lptr, lend, lnew, near, &
      next, dist, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPACK_RANDOM - Fatal error!'
      write ( *, '(a,i8)' ) '  TRMESH returns IER = ', ier
      stop 1
    end if
!
!  Create the triangle list.
!
    call trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPACK_BENCH - Fatal error!'
      write ( *, '(a,i8)' ) '  TRLIST returns IER = ', ier
      stop 1
    end if

    call cpu_time ( time2 )

    time3 = time2 - time1

    write ( *, '(2x,i8,2x,g14.6)') n, time3

    deallocate ( dist )
    deallocate ( lend )
    deallocate ( list )
    deallocate ( lptr )
    deallocate ( ltri )
    deallocate ( near )
    deallocate ( next )
    deallocate ( xyz )

    n = n + ( 4 ** bisect ) * 30

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_BENCH:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
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
!    31 May 2007
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
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

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
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range 
!    of entries of X that we need to compute.  This starts off as 1:N, but 
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: two = 2
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, two ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
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
!    Volume 8, 1969, pages 136-143.
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
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

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
subroutine uniform_on_sphere01_map ( dim_num, n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This procedure is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,n)

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( dim_num, seed, x(1:dim_num,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm

  end do

  return
end
