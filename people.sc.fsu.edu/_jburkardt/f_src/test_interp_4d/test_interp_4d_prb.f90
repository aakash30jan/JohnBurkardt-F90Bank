program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INTERP_4D_PRB.
!
!  Discussion:
!
!    TEST_INTERP_4D_PRB tests the TEST_INTERP_4D library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_4D_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INTERP_4D library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_4D_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 shows how P00_TITLE can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate some of the bookkeeping routines.'
  write ( *, '(a)' ) '  P00_PROB_NUM returns the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns the problem title.'
  write ( *, '(a)' ) '  P00_LIMIT returns the problem limits.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of problems = ', prob_num

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    call p00_title ( prob, title )
    write ( *, '(a)' ) '  Problem TITLE = "' // trim ( title ) // '".'
    call p00_lim_4d ( prob, a, b )
    write ( *, '(a,g14.6)' ) '  Problem lower limit A = ', a
    write ( *, '(a,g14.6)' ) '  Problem upper limit B = ', b

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 shows how P00_STORY can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_STORY prints the problem "story".'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob

    call p00_story ( prob )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses nearest neighbor interpolation on a regular grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) e
  real ( kind = 8 ) fd
  real ( kind = 8 ) fs
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ng
  integer ( kind = 4 ) ng_log
  integer ( kind = 4 ) ngm
  integer ( kind = 4 ) n
  real ( kind = 8 ) p00_fun_4d
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_round
  integer ( kind = 4 ) seed
  character ( len = 80 ) title
  real ( kind = 8 ) xd(m)
  real ( kind = 8 ) xs(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Nearest neighbor interpolation.'
  write ( *, '(a)' ) '  Evaluate the function on an NxN equally spaced grid.'
  write ( *, '(a)' ) '  The interpolated value at any point is the value'
  write ( *, '(a)' ) '  at the nearest interpolating node.'
  write ( *, '(a)' ) '  Estimate the integral of the square of the error'
  write ( *, '(a)' ) '  between the function and the interpolant.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )
    call p00_lim_4d ( prob, a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )
    do i = 1, m
      write ( *, '(g14.6,a,i1,a,g14.6)' ) a(i), ' <= X(', i, ') <= ', b(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     Grid Side  RMS error'
    write ( *, '(a)' ) ' '

    ng = 1

    do ng_log = 0, 4

      n = 50
      seed = 123456789
      e = 0.0D+00
    
      do j = 1, n

        call r8vec_uniform_01 ( m, seed, xs )
!
!  Discretize the sample point XS to a point XD on the grid.
!
        do i = 1, m
          h = ( b(i) - a(i) ) / real ( ng, kind = 8 )
          r = ( xs(i) - a(i) ) / h
          k = r8_round ( r )
          k = min ( k, ng )
          k = max ( k, 0 )
          xd(i) = a(i) + real ( k, kind = 8 ) * h
        end do

        fs = p00_fun_4d ( prob, xs )
        fd = p00_fun_4d ( prob, xd )
 
        e = e + ( fs - fd ) ** 2
      
      end do

      e = sqrt ( e ) * product ( b(1:m) - a(1:m) ) / real ( n, kind = 8 )

      ngm = ng ** m

      write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) ng, ngm, e

      ng = ng * 2

    end do

  end do

  return
end
