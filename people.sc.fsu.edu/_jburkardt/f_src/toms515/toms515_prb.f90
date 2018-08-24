program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS515_PRB.
!
!  Discussion:
!
!    TOMS515_PRB tests the TOMS515 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS515_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS515 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS515_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests COMB by generating all 3-subsets of a 5 set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  parameter ( n = 5 )
  integer ( kind = 4 ) k
  parameter ( k = 3 )

  integer ( kind = 4 ) binom
  integer ( kind = 4 ) c(k)
  logical i4_choose_check
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmax

  lmax = binom ( n, k )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Generate all K-subsets of an N set.'
  write ( *, '(a,i6)' ) '  K = ', k
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i10)' ) '  LMAX = ', lmax

  if ( .not. i4_choose_check ( n, k ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a)' ) '  The binomial coefficient cannot be'
    write ( *, '(a)' ) '  computed in integer arithmetic for'
    write ( *, '(a)' ) '  this choice of parameters.'
    return
  end if

  write ( *, '(a)' ) ' '

  do l = 1, lmax
    call comb ( n, k, l, c )
    write ( *, * ) l, c(1:k)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests COMB by generating 10 random 3-subsets of a 10 set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  parameter ( n = 5 )
  integer ( kind = 4 ) k
  parameter ( k = 3 )

  integer ( kind = 4 ) binom
  integer ( kind = 4 ) c(k)
  integer ( kind = 4 ) i
  logical i4_choose_check
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) seed

  lmax = binom ( n, k )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Generate 10 random K-subsets of an N set.'
  write ( *, '(a,i6)' ) '  K = ', k
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i10)' ) '  LMAX = ', lmax

  if ( .not. i4_choose_check ( n, k ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST02 - Warning!'
    write ( *, '(a)' ) '  The binomial coefficient cannot be'
    write ( *, '(a)' ) '  computed in integer arithmetic for'
    write ( *, '(a)' ) '  this choice of parameters.'
    return
  end if

  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    l = i4_uniform_ab ( 1, lmax, seed )
    call comb ( n, k, l, c )
    write ( *, * ) l, c(1:k)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests COMB by generating 10 random 3-subsets of a 25 set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  parameter ( n = 25 )
  integer ( kind = 4 ) k
  parameter ( k = 3 )

  integer ( kind = 4 ) binom
  integer ( kind = 4 ) c(k)
  integer ( kind = 4 ) i
  logical i4_choose_check
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) seed

  lmax = binom ( n, k )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Generate 10 random K-subsets of an N set.'
  write ( *, '(a,i6)' ) '  K = ', k
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i10)' ) '  LMAX = ', lmax

  if ( .not. i4_choose_check ( n, k ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST03 - Warning!'
    write ( *, '(a)' ) '  The binomial coefficient cannot be'
    write ( *, '(a)' ) '  computed in integer arithmetic for'
    write ( *, '(a)' ) '  this choice of parameters.'
    return
  end if

  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    l = i4_uniform_ab ( 1, lmax, seed )
    call comb ( n, k, l, c )
    write ( *, * ) l, c(1:k)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests COMB by generating 10 random 3-subsets of a 100 set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  parameter ( n = 100 )
  integer ( kind = 4 ) k
  parameter ( k = 3 )

  integer ( kind = 4 ) binom
  integer ( kind = 4 ) c(k)
  integer ( kind = 4 ) i
  logical i4_choose_check
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) seed

  lmax = binom ( n, k )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Generate 10 random K-subsets of an N set.'
  write ( *, '(a,i6)' ) '  K = ', k
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i10)' ) '  LMAX = ', lmax

  if ( .not. i4_choose_check ( n, k ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST04 - Warning!'
    write ( *, '(a)' ) '  The binomial coefficient cannot be'
    write ( *, '(a)' ) '  computed in integer arithmetic for'
    write ( *, '(a)' ) '  this choice of parameters.'
    return
  end if

  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    l = i4_uniform_ab ( 1, lmax, seed )
    call comb ( n, k, l, c )
    write ( *, * ) l, c(1:k)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests COMB by generating 10 random 10-subsets of a 100 set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  parameter ( n = 100 )
  integer ( kind = 4 ) k
  parameter ( k = 10 )

  integer ( kind = 4 ) binom
  integer ( kind = 4 ) c(k)
  integer ( kind = 4 ) i
  logical i4_choose_check
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) seed

  lmax = binom ( n, k )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Generate 10 random K-subsets of an N set.'
  write ( *, '(a,i6)' ) '  K = ', k
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i10)' ) '  LMAX = ', lmax
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that this function is already'
  write ( *, '(a)' ) '  failing because LMAX is negative.'
  write ( *, '(a)' ) '  The combinatorial coefficient C(100,10)'
  write ( *, '(a)' ) '  is too large to store in an integer.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Although the program continues to give'
  write ( *, '(a)' ) '  results, they cannot be relied on!'

  if ( .not. i4_choose_check ( n, k ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a)' ) '  The binomial coefficient cannot be'
    write ( *, '(a)' ) '  computed in integer arithmetic for'
    write ( *, '(a)' ) '  this choice of parameters.'
    return
  end if

  write ( *, '(a)' ) ''

  seed = 123456789

  do i = 1, 10
    l = i4_uniform_ab ( 1, lmax, seed )
    call comb ( n, k, l, c )
    write ( *, * ) l, c(1:k)
  end do

  return
end

