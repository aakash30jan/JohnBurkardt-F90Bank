program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS179_PRB.
!
!  Discussion:
!
!    TOMS179_PRB tests the TOMS179 library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS179_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS179 library.'

  call alogam_test ( )
  call mdbeta_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS179_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine alogam_test ( )

!*****************************************************************************80
!
!! ALOGAM_TEST tests ALOGAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alogam
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALOGAM_TEST'
  write ( *, '(a)' ) '  ALOGAM estimates the logarithm of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X         Exact Value               ', &
    'Computed                Diff'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    fx2 = alogam ( x, ifault )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine mdbeta_test ( )

!*****************************************************************************80
!
!! MDBETA_TEST tests MDBETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MDBETA_TEST'
  write ( *, '(a)' ) '  MDBETA_TEST estimates the modified Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X         P         Q         ', &
    'Exact Value               Computed                Diff'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_cdf_values ( n_data, p, q, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    call mdbeta ( x, p, q, fx2, ier )

    write ( *, &
      '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, p, q, fx, fx2, abs ( fx - fx2 )

  end do

  return
end

