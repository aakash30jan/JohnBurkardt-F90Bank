function r8_invgam_pdf ( beta, alpha, rval )

!*****************************************************************************80
!
!! R8_INVGAM_PDF evaluates the PDF of an inverse gamma distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2013
!
!  Author:
!
!    Original FORTRAN90 version by Guannan Zhang.
!    Modifications by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA, the rate parameter.
!    0.0 < BETA.
!
!    Input, real ( kind = 8 ) ALPHA, the shape parameter.
!    0.0 < ALPHA.
!
!    Input, real ( kind = 8 ) RVAL, the point where the PDF is evaluated.
!
!    Output, real ( kind = 8 ) R8_INVGAM_PDF, the value of the PDF at RVAL.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) r8_invgam_pdf
  real ( kind = 8 ) rval
  real ( kind = 8 ) temp

  if ( alpha <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_INVGAM_PDF - Fatal error!'
    write ( *, '(a)' ) '  Parameter ALPHA is not positive.'
    stop
  end if

  if ( beta <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_INVGAM_PDF - Fatal error!'
    write ( *, '(a)' ) '  Parameter BETA is not positive.'
    stop
  end if

  if ( rval <= 0.0D+00 ) then

    r8_invgam_pdf = 0.0D+00

  else

    temp = alpha * log ( beta ) - ( alpha + 1.0D+00 ) * log ( rval ) &
      - beta / rval - r8_gamma_log ( alpha )

    r8_invgam_pdf = exp ( temp )

  end if

  return
end
