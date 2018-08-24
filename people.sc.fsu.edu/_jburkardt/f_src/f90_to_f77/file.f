      function r8_invgam_pdf ( beta, alpha, rval )

c*****************************************************************************80
c
c! R8_INVGAM_PDF evaluates the PDF of an inverse gamma distribution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Original FORTRAN90 version by Guannan Zhang.
c    Modifications by John Burkardt.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) BETA, the rate parameter.
c    0.0 < BETA.
c
c    Input, real ( kind = 8 ) ALPHA, the shape parameter.
c    0.0 < ALPHA.
c
c    Input, real ( kind = 8 ) RVAL, the point where the PDF is evaluated.
c
c    Output, real ( kind = 8 ) R8_INVGAM_PDF, the value of the PDF at RVAL.
c
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

        temp = alpha * log ( beta ) - ( alpha + 1.0D+00 ) * log ( rval )
     &       - beta / rval - r8_gamma_log ( alpha )

        r8_invgam_pdf = exp ( temp )

      end if

      return
      end
