function c4_abs ( z )

!*****************************************************************************80
!
!! C4_ABS evaluates the absolute value of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the absolute value of a C4 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) C4_ABS, the function value.
!
  implicit none

  real ( kind = 4 ) c4_abs
  complex ( kind = 4 ) z

  c4_abs = sqrt ( ( real ( z, kind = 4 ) )**2 &
                + ( aimag ( z ) )**2 )

  return
end
function c4_acos ( z )

!*****************************************************************************80
!
!! C4_ACOS evaluates the inverse cosine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse cosine of a C4.
!
!    Here we use the relationship:
!
!       C4_ACOS ( Z ) = pi/2 - C4_ASIN ( Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ACOS, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_acos
  complex ( kind = 4 ) c4_asin
  real ( kind = 4 ), parameter :: r4_pi_half = 1.57079632679489661923E+00
  complex ( kind = 4 ) z

  c4_acos = r4_pi_half - c4_asin ( z )

  return
end
function c4_acosh ( z )

!*****************************************************************************80
!
!! C4_ACOSH evaluates the inverse hyperbolic cosine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic cosine of a C4.
!
!    Here we use the relationship:
!
!      C4_ACOSH ( Z ) = i * C4_ACOS ( Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ACOSH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_acos
  complex ( kind = 4 ) c4_acosh
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  c4_acosh = c4_i * c4_acos ( z )

  return
end
function c4_add ( z1, z2 )

!*****************************************************************************80
!
!! C4_ADD adds two C4's.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports addition of C4's with the "+" operator.
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
!    Input, complex ( kind = 4 ) Z1, Z2, the values to add.
!
!    Output, complex ( kind = 4 ) C4_ADD, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_add
  complex ( kind = 4 ) z1
  complex ( kind = 4 ) z2

  c4_add = z1 + z2

  return
end
function c4_arg ( x )

!*****************************************************************************80
!
!! C4_ARG returns the argument of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the argument of a C4.
!
!    By convention, the argument of a C4 is expected to lie between
!    -PI and PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the complex number.
!
!    Output, real ( kind = 4 ) C4_ARG, the function value.
!
  implicit none

  real ( kind = 4 ) c4_arg
  complex ( kind = 4 ) x

  if ( aimag ( x )           == 0.0E+00 .and. &
       real  ( x, kind = 4 ) == 0.0E+00 ) then

    c4_arg = 0.0E+00

  else

    c4_arg = atan2 ( aimag ( x ), real ( x ) )

  end if

  return
end
function c4_asin ( z )

!*****************************************************************************80
!
!! C4_ASIN evaluates the inverse sine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse sine of a C4.
!
!    Here we use the relationship:
!
!      C4_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ASIN, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_asin
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  c4_asin = - c4_i * log ( c4_i * z + sqrt ( 1.0E+00 - z * z ) )

  return
end
function c4_asinh ( z )

!*****************************************************************************80
!
!! C4_ASINH evaluates the inverse hyperbolic sine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic sine of a C4.
!
!    Here we use the relationship:
!
!      C4_ASINH ( Z ) = - i * C4_ASIN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ASINH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_asin
  complex ( kind = 4 ) c4_asinh
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  c4_asinh = - c4_i * c4_asin ( c4_i * z )

  return
end
function c4_atan ( z )

!*****************************************************************************80
!
!! C4_ATAN evaluates the inverse tangent of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse tangent of a C4.
!
!    Here we use the relationship:
!
!      C4_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ATAN, the function value.
!
  implicit none

  complex ( kind = 4 ) arg
  complex ( kind = 4 ) c4_atan
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  arg = ( 1.0E+00 - c4_i * z ) / ( 1.0E+00 + c4_i * z )

  c4_atan = 0.5E+00 * c4_i * log ( arg )

  return
end
function c4_atanh ( z )

!*****************************************************************************80
!
!! C4_ATANH evaluates the inverse hyperbolic tangent of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the inverse hyperbolic tangent of a C4.
!
!    Here we use the relationship:
!
!      C4_ATANH ( Z ) = - i * C4_ATAN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_ATANH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_atan
  complex ( kind = 4 ) c4_atanh
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  c4_atanh = - c4_i * c4_atan ( c4_i * z )

  return
end
function c4_conj ( z )

!*****************************************************************************80
!
!! C4_CONJ evaluates the conjugate of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the conjugate of a C4 with the CONJG function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_CONJ, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_conj
  complex ( kind = 4 ) z

  c4_conj = cmplx ( real ( z, kind = 4 ), - aimag ( z ) )

  return
end
function c4_copy ( z )

!*****************************************************************************80
!
!! C4_COPY copies a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the copy of a C4 with the "=" operator.
!
!    The order of the arguments may seem unnatural, but it is arranged so
!    that the call
!
!      c4_copy ( c1, c2 )
!
!    mimics the assignment
!
!      c1 = c2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_COPY, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_copy
  complex ( kind = 4 ) z

  c4_copy = z

  return
end
function c4_cos ( z )

!*****************************************************************************80
!
!! C4_COS evaluates the cosine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the cosine of a C4 with the COS function.
!
!    We use the relationship:
!
!      C4_COS ( C ) = ( C4_EXP ( i * C ) + C4_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_COS, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_cos
  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00 )

  c4_cos = ( c4_exp ( c4_i * z ) + c4_exp ( - c4_i * z ) ) / 2.0E+00

  return
end
function c4_cosh ( z )

!*****************************************************************************80
!
!! C4_COSH evaluates the hyperbolic cosine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the hyperbolic cosine of a C4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_COSH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_cosh
  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) z

  c4_cosh = ( c4_exp ( z ) + c4_exp ( - z ) ) / 2.0E+00

  return
end
function c4_cube_root ( x )

!*****************************************************************************80
!
!! C4_CUBE_ROOT returns the principal cube root of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the cube root of a C4 through the 
!    "**(1.0/3.0)" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the argument.
!
!    Output, complex ( kind = 4 ) C4_CUBE_ROOT, the function value.
!
  implicit none

  real ( kind = 4 ) arg
  real ( kind = 4 ) c4_arg
  complex ( kind = 4 ) c4_cube_root
  real ( kind = 4 ) c4_mag
  real ( kind = 4 ) mag
  complex ( kind = 4 ) x

  arg = c4_arg ( x )
  mag = c4_mag ( x )

  if ( mag == 0.0E+00 ) then

    c4_cube_root = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  else

    c4_cube_root = mag**( 1.0E+00 / 3.0E+00 ) &
      * cmplx ( cos ( arg / 3.0E+00 ), &
                sin ( arg / 3.0E+00 ), kind = 4 )

  end if

  return
end
function c4_div ( z1, z2 )

!*****************************************************************************80
!
!! C4_DIV divides two C4's.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports division of C4's with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z1, Z2, the arguments.
!
!    Output, complex ( kind = 4 ) C4_DIV, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_div
  complex ( kind = 4 ) z1
  complex ( kind = 4 ) z2

  c4_div = z1 / z2

  return
end
function c4_div_r4 ( z1, r )

!*****************************************************************************80
!
!! C4_DIV_R4 divides a C4 by an R4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!    An R4 is a real ( kind = 4 ) value.
!
!    FORTRAN90 supports division of a C4 by an R4 with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z1, the value to be divided.
!
!    Input, real ( kind = 4 ) R, the divisor.
!
!    Output, complex ( kind = 4 ) C4_DIV_R4, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_div_r4
  real ( kind = 4 ) r
  complex ( kind = 4 ) z1

  c4_div_r4 = z1 / r

  return
end
function c4_exp ( z )

!*****************************************************************************80
!
!! C4_EXP evaluates the exponential of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the exponential of a C4 with the EXP function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double complex Z, the argument.
!
!    Output, double complex C4_EXP, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) z
  real ( kind = 4 ) zi
  real ( kind = 4 ) zr

  zr = real ( z, kind = 4 )
  zi = aimag ( z )

  c4_exp = exp ( zr ) * cmplx ( cos ( zi ), sin ( zi ), kind = 4 )

  return
end
function c4_i ( )

!*****************************************************************************80
!
!! C4_I returns the imaginary unit, i as a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 4 ) C4_I, the value of complex i.
!
  implicit none

  complex ( kind = 4 ) c4_i

  c4_i = cmplx ( 0.0E+00, 1.0E+00, kind = 4 )

  return
end
function c4_imag ( z )

!*****************************************************************************80
!
!! C4_IMAG evaluates the imaginary part of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the imaginary part of a C4 with the AIMAG function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) C4_IMAG, the function value.
!
  implicit none

  real ( kind = 4 ) c4_imag
  complex ( kind = 4 ) z

  c4_imag = aimag ( z )

  return
end
function c4_inv ( z )

!*****************************************************************************80
!
!! C4_INV evaluates the inverse of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the inverse of a C4 with the "1/" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_INV, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_inv
  complex ( kind = 4 ) z
  real ( kind = 4 ) z_imag
  real ( kind = 4 ) z_norm
  real ( kind = 4 ) z_real

  z_real = real ( z, kind = 4 )
  z_imag = aimag ( z )

  z_norm = sqrt ( z_real * z_real + z_imag * z_imag )

  c4_inv = cmplx ( z_real, - z_imag ) / z_norm / z_norm

  return
end
function c4_le_l1 ( x, y )

!*****************************************************************************80
!
!! C4_LE_L1 := X <= Y for C4 values, and the L1 norm.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    The L1 norm can be defined here as:
!
!      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, Y, the values to be compared.
!
!    Output, logical C4_LE_L1, is TRUE if X <= Y.
!
  implicit none

  logical c4_le_l1
  complex ( kind = 4 ) x
  complex ( kind = 4 ) y

  if ( abs ( real ( x, kind = 4 ) ) + abs ( aimag ( x ) ) <= &
       abs ( real ( y, kind = 4 ) ) + abs ( aimag ( y ) ) ) then
    c4_le_l1 = .true.
  else
    c4_le_l1 = .false.
  end if

  return
end
function c4_le_l2 ( x, y )

!*****************************************************************************80
!
!! C4_LE_L2 := X <= Y for C4 values, and the L2 norm.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    The L2 norm can be defined here as:
!
!      C4_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( imag (X) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, Y, the values to be compared.
!
!    Output, logical C4_LE_L2, is TRUE if X <= Y.
!
  implicit none

  logical c4_le_l2
  complex ( kind = 4 ) x
  complex ( kind = 4 ) y

  if ( ( real ( x, kind = 4 ) )**2 + ( aimag ( x ) )**2 <= &
       ( real ( y, kind = 4 ) )**2 + ( aimag ( y ) )**2 ) then
    c4_le_l2 = .true.
  else
    c4_le_l2 = .false.
  end if

  return
end
function c4_le_li ( x, y )

!*****************************************************************************80
!
!! C4_LE_LI := X <= Y for C4 values, and the L Infinity norm.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    The L Infinity norm can be defined here as:
!
!      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, Y, the values to be compared.
!
!    Output, logical C4_LE_LI, is TRUE if X <= Y.
!
  implicit none

  logical c4_le_li
  complex ( kind = 4 ) x
  complex ( kind = 4 ) y

  if ( max ( abs ( real ( x, kind = 4 ) ), abs ( aimag ( x ) ) ) <= &
       max ( abs ( real ( y, kind = 4 ) ), abs ( aimag ( y ) ) ) ) then
    c4_le_li = .true.
  else
    c4_le_li = .false.
  end if

  return
end
function c4_log ( z )

!*****************************************************************************80
!
!! C4_LOG evaluates the logarithm of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the logarithm of a C4 with the LOG function.
!
!    Here we use the relationship:
!
!      C4_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_LOG, the function value.
!
  implicit none

  real ( kind = 4 ) arg
  real ( kind = 4 ) c4_arg
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) c4_log
  real ( kind = 4 ) c4_mag
  real ( kind = 4 ) mag
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00 )

  arg = c4_arg ( z )
  mag = c4_mag ( z )

  c4_log = log ( mag ) + c4_i * arg
 
  return
end
function c4_mag ( x )

!*****************************************************************************80
!
!! C4_MAG returns the magnitude of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the magnitude of a C4 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the argument.
!
!    Output, real ( kind = 4 ) C4_MAG, the function value.
!
  implicit none

  real ( kind = 4 ) c4_mag
  complex ( kind = 4 ) x

  c4_mag = sqrt ( ( real ( x, kind = 4 ) )**2 + ( aimag ( x ) )**2 )

  return
end
function c4_mul ( z1, z2 )

!*****************************************************************************80
!
!! C4_MUL multiplies two C4's.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports multiplication of C4's with the "*" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z1, Z2, the values to multiply.
!
!    Output, complex ( kind = 4 ) C4_MUL, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_mul
  complex ( kind = 4 ) z1
  complex ( kind = 4 ) z2

  c4_mul = z1 * z2

  return
end
function c4_neg ( c1 )

!*****************************************************************************80
!
!! C4_NEG returns the negative of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports negation of a C4 with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) C1, the value to be negated.
!
!    Output, complex ( kind = 4 ) C4_NEG, the function value.
!
  implicit none

  complex ( kind = 4 ) c1
  complex ( kind = 4 ) c4_neg

  c4_neg = - c1

  return
end
function c4_nint ( c1 )

!*****************************************************************************80
!
!! C4_NINT returns the nearest complex integer of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) C1, the value to be NINT'ed.
!
!    Output, complex ( kind = 4 ) C4_NINT, the NINT'ed value.
!
  implicit none

  complex ( kind = 4 ) c1
  complex ( kind = 4 ) c4_nint
  real ( kind = 4 ) r
  real ( kind = 4 ) r_min
  real ( kind = 4 ) r4_floor
  real ( kind = 4 ) x
  real ( kind = 4 ) x_min
  real ( kind = 4 ) xc
  real ( kind = 4 ) y
  real ( kind = 4 ) y_min
  real ( kind = 4 ) yc

  xc = real ( c1 )
  yc = imag ( c1 )
!
!  Lower left.
!
  x = r4_floor ( real ( c1 ) )
  y = r4_floor ( imag ( c1 ) )
  r = ( x - xc )**2 + ( y - yc )**2
  r_min = r
  x_min = x
  y_min = y
!
!  Lower right.
!
  x = r4_floor ( real ( c1 ) ) + 1.0D+00
  y = r4_floor ( imag ( c1 ) )
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if
!
!  Upper right.
!
  x = r4_floor ( real ( c1 ) ) + 1.0D+00
  y = r4_floor ( imag ( c1 ) ) + 1.0D+00
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if
!
!  Upper left.
!
  x = r4_floor ( real ( c1 ) )
  y = r4_floor ( imag ( c1 ) ) + 1.0D+00
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if

  c4_nint = cmplx ( x_min, y_min, kind = 4 )

  return
end
function c4_norm_l1 ( x )

!*****************************************************************************80
!
!! C4_NORM_L1 evaluates the L1 norm of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    Numbers of equal norm lie along diamonds centered at (0,0).
!
!    The L1 norm can be defined here as:
!
!      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 4 ) C4_NORM_L1, the norm of X.
!
  implicit none

  real ( kind = 4 ) c4_norm_l1
  complex ( kind = 4 ) x

  c4_norm_l1 = abs ( real ( x, kind = 4 ) ) + abs ( aimag ( x ) )

  return
end
function c4_norm_l2 ( x )

!*****************************************************************************80
!
!! C4_NORM_L2 evaluates the L2 norm of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    Numbers of equal norm lie on circles centered at (0,0).
!
!    The L2 norm can be defined here as:
!
!      C4_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( imag ( X ) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 4 ) C4_NORM_L2, the 2-norm of X.
!
  implicit none

  real ( kind = 4 ) c4_norm_l2
  complex ( kind = 4 ) x

  c4_norm_l2 = sqrt ( ( real ( x, kind = 4 ) )**2 &
                   + ( aimag ( x ) )**2 )

  return
end
function c4_norm_li ( x )

!*****************************************************************************80
!
!! C4_NORM_LI evaluates the L-infinity norm of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    Numbers of equal norm lie along squares whose centers are at (0,0).
!
!    The L-infinity norm can be defined here as:
!
!      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 4 ) C4_NORM_LI, the infinity norm of X.
!
  implicit none

  real ( kind = 4 ) c4_norm_li
  complex ( kind = 4 ) x

  c4_norm_li = max ( abs ( real ( x, kind = 4 ) ), abs ( aimag ( x ) ) )

  return
end
function c4_normal_01 ( seed )

!*****************************************************************************80
!
!! C4_NORMAL_01 returns a unit pseudonormal C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 4 ) C4_NORMAL_01, a unit pseudornormal value.
!
  implicit none

  complex ( kind = 4 ) c4_normal_01
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2
  real ( kind = 4 ) x_c
  real ( kind = 4 ) x_r

  v1 = r4_uniform_01 ( seed )
  v2 = r4_uniform_01 ( seed )

  x_r = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * r4_pi * v2 )
  x_c = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * r4_pi * v2 )

  c4_normal_01 = cmplx ( x_r, x_c, kind = 4 )

  return
end
function c4_one ( )

!*****************************************************************************80
!
!! C4_ONE returns the value of 1 as a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 4 ) C4_ONE, the value of complex 1.
!
  implicit none

  complex ( kind = 4 ) c4_one

  c4_one = cmplx ( 1.0E+00, 0.0E+00, kind = 4 )

  return
end
subroutine c4_print ( a, title )

!*****************************************************************************80
!
!! C4_PRINT prints a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) A, the value to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  complex ( kind = 4 ) a
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a,2x,a,g14.6,a,g14.6,a)' ) &
      trim ( title ), '(', real ( a ), ',', imag ( a ), ')'
  else
    write ( *, '(a,g14.6,a,g14.6,a)' ) &
      '(', real ( a ), ',', imag ( a ), ')'
  end if

  return
end
function c4_real ( z )

!*****************************************************************************80
!
!! C4_REAL evaluates the real part of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the real part of a C4 with the REAL function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) C4_REAL, the function value.
!
  implicit none

  real ( kind = 4 ) c4_real
  complex ( kind = 4 ) z

  c4_real = real ( z )

  return
end
function c4_sin ( z )

!*****************************************************************************80
!
!! C4_SIN evaluates the sine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the sine of a C4 with the SIN function.
!
!    We use the relationship:
!
!      C4_SIN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_SIN, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) c4_sin
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00 )

  c4_sin = - c4_i * ( c4_exp ( c4_i * z ) - c4_exp ( - c4_i * z ) ) / 2.0E+00

  return
end
function c4_sinh ( z )

!*****************************************************************************80
!
!! C4_SINH evaluates the hyperbolic sine of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the hyperbolic sine of a C4.
!
!    We use the relationship:
!
!      C4_SINH ( C ) = ( C4_EXP ( C ) - C4_EXP ( - C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_SINH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) c4_sinh
  complex ( kind = 4 ) z

  c4_sinh = ( c4_exp ( z ) - c4_exp ( - z ) ) / 2.0E+00

  return
end
function c4_sqrt ( x )

!*****************************************************************************80
!
!! C4_SQRT returns the principal square root of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 supports the square root of a C4 with the SQRT function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) X, the argument.
!
!    Output, complex ( kind = 4 ) C4_SQRT, the function value.
!
  implicit none

  real ( kind = 4 ) arg
  real ( kind = 4 ) c4_arg
  real ( kind = 4 ) c4_mag
  complex ( kind = 4 ) c4_sqrt
  real ( kind = 4 ) mag
  complex ( kind = 4 ) x

  arg = c4_arg ( x )
  mag = c4_mag ( x )

  if ( mag == 0.0E+00 ) then

    c4_sqrt = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  else

    c4_sqrt = sqrt ( mag ) &
      * cmplx ( cos ( arg / 2.0E+00 ), &
                sin ( arg / 2.0E+00 ), kind = 4 )

  end if

  return
end
function c4_sub ( z1, z2 )

!*****************************************************************************80
!
!! C4_SUB subtracts two C4's.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 directly supports C4 subtraction with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z1, Z2, the values to subtract.
!
!    Output, complex ( kind = 4 ) C4_SUB, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_sub
  complex ( kind = 4 ) z1
  complex ( kind = 4 ) z2

  c4_sub = z1 - z2

  return
end
subroutine c4_swap ( x, y )

!*****************************************************************************80
!
!! C4_SWAP swaps two C4's.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex ( kind = 4 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex ( kind = 4 ) x
  complex ( kind = 4 ) y
  complex ( kind = 4 ) z

  z = x
  x = y
  y = z

  return
end
function c4_tan ( z )

!*****************************************************************************80
!
!! C4_TAN evaluates the tangent of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the tangent of a C4.
!
!    We use the relationship:
!
!      C4_TAN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) 
!                         / ( C4_EXP ( I * C ) + C4_EXP ( - i * C ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_TAN, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) c4_i
  complex ( kind = 4 ) c4_tan
  complex ( kind = 4 ) z

  c4_i = cmplx ( 0.0E+00, 1.0E+00 )

  c4_tan =  - c4_i * ( c4_exp ( c4_i * z ) - c4_exp ( - c4_i * z ) ) &
         /           ( c4_exp ( c4_i * z ) + c4_exp ( - c4_i * z ) )

  return
end
function c4_tanh ( z )

!*****************************************************************************80
!
!! C4_TANH evaluates the hyperbolic tangent of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    FORTRAN90 does not support the hyperbolic tangent of a C4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, complex ( kind = 4 ) C4_TANH, the function value.
!
  implicit none

  complex ( kind = 4 ) c4_exp
  complex ( kind = 4 ) c4_tanh
  complex ( kind = 4 ) z

  c4_tanh = ( c4_exp ( z ) - c4_exp ( - z ) ) &
          / ( c4_exp ( z ) + c4_exp ( - z ) )

  return
end
subroutine c4_to_cartesian ( z, x, y )

!*****************************************************************************80
!
!! C4_TO_CARTESIAN converts a C4 to Cartesian form.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) X, Y, the Cartesian form.
!
  implicit none

  real ( kind = 4 ) x
  real ( kind = 4 ) y
  complex ( kind = 4 ) z

  x = real ( z )
  y = aimag ( z )

  return
end
subroutine c4_to_polar ( z, r, theta )

!*****************************************************************************80
!
!! C4_TO_POLAR converts a C4 to polar form.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 4 ) Z, the argument.
!
!    Output, real ( kind = 4 ) R, THETA, the polar form.
!
  implicit none

  real ( kind = 4 ) c4_arg
  real ( kind = 4 ) c4_mag
  real ( kind = 4 ) r
  real ( kind = 4 ) theta
  complex ( kind = 4 ) z

  r = c4_mag ( z )
  theta = c4_arg ( z )

  return
end
function c4_uniform_01 ( seed )

!*****************************************************************************80
!
!! C4_UNIFORM_01 returns a unit pseudorandom C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C4_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

  c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

  return
end
function c4_zero ( )

!*****************************************************************************80
!
!! C4_ZERO returns the value of 0 as a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 4 ) C4_ZERO, the value of complex 0.
!
  implicit none

  complex ( kind = 4 ) c4_zero

  c4_zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  return
end
subroutine c4mat_add ( m, n, alpha, a, beta, b, c )

!*****************************************************************************80
!
!! C4MAT_ADD combines two C4MAT's with complex scalar factors.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, complex ( kind = 4 ) ALPHA, the first scale factor.
!
!    Input, complex ( kind = 4 ) A(M,N), the first matrix.
!
!    Input, complex ( kind = 4 ) BETA, the second scale factor.
!
!    Input, complex ( kind = 4 ) B(M,N), the second matrix.
!
!    Output, complex ( kind = 4 ) C(M,N), the result matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  complex ( kind = 4 ) alpha
  complex ( kind = 4 ) b(m,n)
  complex ( kind = 4 ) beta
  complex ( kind = 4 ) c(m,n)

  c(1:m,1:n) = alpha * a(1:m,1:n) + beta * b(1:m,1:n)

  return
end
subroutine c4mat_add_r4 ( m, n, alpha, a, beta, b, c )

!*****************************************************************************80
!
!! C4MAT_ADD_R4 combines two C4MAT's with real scalar factors.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 4 ) ALPHA, the first scale factor.
!
!    Input, complex ( kind = 4 ) A(M,N), the first matrix.
!
!    Input, real ( kind = 4 ) BETA, the second scale factor.
!
!    Input, complex ( kind = 4 ) B(M,N), the second matrix.
!
!    Output, complex ( kind = 4 ) C(M,N), the result matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  real ( kind = 4 ) alpha
  complex ( kind = 4 ) b(m,n)
  real ( kind = 4 ) beta
  complex ( kind = 4 ) c(m,n)

  c(1:m,1:n) = alpha * a(1:m,1:n) + beta * b(1:m,1:n)

  return
end
subroutine c4mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! C4MAT_COPY copies a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Output, complex ( kind = 4 ) B(M,N), the copied matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  complex ( kind = 4 ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine c4mat_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! C4MAT_FSS factors and solves a system with multiple right hand sides.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2013
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
!    Input/output, complex ( kind = 4 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, complex ( kind = 4 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 4 ) piv
  complex ( kind = 4 ) row(n)
  complex ( kind = 4 ) t(nb)
  complex ( kind = 4 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0E+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C4MAT_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0E+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0E+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0E+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
subroutine c4mat_identity ( n, a )

!*****************************************************************************80
!
!! C4MAT_IDENTITY sets a C4MAT to the identity.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, complex ( kind = 4 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  do i = 1, n
    a(i,i) = cmplx ( 1.0E+00, 0.0E+00, kind = 4 )
  end do

  return
end
subroutine c4mat_indicator ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_INDICATOR returns the C4MAT indicator matrix.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2008
!
!  Author
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Output, complex ( kind = 4 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = cmplx ( i, j, kind = 4 )
    end do
  end do

  return
end
subroutine c4mat_minvm ( n1, n2, a, b, c )

!*****************************************************************************80
!
!! C4MAT_MINVM computes inverse(A) * B for C4MAT's.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the matrices.
!
!    Input, complex ( kind = 4 ) A(N1,N1), B(N1,N2), the matrices.
!
!    Output, complex ( kind = 4 ) C(N1,N2), the result, C = inverse(A) * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  complex ( kind = 4 ) a(n1,n1)
  complex ( kind = 4 ) alu(n1,n1)
  complex ( kind = 4 ) b(n1,n2)
  complex ( kind = 4 ) c(n1,n2)
  integer ( kind = 4 ) info

  alu(1:n1,1:n1) = a(1:n1,1:n1)
  c(1:n1,1:n2) = b(1:n1,1:n2)

  call c4mat_fss ( n1, alu, n2, c, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4MAT_MINVM - Fatal error!'
    write ( *, '(a)' ) '  The matrix A was numerically singular.'
    stop
  end if

  return
end
subroutine c4mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! C4MAT_MM multiplies two C4MAT's.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, define the orders of
!    the matrices.
!
!    Input, complex ( kind = 4 ) A(N1,N2), B(N2,N3), the matrix factors.
!
!    Output, complex ( kind = 4 ) C(N1,N3), the product matrix.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  complex ( kind = 4 ) a(n1,n2)
  complex ( kind = 4 ) b(n2,n3)
  complex ( kind = 4 ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( a(1:n1,1:n2), b(1:n2,1:n3) )

  return
end
subroutine c4mat_nint ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_NINT rounds the entries of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, complex ( kind = 4 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = cmplx ( nint ( real ( a(i,j) ) ), &
                       nint ( imag ( a(i,j) ) ), kind = 4 )
    end do
  end do

  return
end
function c4mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_NORM_FRO returns the Frobenius norm of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      C4MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      c4vec_norm_l2 ( A * x ) <= c4mat_norm_fro ( A ) * c4vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!    Output, real ( kind = 4 ) C4MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) c4mat_norm_fro

  c4mat_norm_fro = &
    sqrt ( &
      sum ( &
        ( &
          abs ( a(1:m,1:n) ) &
        )**2 &
      ) &
    )

  return
end
function c4mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_NORM_L1 returns the matrix L1 norm of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!    The matrix L1 norm is defined as:
!
!      C4MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      c4vec_norm_l1 ( A * x ) <= c4mat_norm_l1 ( A ) * c4vec_norm_l1 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix whose L1 norm is desired.
!
!    Output, real ( kind = 4 ) C4MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  real ( kind = 4 ) c4mat_norm_l1
  real ( kind = 4 ) col_sum
  integer ( kind = 4 ) j

  c4mat_norm_l1 = 0.0E+00

  do j = 1, n
    col_sum = sum ( abs ( a(1:m,j) ) )
    c4mat_norm_l1 = max ( c4mat_norm_l1, col_sum )
  end do

  return
end
function c4mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_NORM_LI returns the matrix L-oo norm of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
!
!    The matrix L-oo norm is defined as:
!
!      C4MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-oo norm is derived from the vector L-oo norm,
!    and satisifies:
!
!      c4vec_norm_li ( A * x ) <= c4mat_norm_li ( A ) * c4vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix whose L-oo
!    norm is desired.
!
!    Output, real ( kind = 4 ) C4MAT_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  real ( kind = 4 ) c4mat_norm_li
  integer ( kind = 4 ) i
  real ( kind = 4 ) row_sum

  c4mat_norm_li = 0.0E+00

  do i = 1, m
    row_sum = sum ( abs ( a(i,1:n) ) )
    c4mat_norm_li = max ( c4mat_norm_li, row_sum )
  end do

  return
end
subroutine c4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! C4MAT_PRINT prints a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call c4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C4MAT_PRINT_SOME prints some of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = 4 ) zero

  zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( imag ( a(i,j) ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 4 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c4mat_scale ( m, n, alpha, a )

!*****************************************************************************80
!
!! C4MAT_SCALE scales a C4MAT by a complex scalar factor.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, complex ( kind = 4 ) ALPHA, the scale factor.
!
!    Input/output, complex ( kind = 4 ) A(M,N), the matrix to be scaled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  complex ( kind = 4 ) alpha

  a(1:m,1:n) = alpha * a(1:m,1:n)

  return
end
subroutine c4mat_scale_r8 ( m, n, alpha, a )

!*****************************************************************************80
!
!! C4MAT_SCALE_R8 scales a C4MAT by a real scalar factor.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 4 ) ALPHA, the scale factor.
!
!    Input/output, complex ( kind = 4 ) A(M,N), the matrix to be scaled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  real ( kind = 4 ) alpha

  a(1:m,1:n) = alpha * a(1:m,1:n)

  return
end
subroutine c4mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

    end do

  end do

  return
end
subroutine c4mat_zero ( m, n, a )

!*****************************************************************************80
!
!! C4MAT_ZERO zeroes a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, complex ( kind = 4 ) A(M,N), the zeroed matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)

  a(1:m,1:n) = 0.0E+00

  return
end
subroutine c4vec_copy ( n, a, b )

!*****************************************************************************80
!
!! C4VEC_COPY copies a C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, complex ( kind = 4 ) A(N), the array to be copied.
!
!    Output, complex ( kind = 4 ) B(N), the copy of the array.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  complex ( kind = 4 ) b(n)

  b(1:n) = a(1:n)

  return
end
subroutine c4vec_indicator ( n, a )

!*****************************************************************************80
!
!! C4VEC_INDICATOR sets a C4VEC to the indicator vector.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = 4 )
  end do

  return
end
subroutine c4vec_nint ( n, a )

!*****************************************************************************80
!
!! C4VEC_NINT rounds the entries of a C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, complex ( kind = 4 ) A(N), the vector to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( nint ( real ( a(i) ) ), &
                   nint ( imag ( a(i) ) ), kind = 4 )
  end do

  return
end
function c4vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! C4VEC_NORM_L2 returns the L2 norm of a C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The vector L2 norm is defined as:
!
!      C4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, complex ( kind = 4 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 4 ) C4VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  real ( kind = 4 ) c4vec_norm_l2

  c4vec_norm_l2 = sqrt ( sum ( conjg ( a(1:n) ) * a(1:n) ) )

  return
end
function c4vec_norm_squared ( n, a )

!*****************************************************************************80
!
!! C4VEC_NORM_SQUARED returns the square of the L2 norm of a C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The square of the vector L2 norm is defined as:
!
!      C4VEC_NORM_SQUARED = sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, complex ( kind = 4 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 4 ) C4VEC_NORM_SQUARED, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  real ( kind = 4 ) c4vec_norm_squared

  c4vec_norm_squared = sum ( conjg ( a(1:n) ) * a(1:n) )

  return
end
subroutine c4vec_print ( n, a, title )

!*****************************************************************************80
!
!! C4VEC_PRINT prints a C4VEC, with an optional title.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a1,1x,2g14.6)' ) i, ':', a(i)
  end do

  return
end
subroutine c4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C4VEC_PRINT_PART prints "part" of a C4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c4vec_print_some ( n, x, i_lo, i_hi, title )

!*****************************************************************************80
!
!! C4VEC_PRINT_SOME prints some of a C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 4 ) X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last entries
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( 1, i_lo ), min ( n, i_hi )
    write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
  end do

  return
end
subroutine c4vec_sort_a_l1 ( n, x )

!*****************************************************************************80
!
!! C4VEC_SORT_A_L1 ascending sorts a C4VEC by L1 norm.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The L1 norm of A+Bi is abs(A) + abs(B).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 4 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c4_le_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 4 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c4_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c4_le_l1 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c4vec_sort_a_l2 ( n, x )

!*****************************************************************************80
!
!! C4VEC_SORT_A_L2 ascending sorts a C4VEC by L2 norm.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 4 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c4_le_l2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 4 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c4_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c4_le_l2 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c4vec_sort_a_li ( n, x )

!*****************************************************************************80
!
!! C4VEC_SORT_A_LI ascending sorts a C4VEC by L-infinity norm.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The L infinity norm of A+Bi is max ( abs ( A ), abs ( B ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 4 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c4_le_li
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 4 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c4_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c4_le_li ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c4vec_spiral ( n, m, c1, c2, c )

!*****************************************************************************80
!
!! C4VEC_SPIRAL returns N points on a spiral between C1 and C2.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
!    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
!    
!    Then the polar form of the I-th point C(I) is:
!
!      R(I) = ( ( N - I     ) * R1 
!             + (     I - 1 ) * R2 ) 
!              / ( N    - 1 )
!
!      T(I) = ( ( N - I     ) * T1 
!             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
!             / ( N     - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on the spiral.
!
!    Input, integer ( kind = 4 ) M, the number of full circuits the 
!    spiral makes.
!
!    Input, complex ( kind = 4 ) C1, C2, the first and last points 
!    on the spiral.
!
!    Output, complex ( kind = 4 ) C(N), the points.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(n)
  complex ( kind = 4 ) c1
  complex ( kind = 4 ) c2
  real ( kind = 4 ) c4_arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) ri
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) t1
  real ( kind = 4 ) t2
  real ( kind = 4 ) ti

  r1 = abs ( c1 )
  r2 = abs ( c2 )

  t1 = c4_arg ( c1 )
  t2 = c4_arg ( c2 )

  if ( m == 0 ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0E+00 * r4_pi
    end if

  else if ( 0 < m ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0E+00 * r4_pi
    end if

    t2 = t2 + real ( m, kind = 4 ) * 2.0E+00 * r4_pi

  else if ( m < 0 ) then

    if ( t1 < t2 ) then
      t2 = t2 - 2.0E+00 * r4_pi
    end if

    t2 = t2 - real ( m, kind = 4 ) * 2.0E+00 * r4_pi

  end if

  do i = 1, n

    ri = ( real ( n - i,     kind = 4 ) * r1 &
         + real (     i - 1, kind = 4 ) * r2 ) &
         / real ( n     - 1, kind = 4 )

    ti = ( real ( n - i,     kind = 4 ) * t1 &
         + real (     i - 1, kind = 4 ) * t2 ) &
         / real ( n     - 1, kind = 4 )

    call polar_to_c4 ( ri, ti, c(i) )

  end do

  return
end
subroutine c4vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

  end do

  return
end
subroutine c4vec_unity ( n, a )

!*****************************************************************************80
!
!! C4VEC_UNITY returns the N roots of unity.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 4 ) A(N), the N roots of unity.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) theta

  do i = 1, n
    theta = r4_pi * real ( 2 * ( i - 1 ), kind = 4 ) / real ( n, kind = 4 )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = 4 )
  end do

  return
end
subroutine cartesian_to_c4 ( x, y, z )

!*****************************************************************************80
!
!! CARTESIAN_TO_C4 converts a Cartesian form to a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, Y, the Cartesian form.
!
!    Output, complex ( kind = 4 ) Z, the complex number.
!
  implicit none

  real ( kind = 4 ) x
  real ( kind = 4 ) y
  complex ( kind = 4 ) z

  z = cmplx ( x, y )

  return
end
subroutine polar_to_c4 ( r, theta, z )

!*****************************************************************************80
!
!! POLAR_TO_C4 converts a polar form to a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) R, THETA, the polar form.
!
!    Output, complex ( kind = 4 ) Z, the complex number.
!
  implicit none

  real ( kind = 4 ) r
  real ( kind = 4 ) theta
  complex ( kind = 4 ) z

  z = r * cmplx ( cos ( theta ), sin ( theta ) )

  return
end
function r4_csqrt ( x )

!*****************************************************************************80
!
!! R4_CSQRT returns the complex square root of an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number whose square root is desired.
!
!    Output, complex ( kind = 4 ) R4_CSQRT, the square root of X:
!
  implicit none

  real ( kind = 4 ) argument
  real ( kind = 4 ) magnitude
  real ( kind = 4 ), parameter :: pi = 3.141592653589793E+00
  complex ( kind = 4 ) r4_csqrt
  real ( kind = 4 ) x

  if ( 0.0E+00 < x ) then
    magnitude = x
    argument = 0.0E+00
  else if ( 0.0E+00 == x ) then
    magnitude = 0.0E+00
    argument = 0.0E+00
  else if ( x < 0.0E+00 ) then
    magnitude = -x
    argument = pi
  end if

  magnitude = sqrt ( magnitude )
  argument = argument / 2.0E+00

  r4_csqrt = magnitude * cmplx ( cos ( argument ), sin ( argument ), kind = 4 )

  return
end
function r4_floor ( r )

!*****************************************************************************80
!
!! R4_FLOOR rounds an R4 "down" (towards -oo) to the nearest integral R4.
!
!  Example:
!
!    R     Value
!
!   -1.1  -2.0
!   -1.0  -1.0
!   -0.9  -1.0
!    0.0   0.0
!    5.0   5.0
!    5.1   5.0
!    5.9   5.0
!    6.0   6.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) R, the value to be rounded down.
!
!    Output, real ( kind = 4 ) R4_FLOOR, the rounded value.
!
  implicit none

  real ( kind = 4 ) r
  real ( kind = 4 ) r4_floor
  real ( kind = 4 ) value

  value = real ( int ( r ), kind = 4 )
  if ( r < value ) then
    value = value - 1.0E+00
  end if

  r4_floor = value

  return
end
function r4_log_2 ( x )

!*****************************************************************************80
!
!! R4_LOG_2 returns the logarithm base 2 of an R4.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 4 ) R4_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2^R4_LOG_2.
!
  implicit none

  real ( kind = 4 ) r4_log_2
  real ( kind = 4 ) x

  if ( x == 0.0E+00 ) then
    r4_log_2 = - huge ( x )
  else
    r4_log_2 = log ( abs ( x ) ) / log ( 2.0E+00 )
  end if

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 4 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
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
!    13 December 2008
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r4poly2_root ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R4POLY2_ROOT returns the two roots of a quadratic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X * X + B * X + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, C, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 4 ) R1, R2, the roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  complex ( kind = 4 ) disc
  complex ( kind = 4 ) q
  complex ( kind = 4 ) r1
  complex ( kind = 4 ) r2

  if ( a == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4POLY2_ROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0E+00 * a * c
  q = -0.5E+00 * ( b + sign ( 1.0E+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine r4poly3_root ( a, b, c, d, r1, r2, r3 )

!*****************************************************************************80
!
!! R4POLY3_ROOT returns the three roots of a cubic polynomial.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A * X^3 + B * X^2 + C * X + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 4 ) R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) d
  complex ( kind = 4 ) i
  complex ( kind = 4 ) one
  real ( kind = 4 ), parameter :: pi = 3.141592653589793E+00
  real ( kind = 4 ) q
  real ( kind = 4 ) r
  complex ( kind = 4 ) r1
  complex ( kind = 4 ) r2
  complex ( kind = 4 ) r3
  real ( kind = 4 ) s1
  real ( kind = 4 ) s2
  real ( kind = 4 ) temp
  real ( kind = 4 ) theta

  if ( a == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4POLY3_ROOT - Fatal error!'
    write ( *, '(a)' ) '  A must not be zero!'
    stop
  end if

  one = cmplx ( 1.0d+00, 0.0E+00, kind = 4 )
  i = sqrt ( -one )

  q = ( ( b / a )**2 - 3.0E+00 * ( c / a ) ) / 9.0E+00

  r = ( 2.0E+00 * ( b / a )**3 - 9.0E+00 * ( b / a ) * ( c / a ) &
      + 27.0E+00 * ( d / a ) ) / 54.0E+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q**3 ) )
    r1 = -2.0E+00 * sqrt ( q ) * cos (   theta                  / 3.0E+00 )
    r2 = -2.0E+00 * sqrt ( q ) * cos ( ( theta + 2.0E+00 * pi ) / 3.0E+00 )
    r3 = -2.0E+00 * sqrt ( q ) * cos ( ( theta + 4.0E+00 * pi ) / 3.0E+00 )

  else if ( q * q * q <= r * r ) then

    temp = -r + sqrt ( r**2 - q**3 )
    s1 = sign ( 1.0E+00, temp ) * ( abs ( temp ) )**(1.0E+00/3.0E+00)

    temp = -r - sqrt ( r**2 - q**3 )
    s2 = sign ( 1.0E+00, temp ) * ( abs ( temp ) )**(1.0E+00/3.0E+00)

    r1 = s1 + s2
    r2 = -0.5E+00 * ( s1 + s2 ) + i * 0.5E+00 * sqrt ( 3.0E+00 ) * ( s1 - s2 )
    r3 = -0.5E+00 * ( s1 + s2 ) - i * 0.5E+00 * sqrt ( 3.0E+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0E+00 * a )
  r2 = r2 - b / ( 3.0E+00 * a )
  r3 = r3 - b / ( 3.0E+00 * a )

  return
end
subroutine r4poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

!*****************************************************************************80
!
!! R4POLY4_ROOT returns the four roots of a quartic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 4 ) R1, R2, R3, R4, the roots of the polynomial.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) a3
  real ( kind = 4 ) a4
  real ( kind = 4 ) b
  real ( kind = 4 ) b3
  real ( kind = 4 ) b4
  real ( kind = 4 ) c
  real ( kind = 4 ) c3
  real ( kind = 4 ) c4
  real ( kind = 4 ) d
  real ( kind = 4 ) d3
  real ( kind = 4 ) d4
  real ( kind = 4 ) e
  complex ( kind = 4 ) p
  complex ( kind = 4 ) q
  complex ( kind = 4 ) r
  complex ( kind = 4 ) r1
  complex ( kind = 4 ) r2
  complex ( kind = 4 ) r3
  complex ( kind = 4 ) r4
  complex ( kind = 4 ) zero

  zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  if ( a == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4POLY4_ROOT - Fatal error!'
    write ( *, '(a)') '  A must not be zero!'
    stop
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0E+00
  b3 = -b4
  c3 = a4 * c4 - 4.0E+00 * d4
  d3 = -a4 * a4 * d4 + 4.0E+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call r4poly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25E+00 * A4**2 - B4 + R1 )
!
  r = sqrt ( 0.25E+00 * a4**2 - b4 + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75E+00 * a4**2 - r**2 - 2.0E+00 * b4 &
        + 0.25E+00 * ( 4.0E+00 * a4 * b4 - 8.0E+00 * c4 - a4**3 ) / r )

    q = sqrt ( 0.75E+00 * a4**2 - r**2 - 2.0E+00 * b4 &
        - 0.25E+00 * ( 4.0E+00 * a4 * b4 - 8.0E+00 * c4 - a4**3 ) / r )

  else

    p = sqrt ( 0.75E+00 * a4**2 - 2.0E+00 * b4 &
      + 2.0E+00 * sqrt ( r1**2 - 4.0E+00 * d4 ) )

    q = sqrt ( 0.75E+00 * a4**2 - 2.0E+00 * b4 &
      - 2.0E+00 * sqrt ( r1**2 - 4.0E+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = -0.25E+00 * a4 + 0.5E+00 * r + 0.5E+00 * p
  r2 = -0.25E+00 * a4 + 0.5E+00 * r - 0.5E+00 * p
  r3 = -0.25E+00 * a4 - 0.5E+00 * r + 0.5E+00 * q
  r4 = -0.25E+00 * a4 - 0.5E+00 * r - 0.5E+00 * q

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
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
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
!    06 August 2005
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
