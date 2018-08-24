function anorm ( arg )

!*****************************************************************************80
!
!! ANORM evaluates the normal distribution function.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!      p(x) = 1/sqrt(2*pi) * integral ( -oo < t < x ) exp(-t^2/2) dt
!    The main computation evaluates near-minimax approximations
!    derived from those in "Rational Chebyshev approximations for
!    the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   XMIN  = the smallest positive floating-point number.
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   EPS   = argument below which anorm(x) may be represented by
!       0.5  and above which  x*x  will not underflow.
!       A conservative value is the largest machine number X
!       such that   1.0 + X = 1.0   to machine precision.
!   XLOW  = the most negative argument for which ANORM does not
!       vanish.  This is the negative of the solution to
!        W(x) * (1-1/x**2) = XMIN,
!       where W(x) = exp(-x*x/2)/[x*sqrt(2*pi)].
!   XUPPR = positive argument beyond which anorm = 1.0.  A
!       conservative value is the solution to the equation
!        exp(-x*x/2) = EPS,
!       i.e., XUPPR = sqrt[-2 ln(eps)].
!
! Error returns
!
!  The program returns  ANORM = 0     for  ARG <= XLOW.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) &
    a,anorm,arg,b,c,d,del,eps,half,p,one,q,result,sixten, &
    sqrpi,thrsh,root32,x,xlow,xden,xnum,y,xsq,xuppr,zero
  dimension a(5),b(4),c(9),d(8),p(6),q(5)
!
!  Mathematical constants
!
!  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
!  THRSH is the argument for which anorm = 0.75.
!
  data one,half,zero,sixten/1.0d0,0.5d0,0.0d0,1.60d1/, &
    sqrpi/3.9894228040143267794d-1/,thrsh/0.66291d0/, &
    root32/5.656854248d0/
!
!  Machine-dependent constants
!
  data eps/1.11d-16/,xlow/-37.519d0/,xuppr/8.572d0/
!
!  Coefficients for approximation in first interval
!
  data a/2.2352520354606839287d00,1.6102823106855587881d02, &
     1.0676894854603709582d03,1.8154981253343561249d04, &
     6.5682337918207449113d-2/
  data b/4.7202581904688241870d01,9.7609855173777669322d02, &
     1.0260932208618978205d04,4.5507789335026729956d04/
!
!  Coefficients for approximation in second interval
!
  data c/3.9894151208813466764d-1,8.8831497943883759412d00, &
     9.3506656132177855979d01,5.9727027639480026226d02, &
     2.4945375852903726711d03,6.8481904505362823326d03, &
     1.1602651437647350124d04,9.8427148383839780218d03, &
     1.0765576773720192317d-8/
  data d/2.2266688044328115691d01,2.3538790178262499861d02, &
     1.5193775994075548050d03,6.4855582982667607550d03, &
     1.8615571640885098091d04,3.4900952721145977266d04, &
     3.8912003286093271411d04,1.9685429676859990727d04/
!
!  Coefficients for approximation in third interval
!
  data p/2.1589853405795699d-1,1.274011611602473639d-1, &
     2.2235277870649807d-2,1.421619193227893466d-3, &
     2.9112874951168792d-5,2.307344176494017303d-2/
  data q/1.28426009614491121d00,4.68238212480865118d-1, &
     6.59881378689285515d-2,3.78239633202758244d-3, &
     7.29751555083966205d-5/

  x = arg
  y = abs(x)

  if (y <= thrsh) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    xsq = zero
    if (y > eps) xsq = x * x
    xnum = a(5)*xsq
    xden = xsq
    do i = 1, 3
      xnum = (xnum + a(i)) * xsq
      xden = (xden + b(i)) * xsq
    end do
    result = x * (xnum + a(4)) / (xden + b(4))
    result = half + result
!
!  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
!
  else if (y <= root32) then

    xnum = c(9)*y
    xden = y
    do i = 1, 7
      xnum = (xnum + c(i)) * y
      xden = (xden + d(i)) * y
    end do
    result = (xnum + c(8)) / (xden + d(8))
    xsq = aint(y*sixten)/sixten
    del = (y-xsq)*(y+xsq)
    result = exp(-xsq*xsq*half)*exp(-del*half)*result
    if (x > zero) result = one - result
!
!  Evaluate  anorm  for |X| > sqrt(32)
!
  else

    result = zero

    if ((x >= xlow) .and. (x < xuppr)) then
      xsq = one / (x * x)
      xnum = p(6)*xsq
      xden = xsq
      do i = 1, 4
        xnum = (xnum + p(i)) * xsq
        xden = (xden + q(i)) * xsq
      end do
      result = xsq *(xnum + p(5)) / (xden + q(5))
      result = (sqrpi -  result) / y
      xsq = aint(x*sixten)/sixten
      del = (x-xsq)*(x+xsq)
      result = exp(-xsq*xsq*half)*exp(-del*half)*result
    end if

    if (x > zero) result = one - result

  end if
!
!  Fix up for negative argument, erf, etc.
!
  anorm = result

  return
end
function besei0 ( x )

!*****************************************************************************80
!
!! BESEI0 evaluates the exponentially scaled Bessel I function of order 0.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the first kind of order zero
!    multiplied by EXP(-ABS(X)), where EXP is the
!    exponential function, ABS is the absolute value, and X
!    is any argument.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besei0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 2
  call calci0 ( x, result, jint )
  besei0 = result

  return
end
function besei1 ( x )

!*****************************************************************************80
!
!! BESEI1 evaluates the exponentially scaled Bessel I function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the first kind of order one
!    multiplied by EXP(-ABS(X)), where EXP is the
!    exponential function, ABS is the absolute value, and X
!    is any argument.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besei1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint=2
  call calci1(x,result,jint)
  besei1=result

  return
end
function besek0 ( x )

!*****************************************************************************80
!
!! BESEK0 evaluates the exponentially scaled Bessel K function of order 0.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    multiplied by the Exponential function, for arguments
!    0.0 < ARG.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody and Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besek0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 2
  call calck0(x,result,jint)
  besek0 = result

  return
end
function besek1 ( x )

!*****************************************************************************80
!
!! BESEK1 evaluates the exponentially scaled Bessel K function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the second kind of order one
!    multiplied by the exponential function, for arguments
!    XLEAST <= ARG <= XMAX.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besek1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 2
  call calck1(x,result,jint)
  besek1 = result

  return
end
function besi0 ( x )

!*****************************************************************************80
!
!! BESI0 evaluates the modified Bessel I function of order 0.
!
!  Discussion:
!
!    This function computes approximate values for
!    modified Bessel functions of the first kind of order zero for
!    arguments ABS(ARG) <= XMAX  (see comments heading CALCI0).
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besi0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 1
  call calci0 ( x, result, jint )
  besi0 = result

  return
end
function besi1 ( x )

!*****************************************************************************80
!
!! BESI1 evaluates the modified Bessel I function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for
!    modified Bessel functions of the first kind of order one for
!    arguments ABS(ARG) <= XMAX  (see comments heading CALCI1).
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besi1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint=1
  call calci1(x,result,jint)
  besi1=result

  return
end
function besj0 ( x )

!*****************************************************************************80
!
!! BESJ0 evaluates the Bessel J function of order 0.
!
!  Discussion:
!
!    This function computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX.
!    See the comments heading CALJY0.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besj0
  integer ( kind = 4 ) jint
  real ( kind = 8 )  x, result

  jint=0
  call caljy0(x,result,jint)
  besj0 = result

  return
end
function besj1 ( x )

!*****************************************************************************80
!
!! BESJ1 evaluates the Bessel J function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX.
!    See the comments heading CALJY1.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besj1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result,x

  jint=0
  call caljy1(x,result,jint)
  besj1 = result

  return
end
function besk0 ( x )

!*****************************************************************************80
!
!! BESK0 evaluates the modified Bessel K function of order 0.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    for arguments 0.0 < ARG <= XMAX.  See comments heading
!    CALCK0.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besk0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 1
  call calck0 ( x, result, jint )
  besk0 = result

  return
end
function besk1 ( x )

!*****************************************************************************80
!
!! BESK1 evaluates the modified Bessel K function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for the
!    modified Bessel function of the second kind of order one
!    for arguments XLEAST <= ARG <= XMAX.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besk1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 1
  call calck1(x,result,jint)
  besk1 = result

  return
end
function besy0 ( x )

!*****************************************************************************80
!
!! BESY0 evaluates the Bessel Y function of order 0.
!
!  Discussion:
!
!    This subprogram computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!    See comments heading CALJY0.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besy0
  integer ( kind = 4 ) jint
  real ( kind = 8 )  x, result

  jint=1
  call caljy0(x,result,jint)
  besy0 = result

  return
end
function besy1 ( x )

!*****************************************************************************80
!
!! BESY1 evaluates the Bessel Y function of order 1.
!
!  Discussion:
!
!    This function computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!    See comments heading CALJY1.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) besy1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result,x

  jint=1
  call caljy1(x,result,jint)
  besy1 = result

  return
end
subroutine calcei ( arg, result, int )

!*****************************************************************************80
!
!! CALCEI computes exponential integrals.
!
!  Discussion:
!
!    This routine computes the exponential integrals Ei(x),
!    E1(x), and  exp(-x)*Ei(x)  for real arguments  x  
!    where, if x > 0,
!      Ei(x) = integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
!    while, if x < 0,
!      Ei(x) = -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
!    and where the first integral is a principal value integral.
!
!    The packet contains three function type subprograms: EI, EONE,
!    and EXPEI;  and one subroutine type subprogram: CALCEI.  The
!    calling statements for the primary entries are:
!      Y = EI(X),    where  X /= 0,
!      Y = EONE(X),      where  X > 0,
!      Y = EXPEI(X),     where  X /= 0,
!
!    and where the entry points correspond to the functions Ei(x),
!    E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
!    is intended for internal packet use only, all computations within
!    the packet being concentrated in this routine.  The function
!    subprograms invoke CALCEI with the statement
!      CALL CALCEI(ARG,RESULT,INT)
!    where the parameter usage is as follows
!
!      Function      Parameters for CALCEI
!      Call           ARG         RESULT         INT
!
!      EI(X)          X /= 0      Ei(X)            1
!      EONE(X)        X > 0      -Ei(-X)           2
!      EXPEI(X)       X /= 0      exp(-X)*Ei(X)    3
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations published in Math. Comp. 22, 641-649 (1968), and
!    Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
!    transportable program is patterned after the machine-dependent
!    FUNPACK packet  NATSEI,  but cannot match that version for
!    efficiency or accuracy.  This version uses rational functions
!    that theoretically approximate the exponential integrals to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta = radix for the floating-point system.
!   minexp = smallest representable power of beta.
!   maxexp = smallest power of beta that overflows.
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XBIG = largest argument acceptable to EONE; solution to
!      equation:
!         exp(-x)/x * (1 + 1/x) = beta ** minexp.
!   XINF = largest positive machine number; approximately
!         beta ** maxexp
!   XMAX = largest argument acceptable to EI; solution to
!      equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
!
! Error returns
!
!  The following table shows the types of error that may be
!  encountered in this routine and the function value supplied
!  in each case.
!
!   Error   Argument   function values for
!        Range     EI  EXPEI     EONE
!
!     UNDERFLOW  (-)X > XBIG     0    -     0
!     OVERFLOW  X >= XMAX    XINF  -     -
!     ILLEGAL X   X = 0   -XINF    -XINF     XINF
!     ILLEGAL X  X < 0   -    -     USE ABS(X)
!
  implicit none

  integer ( kind = 4 ) i,int
  real ( kind = 8 ) &
    a,arg,b,c,d,exp40,e,ei,f,four,fourty,frac,half,one,p, &
    plg,px,p037,p1,p2,q,qlg,qx,q1,q2,r,result,s,six,sump, &
    sumq,t,three,twelve,two,two4,w,x,xbig,xinf,xmax,xmx0, &
    x0,x01,x02,x11,y,ysq,zero
  dimension  a(7),b(6),c(9),d(9),e(10),f(10),p(10),q(10),r(10), &
    s(9),p1(10),q1(9),p2(10),q2(9),plg(4),qlg(4),px(10),qx(10)
!
!  Mathematical constants
!  EXP40 = exp(40)
!  X0 = zero of Ei
!  X01/X11 + X02 = zero of Ei to extra precision
!
  data zero,p037,half,one,two/0.0d0,0.037d0,0.5d0,1.0d0,2.0d0/, &
       three,four,six,twelve,two4/3.0d0,4.0d0,6.0d0,12.d0,24.0d0/, &
       fourty,exp40/40.0d0,2.3538526683701998541d17/, &
       x01,x11,x02/381.5d0,1024.0d0,-5.1182968633365538008d-5/, &
       x0/3.7250741078136663466d-1/
!
!  Machine-dependent constants
!
  data xinf/1.79d+308/,xmax/716.351d0/,xbig/701.84d0/
!
! Coefficients  for -1.0 <= X < 0.0
!
  data a/1.1669552669734461083368d2, 2.1500672908092918123209d3, &
     1.5924175980637303639884d4, 8.9904972007457256553251d4, &
     1.5026059476436982420737d5,-1.4815102102575750838086d5, &
     5.0196785185439843791020d0/
  data b/4.0205465640027706061433d1, 7.5043163907103936624165d2, &
     8.1258035174768735759855d3, 5.2440529172056355429883d4, &
     1.8434070063353677359298d5, 2.5666493484897117319268d5/
!
! Coefficients for -4.0 <= X < -1.0
!
  data c/3.828573121022477169108d-1, 1.107326627786831743809d+1, &
     7.246689782858597021199d+1, 1.700632978311516129328d+2, &
     1.698106763764238382705d+2, 7.633628843705946890896d+1, &
     1.487967702840464066613d+1, 9.999989642347613068437d-1, &
     1.737331760720576030932d-8/
  data d/8.258160008564488034698d-2, 4.344836335509282083360d+0, &
     4.662179610356861756812d+1, 1.775728186717289799677d+2, &
     2.953136335677908517423d+2, 2.342573504717625153053d+2, &
     9.021658450529372642314d+1, 1.587964570758947927903d+1, &
     1.000000000000000000000d+0/
!
! Coefficients for X < -4.0
!
  data e/1.3276881505637444622987d+2,3.5846198743996904308695d+4, &
     1.7283375773777593926828d+5,2.6181454937205639647381d+5, &
     1.7503273087497081314708d+5,5.9346841538837119172356d+4, &
     1.0816852399095915622498d+4,1.0611777263550331766871d03, &
     5.2199632588522572481039d+1,9.9999999999999999087819d-1/
  data f/3.9147856245556345627078d+4,2.5989762083608489777411d+5, &
     5.5903756210022864003380d+5,5.4616842050691155735758d+5, &
     2.7858134710520842139357d+5,7.9231787945279043698718d+4, &
     1.2842808586627297365998d+4,1.1635769915320848035459d+3, &
     5.4199632588522559414924d+1,1.0d0/
!
!  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
       -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
       -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
! Coefficients for  0.0 < X < 6.0,
!  ratio of Chebyshev polynomials
!
  data p/-1.2963702602474830028590d01,-1.2831220659262000678155d03, &
     -1.4287072500197005777376d04,-1.4299841572091610380064d06, &
     -3.1398660864247265862050d05,-3.5377809694431133484800d08, &
      3.1984354235237738511048d08,-2.5301823984599019348858d10, &
      1.2177698136199594677580d10,-2.0829040666802497120940d11/
  data q/ 7.6886718750000000000000d01,-5.5648470543369082846819d03, &
      1.9418469440759880361415d05,-4.2648434812177161405483d06, &
      6.4698830956576428587653d07,-7.0108568774215954065376d08, &
      5.4229617984472955011862d09,-2.8986272696554495342658d10, &
      9.8900934262481749439886d10,-8.9673749185755048616855d10/
!
! J-fraction coefficients for 6.0 <= X < 12.0
!
  data r/-2.645677793077147237806d00,-2.378372882815725244124d00, &
     -2.421106956980653511550d01, 1.052976392459015155422d01, &
      1.945603779539281810439d01,-3.015761863840593359165d01, &
      1.120011024227297451523d01,-3.988850730390541057912d00, &
      9.565134591978630774217d00, 9.981193787537396413219d-1/
  data s/ 1.598517957704779356479d-4, 4.644185932583286942650d00, &
      3.697412299772985940785d02,-8.791401054875438925029d00, &
      7.608194509086645763123d02, 2.852397548119248700147d01, &
      4.731097187816050252967d02,-2.369210235636181001661d02, &
      1.249884822712447891440d00/
!
! J-fraction coefficients for 12.0 <= X < 24.0
!
  data p1/-1.647721172463463140042d00,-1.860092121726437582253d01, &
      -1.000641913989284829961d01,-2.105740799548040450394d01, &
      -9.134835699998742552432d-1,-3.323612579343962284333d01, &
       2.495487730402059440626d01, 2.652575818452799819855d01, &
      -1.845086232391278674524d00, 9.999933106160568739091d-1/
  data q1/ 9.792403599217290296840d01, 6.403800405352415551324d01, &
       5.994932325667407355255d01, 2.538819315630708031713d02, &
       4.429413178337928401161d01, 1.192832423968601006985d03, &
       1.991004470817742470726d02,-1.093556195391091143924d01, &
       1.001533852045342697818d00/
!
! J-fraction coefficients for  X >= 24.0
!
  data p2/ 1.75338801265465972390d02,-2.23127670777632409550d02, &
      -1.81949664929868906455d01,-2.79798528624305389340d01, &
      -7.63147701620253630855d00,-1.52856623636929636839d01, &
      -7.06810977895029358836d00,-5.00006640413131002475d00, &
      -3.00000000320981265753d00, 1.00000000000000485503d00/
  data q2/ 3.97845977167414720840d04, 3.97277109100414518365d00, &
       1.37790390235747998793d02, 1.17179220502086455287d02, &
       7.04831847180424675988d01,-1.20187763547154743238d01, &
      -7.99243595776339741065d00,-2.99999894040324959612d00, &
       1.99999999999048104167d00/

  x = arg

  if (x == zero) then

    ei = -xinf
    if (int == 2) ei = -ei

  else if ((x < zero) .or. (int == 2)) then
!
! Calculate EI for negative argument or for E1.
!
    y = abs(x)

    if (y <= one) then

      sump = a(7) * y + a(1)
      sumq = y + b(1)
      do i = 2, 6
        sump = sump * y + a(i)
        sumq = sumq * y + b(i)
      end do
      ei = log(y) - sump / sumq
      if (int == 3) ei = ei * exp(y)

    else if (y <= four) then

      w = one / y
      sump = c(1)
      sumq = d(1)
      do i = 2, 9
        sump = sump * w + c(i)
        sumq = sumq * w + d(i)
      end do
      ei = - sump / sumq
      if (int /= 3) ei = ei * exp(-y)

    else

      if ((y > xbig) .and. (int < 3)) then
        ei = zero
      else
        w = one / y
        sump = e(1)
        sumq = f(1)
        do i = 2, 10
          sump = sump * w + e(i)
          sumq = sumq * w + f(i)
        end do
        ei = -w * (one - w * sump / sumq )
        if (int /= 3) ei = ei * exp(-y)
      end if

    end if

    if (int == 2) ei = -ei

  else if (x < six) then
!
!  To improve conditioning, rational approximations are expressed
!  in terms of Chebyshev polynomials for 0 <= X < 6, and in
!  continued fraction form for larger X.
!
    t = x + x
    t = t / three - two
    px(1) = zero
    qx(1) = zero
    px(2) = p(1)
    qx(2) = q(1)
    do i = 2, 9
      px(i+1) = t * px(i) - px(i-1) + p(i)
      qx(i+1) = t * qx(i) - qx(i-1) + q(i)
    end do
    sump = half * t * px(10) - px(9) + p(10)
    sumq = half * t * qx(10) - qx(9) + q(10)
    frac = sump / sumq
    xmx0 = (x - x01/x11) - x02

    if (abs(xmx0) >= p037) then

      ei = log(x/x0) + xmx0 * frac

      if (int == 3) ei = exp(-x) * ei

    else
!
!  Special approximation to  ln(X/X0)  for X close to X0
!
      y = xmx0 / (x + x0)
      ysq = y*y
      sump = plg(1)
      sumq = ysq + qlg(1)
      do i = 2, 4
        sump = sump*ysq + plg(i)
        sumq = sumq*ysq + qlg(i)
      end do
      ei = (sump / (sumq*(x+x0)) + frac) * xmx0
      if (int == 3) ei = exp(-x) * ei

    end if

  else if (x < twelve) then

    frac = zero
    do i = 1, 9
      frac = s(i) / (r(i) + x + frac)
    end do
    ei = (r(10) + frac) / x
    if (int /= 3) ei = ei * exp(x)

  else if (x <= two4) then

    frac = zero
    do i = 1, 9
      frac = q1(i) / (p1(i) + x + frac)
    end do
    ei = (p1(10) + frac) / x
    if (int /= 3) ei = ei * exp(x)

  else

    if ((x >= xmax) .and. (int < 3)) then

      ei = xinf

    else

      y = one / x
      frac = zero
      do i = 1, 9
        frac = q2(i) / (p2(i) + x + frac)
      end do
      frac = p2(10) + frac
      ei = y + y * y * frac

      if (int /= 3) then
        if (x <= xmax-two4) then
          ei = ei * exp(x)
        else
!
!  Calculation reformulated to avoid premature overflow
!
          ei = (ei * exp(x-fourty)) * exp40
        end if

      end if

    end if

  end if

  result = ei

  return
end
subroutine calci0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCI0 evaluates modified Bessel I functions of order 0.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the first kind
!   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!   arguments X.  It contains two function type subprograms, BESI0
!   and BESEI0, and one subroutine type subprogram, CALCI0.
!   The calling statements for the primary entries are
!
!       Y=BESI0(X)
!   and
!       Y=BESEI0(X)
!
!   where the entry points correspond to the functions I0(X) and
!   EXP(-ABS(X))*I0(X), respectively.  The routine CALCI0 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCI0 with the statement
!      CALL CALCI0(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!function         Parameters for CALCI0
!   Call      ARG      RESULT      JINT
!
!     BESI0(ARG)    ABS(ARG) <= XMAX    I0(ARG)       1
!     BESEI0(ARG)    any real ARG    EXP(-ABS(ARG))*I0(ARG) 2
!
!   The main computation evaluates slightly modified forms of
!   minimax approximations generated by Blair and Edwards, Chalk
!   River (Atomic Energy of Canada Limited) Report AECL-4928,
!   October, 1974.  This transportable program is patterned after
!   the machine-dependent FUNPACK packet NATSI0, but cannot match
!   that version for efficiency or accuracy.  This version uses
!   rational functions that theoretically approximate I-SUB-0(X)
!   to at least 18 significant decimal digits.  The accuracy
!   achieved depends on the arithmetic system, the compiler, the
!   intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XSMALL = Positive argument such that 1.0 - X = 1.0 to
!    machine precision for all ABS(X) <= XSMALL.
!   XINF =   Largest positive machine number; approximately
!    beta**maxexp
!   XMAX =   Largest argument acceptable to BESI0;  Solution to
!    equation:
!       W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
!    where  W(X) = EXP(X)/SQRT(2*PI*X)
!
! Error returns
!
!  The program returns XINF for BESI0 for ABS(ARG) > XMAX.
!
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
     a,arg,b,exp40,forty,one,one5,p,pp,q,qq,result, &
     rec15,sump,sumq,two25,x,xinf,xmax,xsmall,xx
  dimension p(15),pp(8),q(5),qq(7)
!
!  Mathematical constants
!
  data one/1.0d0/,one5/15.0d0/,exp40/2.353852668370199854d17/, &
       forty/40.0d0/,rec15/6.6666666666666666666d-2/, &
       two25/225.0d0/
!
!  Machine-dependent constants
!
  data xsmall/5.55d-17/,xinf/1.79d308/,xmax/713.986d0/
!
!  Coefficients for XSMALL <= ABS(ARG) < 15.0
!
  data  p/-5.2487866627945699800d-18,-1.5982226675653184646d-14, &
      -2.6843448573468483278d-11,-3.0517226450451067446d-08, &
      -2.5172644670688975051d-05,-1.5453977791786851041d-02, &
      -7.0935347449210549190d+00,-2.4125195876041896775d+03, &
      -5.9545626019847898221d+05,-1.0313066708737980747d+08, &
      -1.1912746104985237192d+10,-8.4925101247114157499d+11, &
      -3.2940087627407749166d+13,-5.5050369673018427753d+14, &
      -2.2335582639474375249d+15/
  data  q/-3.7277560179962773046d+03, 6.5158506418655165707d+06, &
      -6.5626560740833869295d+09, 3.7604188704092954661d+12, &
      -9.7087946179594019126d+14/
!
!  Coefficients for 15.0 <= ABS(ARG)
!
  data pp/-3.9843750000000000000d-01, 2.9205384596336793945d+00, &
      -2.4708469169133954315d+00, 4.7914889422856814203d-01, &
      -3.7384991926068969150d-03,-2.6801520353328635310d-03, &
       9.9168777670983678974d-05,-2.1877128189032726730d-06/
  data qq/-3.1446690275135491500d+01, 8.5539563258012929600d+01, &
      -6.0228002066743340583d+01, 1.3982595353892851542d+01, &
      -1.1151759188741312645d+00, 3.2547697594819615062d-02, &
      -5.5194330231005480228d-04/

  x = abs(arg)

  if (x < xsmall) then

    result = one

  else if (x < one5) then
!
!  XSMALL <=  ABS(ARG)  < 15.0
!
    xx = x * x
    sump = p(1)
    do i = 2, 15
      sump = sump * xx + p(i)
    end do

    xx = xx - two25
    sumq = ((((xx+q(1))*xx+q(2))*xx+q(3))*xx+q(4))*xx+q(5)
    result = sump / sumq
    if (jint == 2) result = result * exp(-x)

  else if (x >= one5) then

    if ((jint == 1) .and. (x > xmax)) then

      result = xinf

    else
!
!  15.0  <=  ABS(ARG)
!
      xx = one / x - rec15
      sump = ((((((pp(1)*xx+pp(2))*xx+pp(3))*xx+pp(4))*xx+ &
        pp(5))*xx+pp(6))*xx+pp(7))*xx+pp(8)
      sumq = ((((((xx+qq(1))*xx+qq(2))*xx+qq(3))*xx+ &
        qq(4))*xx+qq(5))*xx+qq(6))*xx+qq(7)
      result = sump / sumq

      if (jint == 2) then

        result = (result - pp(1)) / sqrt(x)

      else
!
!  Calculation reformulated to avoid premature overflow
!
        if (x <=(xmax-one5)) then
          a = exp(x)
          b = one
        else
          a = exp(x-forty)
          b = exp40
        end if

        result = ((result*a-pp(1)*a)/sqrt(x))*b

      end if

    end if

  end if
!
!  Return for ABS(ARG) < XSMALL
!
  return
end
subroutine calci1 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCI1 evaluates modified Bessel I functions of order 1.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the first kind
!    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
!    arguments X.  It contains two function type subprograms, BESI1
!    and BESEI1, and one subroutine type subprogram, CALCI1.
!    The calling statements for the primary entries are
!
!       Y=BESI1(X)
!   and
!       Y=BESEI1(X)
!
!   where the entry points correspond to the functions I1(X) and
!   EXP(-ABS(X))*I1(X), respectively.  The routine CALCI1 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCI1 with the statement
!      CALL CALCI1(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!function         Parameters for CALCI1
!   Call      ARG      RESULT      JINT
!
!     BESI1(ARG)    ABS(ARG) <= XMAX    I1(ARG)       1
!     BESEI1(ARG)    any real ARG    EXP(-ABS(ARG))*I1(ARG) 2
!
!   The main computation evaluates slightly modified forms of
!   minimax approximations generated by Blair and Edwards, Chalk
!   River (Atomic Energy of Canada Limited) Report AECL-4928,
!   October, 1974.  This transportable program is patterned after
!   the machine-dependent FUNPACK packet NATSI1, but cannot match
!   that version for efficiency or accuracy.  This version uses
!   rational functions that theoretically approximate I-SUB-1(X)
!   to at least 18 significant decimal digits.  The accuracy
!   achieved depends on the arithmetic system, the compiler, the
!   intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XSMALL = Positive argument such that 1.0 - X = 1.0 to
!    machine precision for all ABS(X) <= XSMALL.
!   XINF =   Largest positive machine number; approximately
!    beta**maxexp
!   XMAX =   Largest argument acceptable to BESI1;  Solution to
!    equation:
!       EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
!
! Error returns
!
!  The program returns the value XINF for ABS(ARG) > XMAX.
!
  implicit none

  integer ( kind = 4 ) j,jint
  real ( kind = 8 ) &
      a,arg,b,exp40,forty,half,one,one5,p,pbar,pp,q,qq,rec15, &
      result,sump,sumq,two25,x,xinf,xmax,xsmall,xx,zero
  dimension p(15),pp(8),q(5),qq(6)
!
!  Mathematical constants
!
  data one/1.0d0/,one5/15.0d0/,exp40/2.353852668370199854d17/, &
       forty/40.0d0/,rec15/6.6666666666666666666d-2/, &
       two25/225.0d0/,half/0.5d0/,zero/0.0d0/
!
!  Machine-dependent constants
!
  data xsmall/5.55d-17/,xinf/1.79d308/,xmax/713.987d0/
!
!  Coefficients for XSMALL <= ABS(ARG) < 15.0
!
  data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16, &
     -1.1928788903603238754d-12,-1.4831904935994647675d-09, &
     -1.3466829827635152875d-06,-9.1746443287817501309d-04, &
     -4.7207090827310162436d-01,-1.8225946631657315931d+02, &
     -5.1894091982308017540d+04,-1.0588550724769347106d+07, &
     -1.4828267606612366099d+09,-1.3357437682275493024d+11, &
     -6.9876779648010090070d+12,-1.7732037840791591320d+14, &
     -1.4577180278143463643d+15/
  data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06, &
     -8.0059518998619764991d+09, 4.8544714258273622913d+12, &
     -1.3218168307321442305d+15/
!
!  Coefficients for 15.0 <= ABS(ARG)
!
  data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01, &
      -4.2843766903304806403d-01, 9.7356000150886612134d-02, &
      -3.2457723974465568321d-03,-3.6395264712121795296d-04, &
       1.6258661867440836395d-05,-3.6347578404608223492d-07/
  data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00, &
      -8.5017476463217924408d-01, 7.4212010813186530069d-02, &
      -2.2835624489492512649d-03, 3.7510433111922824643d-05/
  data pbar/3.98437500d-01/

  x = abs(arg)

  if (x < xsmall) then
!
!  Return for ABS(ARG) < XSMALL
!
    result = half * x

  else if (x < one5) then
!
!  XSMALL <= ABS(ARG) < 15.0
!
    xx = x * x

    sump = p(1)
    do j = 2, 15
      sump = sump * xx + p(j)
    end do

    xx = xx - two25
    sumq = ((((xx+q(1))*xx+q(2))*xx+q(3))*xx+q(4)) &
      * xx + q(5)
    result = (sump / sumq) * x
    if (jint == 2) result = result * exp(-x)

  else if ((jint == 1) .and. (x > xmax)) then

    result = xinf

  else
!
!  15.0 <= ABS(ARG)
!
    xx = one / x - rec15
    sump = ((((((pp(1)*xx+pp(2))*xx+pp(3))*xx+ &
      pp(4))*xx+pp(5))*xx+pp(6))*xx+pp(7))*xx+pp(8)
    sumq = (((((xx+qq(1))*xx+qq(2))*xx+qq(3))*xx+ &
      qq(4))*xx+qq(5))*xx+qq(6)
    result = sump / sumq

    if (jint /= 1) then

      result = (result + pbar) / sqrt(x)

    else
!
!  Calculation reformulated to avoid premature overflow
!
      if (x > xmax-one5) then
        a = exp(x-forty)
        b = exp40
         else
        a = exp(x)
        b = one
      end if

      result = ((result * a + pbar * a) / sqrt(x)) * b
!
!  Error return for ABS(ARG) > XMAX
!
    end if

  end if

  if (arg < zero) then
    result = -result
  end if

  return
end
subroutine calck0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCK0 evaluates modified Bessel K functions of order 0.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.  It contains two function type subprograms, BESK0
!    and BESEK0, and one subroutine type subprogram, CALCK0.
!    the calling statements for the primary entries are
!
!       Y=BESK0(X)
!   and
!       Y=BESEK0(X)
!
!   where the entry points correspond to the functions K0(X) and
!   EXP(X)*K0(X), respectively.  The routine CALCK0 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCK0 with the statement
!      CALL CALCK0(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!  function         Parameters for CALCK0
!   Call      ARG      RESULT      JINT
!
!     BESK0(ARG)   0 < ARG <= XMAX   K0(ARG)       1
!     BESEK0(ARG)     0 < ARG       EXP(ARG)*K0(ARG)     2
!
!   The main computation evaluates slightly modified forms of near
!   minimax rational approximations generated by Russon and Blair,
!   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!   1969.  This transportable program is patterned after the
!   machine-dependent FUNPACK packet NATSK0, but cannot match that
!   version for efficiency or accuracy.  This version uses rational
!   functions that theoretically approximate K-SUB-0(X) to at
!   least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XSMALL = Argument below which BESK0 and BESEK0 may
!    each be represented by a constant and a log.
!    largest X such that  1.0 + X = 1.0  to machine
!    precision.
!   XINF   = Largest positive machine number; approximately
!    beta**maxexp
!   XMAX   = Largest argument acceptable to BESK0;  Solution to
!    equation:
!       W(X) * (1-1/8X+9/128X**2) = beta**minexp
!    where  W(X) = EXP(-X)*SQRT(PI/2X)
!
! Error returns
!
!  The program returns the value XINF for ARG <= 0.0, and the
!  BESK0 entry returns the value 0.0 for ARG > XMAX.
!
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
      arg,f,g,one,p,pp,q,qq,result,sumf,sumg,sump,sumq,temp, &
      x,xinf,xmax,xsmall,xx,zero
  dimension p(6),q(2),pp(10),qq(10),f(4),g(3)
!
!  Mathematical constants
!
  data one/1.0d0/,zero/0.0d0/
!
!  Machine-dependent constants
!
  data xsmall/1.11d-16/,xinf/1.79d+308/,xmax/705.342d0/
!
!  Coefficients for XSMALL <=  ARG  <= 1.0
!
  data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
        1.1999463724910714109d+01, 4.6850901201934832188d+02, &
        5.9169059852270512312d+03, 2.4708152720399552679d+03/
  data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
  data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
       -1.7733784684952985886d+04,-4.0320340761145482298d+05/
  data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
       -1.6128136304458193998d+06/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
        3.1075408980684392399d+04, 1.0577068948034021957d+05, &
        1.7398867902565686251d+05, 1.5097646353289914539d+05, &
        7.1557062783764037541d+04, 1.8321525870183537725d+04, &
        2.3444738764199315021d+03, 1.1600249425076035558d+02/
  data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
        3.1474655750295278825d+04, 9.7418829762268075784d+04, &
        1.5144644673520157801d+05, 1.2689839587977598727d+05, &
        5.8824616785857027752d+04, 1.4847228371802360957d+04, &
        1.8821890840982713696d+03, 9.2556599177304839811d+01/

  x = arg

  if (x > zero) then

    if (x <= one) then
!
!  0.0 <  ARG  <= 1.0
!
      temp = log(x)

      if (x < xsmall) then
!
!  Return for small ARG
!
        result = p(6)/q(2) - temp

      else

        xx = x * x
        sump = ((((p(1)*xx + p(2))*xx + p(3))*xx + &
          p(4))*xx + p(5))*xx + p(6)
        sumq = (xx + q(1))*xx + q(2)
        sumf = ((f(1)*xx + f(2))*xx + f(3))*xx + f(4)
        sumg = ((xx + g(1))*xx + g(2))*xx + g(3)
        result = sump/sumq - xx*sumf*temp/sumg - temp
        if (jint == 2) result = result * exp(x)

      end if

    else if ((jint == 1) .and. (x > xmax)) then
!
!  Error return for ARG > XMAX
!
      result = zero

    else
!
!  1.0 < ARG
!
      xx = one / x
      sump = pp(1)
      do i = 2, 10
        sump = sump*xx + pp(i)
      end do
      sumq = xx
      do i = 1, 9
        sumq = (sumq + qq(i))*xx
      end do
      sumq = sumq + qq(10)
      result = sump / sumq / sqrt(x)
      if (jint == 1) result = result * exp(-x)

    end if

  else
!
!  Error return for ARG <= 0.0
!
    result = xinf

  end if

  return
end
subroutine calck1 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCK1 evaluates modifies Bessel K functions of order 1.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one,  K1(X)  and  EXP(X)*K1(X), for real arguments X.
!    It contains two function type subprograms, BESK1  and  BESEK1,
!    and one subroutine type subprogram, CALCK1.  The calling
!    statements for the primary entries are
!
!       Y=BESK1(X)
!   and
!       Y=BESEK1(X)
!
!   where the entry points correspond to the functions K1(X) and
!   EXP(X)*K1(X), respectively.  The routine CALCK1 is intended
!   for internal packet use only, all computations within the
!   packet being concentrated in this routine.  The function
!   subprograms invoke CALCK1 with the statement
!      CALL CALCK1(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!    function          Parameters for CALCK1
!    Call     ARG      RESULT      JINT
!
!     BESK1(ARG)  XLEAST < ARG < XMAX    K1(ARG)      1
!     BESEK1(ARG)     XLEAST < ARG   EXP(ARG)*K1(ARG)    2
!
!   The main computation evaluates slightly modified forms of near
!   minimax rational approximations generated by Russon and Blair,
!   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!   1969.  This transportable program is patterned after the
!   machine-dependent FUNPACK packet NATSK1, but cannot match that
!   version for efficiency or accuracy.  This version uses rational
!   functions that theoretically approximate K-SUB-1(X) to at
!   least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XLEAST = Smallest acceptable argument, i.e., smallest machine
!    number X such that 1/X is machine representable.
!   XSMALL = Argument below which BESK1(X) and BESEK1(X) may
!    each be represented by 1/X.  A safe value is the
!    largest X such that  1.0 + X = 1.0  to machine
!    precision.
!   XINF   = Largest positive machine number; approximately
!    beta**maxexp
!   XMAX   = Largest argument acceptable to BESK1;  Solution to
!    equation:
!       W(X) * (1+3/8X-15/128X**2) = beta**minexp
!    where  W(X) = EXP(-X)*SQRT(PI/2X)
!
! Error returns
!
!  The program returns the value XINF for ARG <= 0.0 and the
!   BESK1 entry returns the value 0.0 for ARG > XMAX.
!
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
      arg,f,g,one,p,pp,q,qq,result,sumf,sumg, &
      sump,sumq,x,xinf,xmax,xleast,xsmall,xx,zero
  dimension p(5),q(3),pp(11),qq(9),f(5),g(3)
!
!  Mathematical constants
!
  data one/1.0d0/,zero/0.0d0/
!
!  Machine-dependent constants
!
  data xleast/2.23d-308/,xsmall/1.11d-16/,xinf/1.79d+308/, &
       xmax/705.343d+0/
!
!  Coefficients for  XLEAST <=  ARG  <= 1.0
!
  data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1, &
        7.1885382604084798576d+3, 1.7733324035147015630d+5, &
        7.1938920065420586101d+5/
  data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, &
       -2.2149374878243304548d+6/
  data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1, &
       -4.5051623763436087023d+3,-1.4758069205414222471d+5, &
       -1.3531161492785421328d+6/
  data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4, &
       -2.7062322985570842656d+6/
!
!  Coefficients for  1.0 <  ARG
!
  data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0, &
        1.3182609918569941308d+2, 8.1094256146537402173d+2, &
        2.3123742209168871550d+3, 3.4540675585544584407d+3, &
        2.8590657697910288226d+3, 1.3319486433183221990d+3, &
        3.4122953486801312910d+2, 4.4137176114230414036d+1, &
        2.2196792496874548962d+0/
  data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2, &
        1.2082692316002348638d+3, 2.1181000487171943810d+3, &
        1.9448440788918006154d+3, 9.6929165726802648634d+2, &
        2.5951223655579051357d+2, 3.4552228452758912848d+1, &
        1.7710478032601086579d+0/

  x = arg

  if (x < xleast) then
!
!  Error return for  ARG  < XLEAST
!
    result = xinf

  else if (x <= one) then
!
!  XLEAST <=  ARG  <= 1.0
!
    if (x < xsmall) then
!
!  Return for small ARG
!
      result = one / x

    else

      xx = x * x
      sump = ((((p(1)*xx + p(2))*xx + p(3))*xx + p(4))*xx + p(5))*xx + q(3)
      sumq = ((xx + q(1))*xx + q(2))*xx + q(3)
      sumf = (((f(1)*xx + f(2))*xx + f(3))*xx + f(4))*xx + f(5)
      sumg = ((xx + g(1))*xx + g(2))*xx + g(3)
      result = (xx * log(x) * sumf/sumg + sump/sumq) / x
      if (jint == 2) result = result * exp(x)

    end if

  else if ((jint == 1) .and. (x > xmax)) then
!
!  Error return for  ARG  > XMAX
!
    result = zero

  else
!
!  1.0 <  ARG
!
    xx = one / x

    sump = pp(1)
    do i = 2, 11
      sump = sump * xx + pp(i)
    end do

    sumq = xx
    do i = 1, 8
      sumq = (sumq + qq(i)) * xx
    end do
    sumq = sumq + qq(9)

    result = sump / sumq / sqrt(x)
    if (jint == 1) result = result * exp(-x)

  end if

  return
end
subroutine caljy0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALJY0 evaluates Bessel J and Y functions of order 0.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.  It contains two function-type
!    subprograms,  BESJ0  and  BESY0,  and one subroutine-type
!    subprogram,  CALJY0.  The calling statements for the primary
!    entries are:
!
!       Y = BESJ0(X)
!   and
!       Y = BESY0(X),
!
!   where the entry points correspond to the functions J0(X) and Y0(X),
!   respectively.  The routine  CALJY0  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY0  with
!   the statement
!       CALL CALJY0(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!function      Parameters for CALJY0
!   call      ARG     RESULT      JINT
!
!     BESJ0(ARG)     |ARG| <= XMAX   J0(ARG)      0
!     BESY0(ARG)   0 < ARG <= XMAX    Y0(ARG)      1
!
!   The main computation uses unpublished minimax rational
!   approximations for X <= 8.0, and an approximation from the
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ0(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!    and COS must perform properly for  ABS(X) <= XMAX.
!    We recommend that XMAX be a small integer multiple of
!    sqrt(1/eps), where eps is the smallest positive number
!    such that  1+eps > 1.
!   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
!    to machine precision for all  ABS(X) <= XSMALL.
!    We recommend that  XSMALL < sqrt(eps)/beta, where beta
!    is the floating-point radix (usually 2 or 16).
!
! Error Returns
!
!  The program returns the value zero for  X > XMAX, and returns
!    -XINF when BESLY0 is called with a negative or zero argument.
!
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
     arg,ax,cons,down,eight,five5,four,one,oneov8,pi2,pj0, &
     pj1,plg,prod,py0,py1,py2,p0,p1,p17,qj0,qj1,qlg,qy0,qy1, &
     qy2,q0,q1,resj,result,r0,r1,sixty4,three,twopi,twopi1, &
     twopi2,two56,up,w,wsq,xden,xinf,xmax,xnum,xsmall,xj0, &
     xj1,xj01,xj02,xj11,xj12,xy,xy0,xy01,xy02,xy1,xy11,xy12, &
     xy2,xy21,xy22,z,zero,zsq
  dimension pj0(7),pj1(8),plg(4),py0(6),py1(7),py2(8),p0(6),p1(6), &
        qj0(5),qj1(7),qlg(4),qy0(5),qy1(6),qy2(7),q0(5),q1(5)
!
!  Mathematical constants
!  CONS = ln(.5) + Euler's gamma
!
  data zero,one,three,four,eight/0.0d0,1.0d0,3.0d0,4.0d0,8.0d0/, &
       five5,sixty4,oneov8,p17/5.5d0,64.0d0,0.125d0,1.716d-1/, &
       two56,cons/256.0d0,-1.1593151565841244881d-1/, &
       pi2,twopi/6.3661977236758134308d-1,6.2831853071795864769d0/, &
       twopi1,twopi2/6.28125d0,1.9353071795864769253d-3/
!
!  Machine-dependent constants
!
  data xmax/2.68d+08/,xsmall/3.72d-09/,xinf/1.79d+308/
!
!  Zeroes of Bessel functions
!
  data xj0/2.4048255576957727686d+0/,xj1/5.5200781102863106496d+0/, &
       xy0/8.9357696627916752158d-1/,xy1/3.9576784193148578684d+0/, &
       xy2/7.0860510603017726976d+0/, &
       xj01/ 616.0d+0/, xj02/-1.4244423042272313784d-03/, &
       xj11/1413.0d+0/, xj12/ 5.4686028631064959660d-04/, &
       xy01/ 228.0d+0/, xy02/ 2.9519662791675215849d-03/, &
       xy11/1013.0d+0/, xy12/ 6.4716931485786837568d-04/, &
       xy21/1814.0d+0/, xy22/ 1.1356030177269762362d-04/
!
!  Coefficients for rational approximation to ln(x/a)
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
       -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
       -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
!
  data pj0/6.6302997904833794242d+06,-6.2140700423540120665d+08, &
       2.7282507878605942706d+10,-4.1298668500990866786d+11, &
      -1.2117036164593528341d-01, 1.0344222815443188943d+02, &
      -3.6629814655107086448d+04/
  data qj0/4.5612696224219938200d+05, 1.3985097372263433271d+08, &
       2.6328198300859648632d+10, 2.3883787996332290397d+12, &
       9.3614022392337710626d+02/
!
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
!
  data pj1/4.4176707025325087628d+03, 1.1725046279757103576d+04, &
       1.0341910641583726701d+04,-7.2879702464464618998d+03, &
      -1.2254078161378989535d+04,-1.8319397969392084011d+03, &
       4.8591703355916499363d+01, 7.4321196680624245801d+02/
  data qj1/3.3307310774649071172d+02,-2.9458766545509337327d+03, &
       1.8680990008359188352d+04,-8.4055062591169562211d+04, &
       2.4599102262586308984d+05,-3.5783478026152301072d+05, &
      -2.5258076240801555057d+01/
!
!  Coefficients for rational approximation of
!  (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!  XSMALL  <  |X|  <=  3.0
!
  data py0/1.0102532948020907590d+04,-2.1287548474401797963d+06, &
       2.0422274357376619816d+08,-8.3716255451260504098d+09, &
       1.0723538782003176831d+11,-1.8402381979244993524d+01/
  data qy0/6.6475986689240190091d+02, 2.3889393209447253406d+05, &
       5.5662956624278251596d+07, 8.1617187777290363573d+09, &
       5.8873865738997033405d+11/
!
!  Coefficients for rational approximation of
!  (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!  3.0  <  |X|  <=  5.5
!
  data py1/-1.4566865832663635920d+04, 4.6905288611678631510d+06, &
       -6.9590439394619619534d+08, 4.3600098638603061642d+10, &
       -5.5107435206722644429d+11,-2.2213976967566192242d+13, &
        1.7427031242901594547d+01/
  data qy1/ 8.3030857612070288823d+02, 4.0669982352539552018d+05, &
        1.3960202770986831075d+08, 3.4015103849971240096d+10, &
        5.4266824419412347550d+12, 4.3386146580707264428d+14/
!
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!    5.5  <  |X|  <=  8.0
!
  data py2/ 2.1363534169313901632d+04,-1.0085539923498211426d+07, &
        2.1958827170518100757d+09,-1.9363051266772083678d+11, &
       -1.2829912364088687306d+11, 6.7016641869173237784d+14, &
       -8.0728726905150210443d+15,-1.7439661319197499338d+01/
  data qy2/ 8.7903362168128450017d+02, 5.3924739209768057030d+05, &
        2.4727219475672302327d+08, 8.6926121104209825246d+10, &
        2.2598377924042897629d+13, 3.9272425569640309819d+15, &
        3.4563724628846457519d+17/
!
!  Coefficients for Hart,s approximation,  |X| > 8.0
!
  data p0/3.4806486443249270347d+03, 2.1170523380864944322d+04, &
      4.1345386639580765797d+04, 2.2779090197304684302d+04, &
      8.8961548424210455236d-01, 1.5376201909008354296d+02/
  data q0/3.5028735138235608207d+03, 2.1215350561880115730d+04, &
      4.1370412495510416640d+04, 2.2779090197304684318d+04, &
      1.5711159858080893649d+02/
  data p1/-2.2300261666214198472d+01,-1.1183429920482737611d+02, &
      -1.8591953644342993800d+02,-8.9226600200800094098d+01, &
      -8.8033303048680751817d-03,-1.2441026745835638459d+00/
  data q1/1.4887231232283756582d+03, 7.2642780169211018836d+03, &
      1.1951131543434613647d+04, 5.7105024128512061905d+03, &
      9.0593769594993125859d+01/
!
!  Check for error conditions
!
  ax = abs(arg)

  if ((jint == 1) .and. (arg <= zero)) then
    result = -xinf
    return
  else if (ax > xmax) then
    result = zero
    return
  end if
!
!  Calculate J0 or Y0 for |ARG|  >  8.0
!
  if (ax > eight) then

    z = eight / ax
    w = ax / twopi
    w = aint(w) + oneov8
    w = (ax - w * twopi1) - w * twopi2
    zsq = z * z
    xnum = p0(5) * zsq + p0(6)
    xden = zsq + q0(5)
    up = p1(5) * zsq + p1(6)
    down = zsq + q1(5)

    do i = 1, 4
      xnum = xnum * zsq + p0(i)
      xden = xden * zsq + q0(i)
      up = up * zsq + p1(i)
      down = down * zsq + q1(i)
    end do

    r0 = xnum / xden
    r1 = up / down

    if (jint == 0) then
      result = sqrt(pi2/ax) * (r0*cos(w) - z*r1*sin(w))
    else
      result = sqrt(pi2/ax) * (r0*sin(w) + z*r1*cos(w))
    end if

    return

  end if

  if (ax <= xsmall) then
    if (jint == 0) then
      result = one
    else
      result = pi2 * (log(ax) + cons)
    end if
    return
  end if
!
!  Calculate J0 for appropriate interval, preserving
!  accuracy near the zero of J0
!
  zsq = ax * ax

  if (ax <= four) then

    xnum = (pj0(5) * zsq + pj0(6)) * zsq + pj0(7)
    xden = zsq + qj0(5)
    do i = 1, 4
      xnum = xnum * zsq + pj0(i)
      xden = xden * zsq + qj0(i)
    end do
    prod = ((ax - xj01/two56) - xj02) * (ax + xj0)

  else

    wsq = one - zsq / sixty4
    xnum = pj1(7) * wsq + pj1(8)
    xden = wsq + qj1(7)
    do i = 1, 6
      xnum = xnum * wsq + pj1(i)
      xden = xden * wsq + qj1(i)
    end do
    prod = (ax + xj1) * ((ax - xj11/two56) - xj12)

  end if

  result = prod * xnum / xden

  if (jint == 0) then
    return
  end if
!
!  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!  where xn is a zero of Y0
!
  if (ax <= three) then

    up = (ax-xy01/two56)-xy02
    xy = xy0

  else if (ax <= five5) then

    up = (ax-xy11/two56)-xy12
    xy = xy1

  else
    up = (ax-xy21/two56)-xy22
    xy = xy2
  end if

  down = ax + xy

  if (abs(up) < p17*down) then
    w = up/down
    wsq = w*w
    xnum = plg(1)
    xden = wsq + qlg(1)
    do i = 2, 4
      xnum = xnum*wsq + plg(i)
      xden = xden*wsq + qlg(i)
    end do
    resj = pi2 * result * w * xnum/xden
  else
    resj = pi2 * result * log(ax/xy)
  end if
!
!  Now calculate Y0 for appropriate interval, preserving
!  accuracy near the zero of Y0
!
  if (ax <= three) then

    xnum = py0(6) * zsq + py0(1)
    xden = zsq + qy0(1)

    do i = 2, 5
      xnum = xnum * zsq + py0(i)
      xden = xden * zsq + qy0(i)
    end do

  else if (ax <= five5) then

    xnum = py1(7) * zsq + py1(1)
    xden = zsq + qy1(1)

    do i = 2, 6
      xnum = xnum * zsq + py1(i)
      xden = xden * zsq + qy1(i)
    end do

  else

    xnum = py2(8) * zsq + py2(1)
    xden = zsq + qy2(1)

    do i = 2, 7
      xnum = xnum * zsq + py2(i)
      xden = xden * zsq + qy2(i)
    end do

  end if

  result = resj + up * down * xnum / xden

  return
end
subroutine caljy1 ( arg, result, jint )

!*****************************************************************************80
!
!! CALJY1 evaluates Bessel J and Y functions of order 1.
!
!  Discussion:
!
!    This routine computes first-order Bessel functions of the first and
!    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!    for Y1, and |X| <= XMAX for J1.  It contains two function-type
!    subprograms,  BESJ1  and  BESY1,  and one subroutine-type
!    subprogram,  CALJY1.  The calling statements for the primary
!    entries are:
!
!       Y = BESJ1(X)
!   and
!       Y = BESY1(X),
!
!   where the entry points correspond to the functions J1(X) and Y1(X),
!   respectively.  The routine  CALJY1  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY1  with
!   the statement
!       CALL CALJY1(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!function      Parameters for CALJY1
!   call      ARG     RESULT      JINT
!
!     BESJ1(ARG)     |ARG| <= XMAX   J1(ARG)      0
!     BESY1(ARG)   0 < ARG <= XMAX    Y1(ARG)      1
!
!   The main computation uses unpublished minimax rational
!   approximations for X <= 8.0, and an approximation from the
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ1(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!    and COS must perform properly for  ABS(X) <= XMAX.
!    We recommend that XMAX be a small integer multiple of
!    sqrt(1/eps), where eps is the smallest positive number
!    such that  1+eps > 1.
!   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
!    to machine precision for all  ABS(X) <= XSMALL.
!    We recommend that  XSMALL < sqrt(eps)/beta, where beta
!    is the floating-point radix (usually 2 or 16).
!
! Error Returns
!
!  The program returns the value zero for  X > XMAX, and returns
!    -XINF when BESLY1 is called with a negative or zero argument.
!
  implicit none

  integer ( kind = 4 ) i,jint
  dimension pj0(7),pj1(8),plg(4),py0(7),py1(9),p0(6),p1(6), &
        qj0(5),qj1(7),qlg(4),qy0(6),qy1(8),q0(6),q1(6)
  real ( kind = 8 ) &
     arg,ax,down,eight,four,half,pi2,pj0,pj1,plg,prod,py0, &
     py1,p0,p1,p17,qj0,qj1,qlg,qy0,qy1,q0,q1,resj,result, &
     rtpi2,r0,r1,throv8,twopi,twopi1,twopi2,two56,up,w,wsq, &
     xden,xinf,xmax,xnum,xsmall,xj0,xj1,xj01,xj02,xj11,xj12, &
     xy,xy0,xy01,xy02,xy1,xy11,xy12,z,zero,zsq
!
!  Mathematical constants
!
  data eight/8.0d0/, &
       four/4.0d0/,half/0.5d0/,throv8/0.375d0/, &
       pi2/6.3661977236758134308d-1/,p17/1.716d-1/ &
       twopi/6.2831853071795864769d+0/,zero/0.0d0/, &
       twopi1/6.28125d0/,twopi2/1.9353071795864769253d-03/ &
       two56/256.0d+0/,rtpi2/7.9788456080286535588d-1/
!
!  Machine-dependent constants
!
  data xmax/2.68d+08/,xsmall/3.72d-09/,xinf/1.79d+308/
!
!  Zeroes of Bessel functions
!
  data xj0/3.8317059702075123156d+0/,xj1/7.0155866698156187535d+0/, &
       xy0/2.1971413260310170351d+0/,xy1/5.4296810407941351328d+0/, &
       xj01/ 981.0d+0/, xj02/-3.2527979248768438556d-04/, &
       xj11/1796.0d+0/, xj12/-3.8330184381246462950d-05/, &
       xy01/ 562.0d+0/, xy02/ 1.8288260310170351490d-03/, &
       xy11/1390.0d+0/, xy12/-6.4592058648672279948d-06/
!
!  Coefficients for rational approximation to ln(x/a)
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
       -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
       -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
!
  data pj0/9.8062904098958257677d+05,-1.1548696764841276794d+08, &
     6.6781041261492395835d+09,-1.4258509801366645672d+11, &
    -4.4615792982775076130d+03, 1.0650724020080236441d+01, &
    -1.0767857011487300348d-02/
  data qj0/5.9117614494174794095d+05, 2.0228375140097033958d+08, &
     4.2091902282580133541d+10, 4.1868604460820175290d+12, &
     1.0742272239517380498d+03/
!
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
!
  data pj1/4.6179191852758252280d+00,-7.1329006872560947377d+03, &
     4.5039658105749078904d+06,-1.4437717718363239107d+09, &
     2.3569285397217157313d+11,-1.6324168293282543629d+13, &
     1.1357022719979468624d+14, 1.0051899717115285432d+15/
  data qj1/1.1267125065029138050d+06, 6.4872502899596389593d+08, &
     2.7622777286244082666d+11, 8.4899346165481429307d+13, &
     1.7128800897135812012d+16, 1.7253905888447681194d+18, &
     1.3886978985861357615d+03/
!
!  Coefficients for rational approximation of
!  (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!  XSMALL  <  |X|  <=  4.0
!
  data py0/2.2157953222280260820d+05,-5.9157479997408395984d+07, &
       7.2144548214502560419d+09,-3.7595974497819597599d+11, &
       5.4708611716525426053d+12, 4.0535726612579544093d+13, &
      -3.1714424660046133456d+02/
  data qy0/8.2079908168393867438d+02, 3.8136470753052572164d+05, &
       1.2250435122182963220d+08, 2.7800352738690585613d+10, &
       4.1272286200406461981d+12, 3.0737873921079286084d+14/
!
!  Coefficients for rational approximation of
!  (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!  4.0  <  |X|  <=  8.0
!
  data py1/ 1.9153806858264202986d+06,-1.1957961912070617006d+09, &
        3.7453673962438488783d+11,-5.9530713129741981618d+13, &
        4.0686275289804744814d+15,-2.3638408497043134724d+16, &
       -5.6808094574724204577d+18, 1.1514276357909013326d+19, &
       -1.2337180442012953128d+03/
  data qy1/ 1.2855164849321609336d+03, 1.0453748201934079734d+06, &
        6.3550318087088919566d+08, 3.0221766852960403645d+11, &
        1.1187010065856971027d+14, 3.0837179548112881950d+16, &
        5.6968198822857178911d+18, 5.3321844313316185697d+20/
!
!  Coefficients for Hart's approximation,  |X| > 8.0
!
  data p0/-1.0982405543459346727d+05,-1.5235293511811373833d+06, &
       -6.6033732483649391093d+06,-9.9422465050776411957d+06, &
       -4.4357578167941278571d+06,-1.6116166443246101165d+03/
  data q0/-1.0726385991103820119d+05,-1.5118095066341608816d+06, &
       -6.5853394797230870728d+06,-9.9341243899345856590d+06, &
       -4.4357578167941278568d+06,-1.4550094401904961825d+03/
  data p1/ 1.7063754290207680021d+03, 1.8494262873223866797d+04, &
        6.6178836581270835179d+04, 8.5145160675335701966d+04, &
        3.3220913409857223519d+04, 3.5265133846636032186d+01/
  data q1/ 3.7890229745772202641d+04, 4.0029443582266975117d+05, &
        1.4194606696037208929d+06, 1.8194580422439972989d+06, &
        7.0871281941028743574d+05, 8.6383677696049909675d+02/
!
!  Check for error conditions
!
  ax = abs(arg)

  if ((jint == 1) .and. ((arg <= zero) .or. &
     ((arg < half) .and. (ax*xinf < pi2)))) then
    result = -xinf
    return
  else if (ax > xmax) then
    result = zero
    return
  end if
!
!  Calculate J1 or Y1 for |ARG|  >  8.0
!
  if (ax > eight) then

    z = eight / ax
    w = aint(ax/twopi) + throv8
    w = (ax - w * twopi1) - w * twopi2
    zsq = z * z
    xnum = p0(6)
    xden = zsq + q0(6)
    up = p1(6)
    down = zsq + q1(6)
    do i = 1, 5
      xnum = xnum * zsq + p0(i)
      xden = xden * zsq + q0(i)
      up = up * zsq + p1(i)
      down = down * zsq + q1(i)
    end do
    r0 = xnum / xden
    r1 = up / down

    if (jint == 0) then
      result = (rtpi2/sqrt(ax)) * (r0*cos(w) - z*r1*sin(w))
    else
      result = (rtpi2/sqrt(ax)) * (r0*sin(w) + z*r1*cos(w))
    end if

    if ((jint == 0) .and. (arg < zero)) then
      result = -result
    end if
 
    return

  else if (ax <= xsmall) then

    if (jint == 0) then
      result = arg * half
    else
      result = -pi2 / ax
    end if
    return
  end if
!
!  Calculate J1 for appropriate interval, preserving
!  accuracy near the zero of J1
!
  zsq = ax * ax

  if (ax <= four) then
    xnum = (pj0(7) * zsq + pj0(6)) * zsq + pj0(5)
    xden = zsq + qj0(5)
    do i = 1, 4
      xnum = xnum * zsq + pj0(i)
      xden = xden * zsq + qj0(i)
    end do
    prod = arg * ((ax - xj01/two56) - xj02) * (ax + xj0)
  else
    xnum = pj1(1)
    xden = (zsq + qj1(7)) * zsq + qj1(1)
    do i = 2, 6
      xnum = xnum * zsq + pj1(i)
      xden = xden * zsq + qj1(i)
    end do
    xnum = xnum * (ax - eight) * (ax + eight) + pj1(7)
    xnum = xnum * (ax - four) * (ax + four) + pj1(8)
    prod = arg * ((ax - xj11/two56) - xj12) * (ax + xj1)
  end if

  result = prod * (xnum / xden)

  if (jint == 0) then
    return
  end if
!
!  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
!  where xn is a zero of Y1
!
  if (ax <= four) then
    up = (ax-xy01/two56)-xy02
    xy = xy0
  else
    up = (ax-xy11/two56)-xy12
    xy = xy1
  end if

  down = ax + xy

  if (abs(up) < p17*down) then
    w = up/down
    wsq = w*w
    xnum = plg(1)
    xden = wsq + qlg(1)
    do i = 2, 4
      xnum = xnum*wsq + plg(i)
      xden = xden*wsq + qlg(i)
    end do
    resj = pi2 * result * w * xnum/xden
  else
    resj = pi2 * result * log(ax/xy)
  end if
!
!  Now calculate Y1 for appropriate interval, preserving
!  accuracy near the zero of Y1
!
  if (ax <= four) then
    xnum = py0(7) * zsq + py0(1)
    xden = zsq + qy0(1)
    do i = 2, 6
      xnum = xnum * zsq + py0(i)
      xden = xden * zsq + qy0(i)
    end do
  else
    xnum = py1(9) * zsq + py1(1)
    xden = zsq + qy1(1)
    do i = 2, 8
      xnum = xnum * zsq + py1(i)
      xden = xden * zsq + qy1(i)
    end do
  end if

  result = resj + (up*down/ax) * xnum / xden

  return
end
subroutine calerf ( arg, result, jint )

!*****************************************************************************80
!
!! CALERF evaluates the error function and related quantities.
!
!  Discussion:
!
!    This routine evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!    for a real argument  x.  It contains three FUNCTION type
!    subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!    and one SUBROUTINE type subprogram, CALERF.  The calling
!    statements for the primary entries are:
!
!       Y=ERF(X)     (or   Y=DERF(X)),
!
!       Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!       Y=ERFCX(X)   (or   Y=DERFCX(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!      CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!  function         Parameters for CALERF
!   call      ARG      Result      JINT
!
!     ERF(ARG)  ANY REAL ARGUMENT     ERF(ARG)      0
!     ERFC(ARG)     ABS(ARG) < XBIG    ERFC(ARG)     1
!     ERFCX(ARG)    XNEG < ARG < XMAX   ERFCX(ARG)    2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   XMIN   = the smallest positive floating-point number.
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!    the negative of the solution to the equation
!    2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!    2*x/sqrt(pi)  and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!    the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!    W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!    machine precision.  A conservative value is
!    1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!    of XINF and 1/[sqrt(pi)*XMIN].
!
! Error returns
!
!  The program returns  ERFC = 0  for  ARG >= XBIG;
!
!           ERFCX = XINF  for  ARG < XNEG;
!  and
!           ERFCX = 0     for  ARG >= XMAX.
!
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
       a,arg,b,c,d,del,four,half,p,one,q,result,sixten,sqrpi, &
       two,thresh,x,xbig,xden,xhuge,xinf,xmax,xneg,xnum,xsmall, &
       y,ysq,zero
  dimension a(5),b(4),c(9),d(8),p(6),q(5)
!
!  Mathematical constants
!
  data four,one,half,two,zero/4.0d0,1.0d0,0.5d0,2.0d0,0.0d0/, &
       sqrpi/5.6418958354775628695d-1/,thresh/0.46875d0/, &
       sixten/16.0d0/
!
!  Machine-dependent constants
!
  data xinf,xneg,xsmall/1.79d308,-26.628d0,1.11d-16/, &
       xbig,xhuge,xmax/26.543d0,6.71d7,2.53d307/
!
!  Coefficients for approximation to  erf  in first interval
!
  data a/3.16112374387056560d00,1.13864154151050156d02, &
     3.77485237685302021d02,3.20937758913846947d03, &
     1.85777706184603153d-1/
  data b/2.36012909523441209d01,2.44024637934444173d02, &
     1.28261652607737228d03,2.84423683343917062d03/
!
!  Coefficients for approximation to  erfc  in second interval
!
  data c/5.64188496988670089d-1,8.88314979438837594d0, &
     6.61191906371416295d01,2.98635138197400131d02, &
     8.81952221241769090d02,1.71204761263407058d03, &
     2.05107837782607147d03,1.23033935479799725d03, &
     2.15311535474403846d-8/
  data d/1.57449261107098347d01,1.17693950891312499d02, &
     5.37181101862009858d02,1.62138957456669019d03, &
     3.29079923573345963d03,4.36261909014324716d03, &
     3.43936767414372164d03,1.23033935480374942d03/
!
!  Coefficients for approximation to  erfc  in third interval
!
  data p/3.05326634961232344d-1,3.60344899949804439d-1, &
     1.25781726111229246d-1,1.60837851487422766d-2, &
     6.58749161529837803d-4,1.63153871373020978d-2/
  data q/2.56852019228982242d00,1.87295284992346047d00, &
     5.27905102951428412d-1,6.05183413124413191d-2, &
     2.33520497626869185d-3/

  x = arg
  y = abs(x)

  if (y <= thresh) then
!
!  Evaluate  erf  for  |X| <= 0.46875
!
    ysq = zero
    if (y > xsmall) ysq = y * y
    xnum = a(5)*ysq
    xden = ysq
    do i = 1, 3
      xnum = (xnum + a(i)) * ysq
      xden = (xden + b(i)) * ysq
    end do
    result = x * (xnum + a(4)) / (xden + b(4))
    if (jint /= 0) result = one - result
    if (jint == 2) result = exp(ysq) * result

    return
!
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!
  else if (y <= four) then

    xnum = c(9)*y
    xden = y
    do i = 1, 7
      xnum = (xnum + c(i)) * y
      xden = (xden + d(i)) * y
    end do

    result = (xnum + c(8)) / (xden + d(8))

    if (jint /= 2) then
      ysq = aint(y*sixten)/sixten
      del = (y-ysq)*(y+ysq)
      result = exp(-ysq*ysq) * exp(-del) * result
    end if
!
!  Evaluate  erfc  for |X| > 4.0
!
  else

    result = zero

    if (y >= xbig) then

      if ((jint /= 2) .or. (y >= xmax)) then
        go to 300
      end if

      if (y >= xhuge) then
        result = sqrpi / y
        go to 300
      end if

    end if

    ysq = one / (y * y)
    xnum = p(6)*ysq
    xden = ysq
    do i = 1, 4
      xnum = (xnum + p(i)) * ysq
      xden = (xden + q(i)) * ysq
    end do

    result = ysq *(xnum + p(5)) / (xden + q(5))
    result = (sqrpi -  result) / y

    if (jint /= 2) then
      ysq = aint(y*sixten)/sixten
      del = (y-ysq)*(y+ysq)
      result = exp(-ysq*ysq) * exp(-del) * result
    end if

  end if
!
!  Fix up for negative argument, erf, etc.
!
300 continue

  if (jint == 0) then

    result = (half - result) + half
    if (x < zero) result = -result

  else if (jint == 1) then

    if (x < zero) result = two - result

  else

    if (x < zero) then
      if (x < xneg) then
        result = xinf
      else
        ysq = aint(x*sixten)/sixten
        del = (x-ysq)*(x+ysq)
        y = exp(ysq*ysq) * exp(del)
        result = (y+y) - result
      end if
    end if
  end if

  return
end
function daw ( xx )

!*****************************************************************************80
!
!! DAW evaluates Dawson's integral
!
!  Discussion:
!
!    This function evaluates Dawson's integral,
!      f(x) = exp ( -x^2 ) * integral ( 0 <= t <= x ) exp ( t^2 ) dt
!    for a real argument x.
!
!    The calling sequence for this function is
!      Y=DAW(X)
!    The main computation uses rational Chebyshev approximations
!    published in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek
!    and Thacher.  This transportable program is patterned after the
!    machine-dependent FUNPACK program DDAW(X), but cannot match that
!    version for efficiency or accuracy.  This version uses rational
!    approximations that are theoretically accurate to about 19
!    significant decimal digits.  The accuracy achieved depends on the
!    arithmetic system, the compiler, the intrinsic functions, and
!    proper selection of the machine-dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   XINF   = largest positive machine number
!   XMIN   = the smallest positive machine number.
!   EPS    = smallest positive number such that 1+eps > 1.
!    Approximately  beta**(-p), where beta is the machine
!    radix and p is the number of significant base-beta
!    digits in a floating-point number.
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XMAX   = absolute argument beyond which DAW(X) underflows.
!    XMAX = min(0.5/xmin, xinf).
!   XSMALL = absolute argument below DAW(X)  may be represented
!    by X.  We recommend XSMALL = sqrt(eps).
!   XLARGE = argument beyond which DAW(X) may be represented by
!    1/(2x).  We recommend XLARGE = 1/sqrt(eps).
!
! Error Returns
!
!  The program returns 0.0 for |X| > XMAX.
!
  implicit none

  real ( kind = 8 ) daw
  integer ( kind = 4 ) i
  real ( kind = 8 ) &
       frac,half,one,one225,p1,p2,p3,p4,q1,q2,q3,q4,six25, &
       sump,sumq,two5,w2,x,xx,y,xlarge,xmax,xsmall,zero
  dimension p1(10),p2(10),p3(10),p4(10),q1(10),q2(9),q3(9),q4(9)
!
!  Mathematical constants.
!
  data zero,half,one/0.0d0,0.5d0,1.0d0/, &
       six25,one225,two5/6.25d0,12.25d0,25.0d0/
!
!  Machine-dependent constants
!
  data xsmall/1.05d-08/, xlarge/9.49d+07/, xmax/2.24d+307/
!
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!
  data p1/-2.69020398788704782410d-12, 4.18572065374337710778d-10, &
      -1.34848304455939419963d-08, 9.28264872583444852976d-07, &
      -1.23877783329049120592d-05, 4.07205792429155826266d-04, &
      -2.84388121441008500446d-03, 4.70139022887204722217d-02, &
      -1.38868086253931995101d-01, 1.00000000000000000004d+00/
  data q1/ 1.71257170854690554214d-10, 1.19266846372297253797d-08, &
       4.32287827678631772231d-07, 1.03867633767414421898d-05, &
       1.78910965284246249340d-04, 2.26061077235076703171d-03, &
       2.07422774641447644725d-02, 1.32212955897210128811d-01, &
       5.27798580412734677256d-01, 1.00000000000000000000d+00/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [2.5, 3.5)
!
  data p2/-1.70953804700855494930d+00,-3.79258977271042880786d+01, &
       2.61935631268825992835d+01, 1.25808703738951251885d+01, &
      -2.27571829525075891337d+01, 4.56604250725163310122d+00, &
      -7.33080089896402870750d+00, 4.65842087940015295573d+01, &
      -1.73717177843672791149d+01, 5.00260183622027967838d-01/
  data q2/ 1.82180093313514478378d+00, 1.10067081034515532891d+03, &
      -7.08465686676573000364d+00, 4.53642111102577727153d+02, &
       4.06209742218935689922d+01, 3.02890110610122663923d+02, &
       1.70641269745236227356d+02, 9.51190923960381458747d+02, &
       2.06522691539642105009d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [3.5, 5.0]
!
  data p3/-4.55169503255094815112d+00,-1.86647123338493852582d+01, &
      -7.36315669126830526754d+00,-6.68407240337696756838d+01, &
       4.84507265081491452130d+01, 2.69790586735467649969d+01, &
      -3.35044149820592449072d+01, 7.50964459838919612289d+00, &
      -1.48432341823343965307d+00, 4.99999810924858824981d-01/
  data q3/ 4.47820908025971749852d+01, 9.98607198039452081913d+01, &
       1.40238373126149385228d+01, 3.48817758822286353588d+03, &
      -9.18871385293215873406d+00, 1.24018500009917163023d+03, &
      -6.88024952504512254535d+01,-2.31251575385145143070d+00, &
       2.50041492369922381761d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  |x| > 5.0
!
  data p4/-8.11753647558432685797d+00,-3.84043882477454453430d+01, &
      -2.23787669028751886675d+01,-2.88301992467056105854d+01, &
      -5.99085540418222002197d+00,-1.13867365736066102577d+01, &
      -6.52828727526980741590d+00,-4.50002293000355585708d+00, &
      -2.50000000088955834952d+00, 5.00000000000000488400d-01/
  data q4/ 2.69382300417238816428d+02, 5.04198958742465752861d+01, &
       6.11539671480115846173d+01, 2.08210246935564547889d+02, &
       1.97325365692316183531d+01,-1.22097010558934838708d+01, &
      -6.99732735041547247161d+00,-2.49999970104184464568d+00, &
       7.49999999999027092188d-01/

  x = xx

  if (abs(x) > xlarge) then

    if (abs(x) <= xmax) then
      daw = half / x
    else
      daw = zero
    end if

  else if (abs(x) < xsmall) then

    daw = x

  else

    y = x * x

    if (y < six25) then
!
!  ABS(X) < 2.5
!
      sump = p1(1)
      sumq = q1(1)
      do i = 2, 10
        sump = sump * y + p1(i)
        sumq = sumq * y + q1(i)
      end do
      daw = x * sump / sumq

    else if (y < one225) then
!
!  2.5 <= ABS(X) < 3.5
!
      frac = zero
      do i = 1, 9
        frac = q2(i) / (p2(i) + y + frac)
      end do
      daw = (p2(10) + frac) / x

    else if (y < two5) then
!
!  3.5 <= ABS(X) < 5.0
!
      frac = zero
      do i = 1, 9
        frac = q3(i) / (p3(i) + y + frac)
      end do
      daw = (p3(10) + frac) / x

    else
!
!  5.0 <= ABS(X) <= XLARGE
!
      w2 = one / x / x
      frac = zero
      do i = 1, 9
        frac = q4(i) / (p4(i) + y + frac)
      end do
      frac = p4(10) + frac
      daw = (half + half * w2 * frac) / x

    end if

  end if

  return
end
function derf ( x )

!*****************************************************************************80
!
!! DERF evaluates the error function.
!
!  Discussion:
!
!    This function computes approximate values for erf(x).
!    See comments heading CALERF.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) derf
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 0
  call calerf(x,result,jint)
  derf = result

  return
end
function derfc ( x )

!*****************************************************************************80
!
!! DERFC evaluates the complementary error function.
!
!  Discussion:
!
!    This function computes approximate values for erfc(x).
!    See comments heading CALERF.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) derfc
  integer ( kind = 4 ) jint
  real ( kind = 8 ) x, result

  jint = 1
  call calerf(x,result,jint)
  derfc = result

  return
end
function derfcx ( x )

!*****************************************************************************80
!
!! DERFCX evaluates exp(x^2) * erfc(x).
!
!  Discussion:
!
!    This function computes approximate values for exp(x*x) * erfc(x).
!    See comments heading CALERF.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) derfcx
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 2
  call calerf ( x, result, jint )
  derfcx = result

  return
end
function dgamma ( x )

!*****************************************************************************80
!
!! DGAMMA evaluates the gamma function.
!
!  Discussion:
!
!    This function calculates the GAMMA function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for X >= 12 are from reference 2.
!    The accuracy achieved depends on the arithmetic system, the
!    compiler, the intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
! XBIG   - the largest argument for which GAMMA(X) is representable
!      in the machine, i.e., the solution to the equation
!      GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!      approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!      1.0+EPS > 1.0
! XMININ - the smallest positive floating-point number such that
!      1/XMININ is machine representable
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!  Reference:
!
!    An Overview of Software Development for Special
!    functions", W. J. Cody, Lecture Notes in Mathematics,
!    506, Numerical Analysis Dundee, 1975, G. A. Watson
!    (ed.), Springer Verlag, Berlin, 1976.
!
!    Computer Approximations, Hart, Et. Al., Wiley and
!    sons, New York, 1968.
!
  implicit none

  real ( kind = 8 ) dgamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  logical ( kind = 4 ) parity
  real ( kind = 8 ) &
      c,conv,eps,fact,half,one,p,pi,q,res,sqrtpi,sum,twelve, &
      two,x,xbig,xden,xinf,xminin,xnum,y,y1,ysq,z,zero
  dimension c(7),p(8),q(8)
!
!  Mathematical constants
!
  data one,half,twelve,two,zero/1.0d0,0.5d0,12.0d0,2.0d0,0.0d0/, &
       sqrtpi/0.9189385332046727417803297d0/, &
       pi/3.1415926535897932384626434d0/
!
!  Machine dependent parameters
!
  data xbig,xminin,eps/171.624d0,2.23d-308,2.22d-16/, &
       xinf/1.79d308/
!
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!
  data p/-1.71618513886549492533811d+0,2.47656508055759199108314d+1, &
     -3.79804256470945635097577d+2,6.29331155312818442661052d+2, &
     8.66966202790413211295064d+2,-3.14512729688483675254357d+4, &
     -3.61444134186911729807069d+4,6.64561438202405440627855d+4/
  data q/-3.08402300119738975254353d+1,3.15350626979604161529144d+2, &
    -1.01515636749021914166146d+3,-3.10777167157231109440444d+3, &
      2.25381184209801510330112d+4,4.75584627752788110767815d+3, &
    -1.34659959864969306392456d+5,-1.15132259675553483497211d+5/
!
!  Coefficients for minimax approximation over (12, INF).
!
  data c/-1.910444077728d-03,8.4171387781295d-04, &
       -5.952379913043012d-04,7.93650793500350248d-04, &
       -2.777777777777681622553d-03,8.333333333333333331554247d-02, &
    5.7083835261d-03/
!
!  Statement functions for conversion between integer and float
!
  conv(i) = dble(i)
  parity = .false.
  fact = one
  n = 0
  y = x

  if (y <= zero) then
!
!  Argument is negative
!
    y = -x
    y1 = aint(y)
    res = y - y1

    if (res /= zero) then
      if (y1 /= aint(y1*half)*two) parity = .true.
      fact = -pi / sin(pi*res)
      y = y + one
    else
      res = xinf
      dgamma = res
      return
    end if

  end if
!
!  Argument is positive
!
  if (y < eps) then
!
!  Argument < EPS
!
    if (y >= xminin) then
      res = one / y
    else
      res = xinf
      dgamma = res
      return
    end if

  else if (y < twelve) then

    y1 = y
    if (y < one) then
!
!  0.0 < argument < 1.0
!
      z = y
      y = y + one

    else
!
!  1.0 < argument < 12.0, reduce argument if necessary
!
      n = int(y) - 1
      y = y - conv(n)
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = (xnum + p(i)) * z
      xden = xden * z + q(i)
    end do
    res = xnum / xden + one

    if (y1 < y) then
!
!  Adjust result for case  0.0 < argument < 1.0
!
      res = res / y1

    else if (y1 > y) then
!
!  Adjust result for case  2.0 < argument < 12.0
!
      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for argument >= 12.0,
!
    if (y <= xbig) then
      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum/y - y + sqrtpi
      sum = sum + (y-half)*log(y)
      res = exp(sum)

    else

      res = xinf
      dgamma = res
      return

    end if

  end if
!
!  Final adjustments and return
!
  if (parity) res = -res
  if (fact /= one) res = fact / res

  dgamma = res

  return
end
function dlgama ( x )

!*****************************************************************************80
!
!! DLGAMA calculates the logarithm of the gamma function.
!
!  Discussion:
!
!    This function calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.  The accuracy achieved depends
!    on the arithmetic system, the compiler, the intrinsic functions,
!    and proper selection of the machine-dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - largest argument for which LN(GAMMA(X)) is representable
!      in the machine, i.e., the solution to the equation
!      LN(GAMMA(XBIG)) = beta**maxexp
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
! XINF   - largest machine representable floating-point number;
!      approximately beta**maxexp.
! EPS    - The smallest positive floating-point number such that
!      1.0+EPS > 1.0
! FRTBIG - Rough estimate of the fourth root of XBIG
!
! Error returns
!
!  The program returns the value XINF for X <= 0.0 or when
!     overflow would occur.  The computation is believed to
!     be free of underflow and overflow.
!
!  Reference:
!
!    W. J. Cody and K. E. Hillstrom, 
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp. 
!    21, 1967, pp. 198-203.
!
!  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
!     1969.
!
!  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
!     York, 1968.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) dlgama
  real ( kind = 8 ) &
      c,corr,d1,d2,d4,eps,frtbig,four,half,one,pnt68,p1,p2,p4, &
      q1,q2,q4,res,sqrtpi,thrhal,twelve,two,x,xbig,xden,xinf, &
      xm1,xm2,xm4,xnum,y,ysq,zero
  dimension c(7),p1(8),p2(8),p4(8),q1(8),q2(8),q4(8)
!
!  Mathematical constants
!
  data one,half,twelve,zero/1.0d0,0.5d0,12.0d0,0.0d0/, &
       four,thrhal,two,pnt68/4.0d0,1.5d0,2.0d0,0.6796875d0/, &
       sqrtpi/0.9189385332046727417803297d0/
!
!  Machine dependent parameters
!
  data xbig,xinf,eps,frtbig/2.55d305,1.79d308,2.22d-16,2.25d76/
!
!  Numerator and denominator coefficients for rational minimax
!     approximation over (0.5,1.5).
!
  data d1/-5.772156649015328605195174d-1/
  data p1/4.945235359296727046734888d0,2.018112620856775083915565d2, &
      2.290838373831346393026739d3,1.131967205903380828685045d4, &
      2.855724635671635335736389d4,3.848496228443793359990269d4, &
      2.637748787624195437963534d4,7.225813979700288197698961d3/
  data q1/6.748212550303777196073036d1,1.113332393857199323513008d3, &
      7.738757056935398733233834d3,2.763987074403340708898585d4, &
      5.499310206226157329794414d4,6.161122180066002127833352d4, &
      3.635127591501940507276287d4,8.785536302431013170870835d3/
!
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (1.5,4.0).
!
  data d2/4.227843350984671393993777d-1/
  data p2/4.974607845568932035012064d0,5.424138599891070494101986d2, &
      1.550693864978364947665077d4,1.847932904445632425417223d5, &
      1.088204769468828767498470d6,3.338152967987029735917223d6, &
      5.106661678927352456275255d6,3.074109054850539556250927d6/
  data q2/1.830328399370592604055942d2,7.765049321445005871323047d3, &
      1.331903827966074194402448d5,1.136705821321969608938755d6, &
      5.267964117437946917577538d6,1.346701454311101692290052d7, &
      1.782736530353274213975932d7,9.533095591844353613395747d6/
!
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (4.0,12.0).
!
  data d4/1.791759469228055000094023d0/
  data p4/1.474502166059939948905062d4,2.426813369486704502836312d6, &
      1.214755574045093227939592d8,2.663432449630976949898078d9, &
    2.940378956634553899906876d10,1.702665737765398868392998d11, &
    4.926125793377430887588120d11,5.606251856223951465078242d11/
  data q4/2.690530175870899333379843d3,6.393885654300092398984238d5, &
      4.135599930241388052042842d7,1.120872109616147941376570d9, &
    1.488613728678813811542398d10,1.016803586272438228077304d11, &
    3.417476345507377132798597d11,4.463158187419713286462081d11/
!
!  Coefficients for minimax approximation over (12, INF).
!
  data c/-1.910444077728d-03,8.4171387781295d-04, &
       -5.952379913043012d-04,7.93650793500350248d-04, &
       -2.777777777777681622553d-03,8.333333333333333331554247d-02, &
    5.7083835261d-03/

  y = x

  if ((y > zero) .and. (y <= xbig)) then

    if (y <= eps) then

      res = -log(y)

    else if (y <= thrhal) then
!
!  EPS < X <= 1.5
!
      if (y < pnt68) then
        corr = -log(y)
        xm1 = y
      else
        corr = zero
        xm1 = (y - half) - half
      end if

      if ((y <= half) .or. (y >= pnt68)) then
        xden = one
        xnum = zero
        do i = 1, 8
          xnum = xnum*xm1 + p1(i)
          xden = xden*xm1 + q1(i)
        end do
        res = corr + (xm1 * (d1 + xm1*(xnum/xden)))
      else
        xm2 = (y - half) - half
        xden = one
        xnum = zero
        do i = 1, 8
          xnum = xnum*xm2 + p2(i)
          xden = xden*xm2 + q2(i)
        end do
        res = corr + xm2 * (d2 + xm2*(xnum/xden))
      end if

    else if (y <= four) then
!
!  1.5 < X <= 4.0
!
      xm2 = y - two
      xden = one
      xnum = zero
      do i = 1, 8
        xnum = xnum*xm2 + p2(i)
        xden = xden*xm2 + q2(i)
      end do
      res = xm2 * (d2 + xm2*(xnum/xden))

    else if (y <= twelve) then
!
!  4.0 < X <= 12.0
!
      xm4 = y - four
      xden = -one
      xnum = zero
      do i = 1, 8
        xnum = xnum*xm4 + p4(i)
        xden = xden*xm4 + q4(i)
      end do
      res = d4 + xm4*(xnum/xden)

    else
!
!  Evaluate for argument >= 12.0,
!
      res = zero

      if (y <= frtbig) then

        res = c(7)
        ysq = y * y
        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if
      res = res/y
      corr = log(y)
      res = res + sqrtpi - half*corr
      res = res + y*(corr-one)

    end if

  else
!
!  Return for bad arguments
!
    res = xinf

  end if

  dlgama = res

  return
end
subroutine dsubn ( x, nmax, xmax, d )

!*****************************************************************************80
!
!! DSUBN evaluates derivatives of Ei(X).
!
!  Discussion:
!
!    This routine is a translation of Gautschi's ACM Algorithm 282 for
!    derivatives of Ei(x).
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  logical ( kind = 4 ) bool1, bool2
  integer ( kind = 4 ) j,nmax,n0,mini,n,n1,lim
  real ( kind = 8 ) &
       b0,b1,b2,b3,b4,b5,b6,c0,c1,d,e,en,one,p,q,t,ten, &
       two,x,xmax,x1,z,zero
  dimension d(0:nmax)

  data zero/0.0d0/,one/1.0d0/,two/2.0d0/,ten/10.0d0/
  data c0/2.7183d0/,c1/4.67452d1/
  data b0/5.7941d-5/,b1/-1.76148d-3/,b2/2.08645d-2/, &
       b3/-1.29013d-1/,b4/8.5777d-1/,b5/1.0125d0/,b6/7.75d-1/

  x1 = abs(x)
  n0 = int(x1)
  e = exp(x)
  d(0) = e/x
  bool1 = (x < zero) .or. (x1 <= two)
  bool2 = n0 < nmax
  mini = min(n0,nmax)

  if (bool1) then
    lim = nmax
  else
    lim = mini
  end if

  n = 1
  en = one

50 continue

  d(n) = (e - en*d(n-1))/x
  n = n +1
  en = en + one

  if (x1 < one) then
    if ((abs(d(n-1)) < abs(xmax*x/en)) .and. (n <= lim)) then
      go to 50
    end if
  else
    if ((abs(d(n-1)/x) < xmax/en) .and. (n <= lim)) then
      go to 50
    end if
  end if

  do j = n, lim
    d(n) = zero
  end do

  if ((.not. bool1) .and. bool2) then

    t = (x1+c1)/(c0*x1)

    if (t < ten) then
      t = ((((b0*t + b1)*t + b2)*t + b3)*t + b4)*t + b5
    else
      z = log(t) - b6
      p = (b6-log(z))/(one+z)
      p = one/(one+p)
      t = t*p/z
    end if

    n1 = c0*x1*t - one
    if (n1 < nmax) n1 = nmax
    q = one/x
    en = one

    do n = 1,n1+1
      q = -en*q/x
      en = en+one
    end do

    do n = n1,n0+1,-1
      en = en - one
      q = (e-x*q)/en
      if (n <= nmax) d(n) = q
    end do

  end if

  return
end
function ei ( x )

!*****************************************************************************80
!
!! EI evaluates the exponential integral Ei(x).
!
!  Discussion:
!
!    This function computes approximate values for the
!    exponential integral Ei(x), where x is real.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) ei
  integer ( kind = 4 ) int
  real ( kind = 8 )  x, result

  int = 1
  call calcei ( x, result, int )
  ei = result

  return
end
function eone ( x )

!*****************************************************************************80
!
!! EONE evaluates the exponential integral E1(x).
!
!  Discussion:
!
!    This function computes approximate values for the
!    exponential integral E1(x), where  x  is real.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) eone
  integer ( kind = 4 ) int
  real ( kind = 8 )  x, result

  int = 2
  call calcei(x,result,int)
  eone = result

  return
end
function expei ( x )

!*****************************************************************************80
!
!! EXPEI evaluates exp(-x) * Ei(x).
!
!  Discussion:
!
!    This function computes approximate values for the
!    function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
!    integral, and  x  is real.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real ( kind = 8 ) expei
  integer ( kind = 4 ) int
  real ( kind = 8 ) x, result

  int = 3
  call calcei(x,result,int)
  expei = result

  return
end
subroutine machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
  maxexp, eps, epsneg, xmin, xmax )

!*****************************************************************************80
!
!! MACHAR determines the machine arithmetic parameters.
!
!  Discussion:
!
!    This routine is intended to determine the parameters
!    of the floating-point arithmetic system specified below.  The
!    determination of the first three uses an extension of an algorithm
!    due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
!    but not all, of the improvements suggested by M. Gentleman and S.
!    Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
!    program was published in the book Software Manual for the
!    Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
!    Englewood Cliffs, NJ, 1980.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameter values reported are as follows:
!
!   IBETA   - the radix for the floating-point representation
!   IT  - the number of base IBETA digits in the floating-point
!         significand
!   IRND    - 0 if floating-point addition chops
!         1 if floating-point addition rounds, but not in the
!       IEEE style
!         2 if floating-point addition rounds in the IEEE style
!         3 if floating-point addition chops, and there is
!       partial underflow
!         4 if floating-point addition rounds, but not in the
!       IEEE style, and there is partial underflow
!         5 if floating-point addition rounds in the IEEE style,
!       and there is partial underflow
!   NGRD    - the number of guard digits for multiplication with
!         truncating arithmetic.  It is
!         0 if floating-point arithmetic rounds, or if it
!       truncates and only  IT  base  IBETA digits
!       participate in the post-normalization shift of the
!       floating-point significand in multiplication;
!         1 if floating-point arithmetic truncates and more
!       than  IT  base  IBETA  digits participate in the
!       post-normalization shift of the floating-point
!       significand in multiplication.
!   MACHEP  - the largest negative integer such that
!         1.0+FLOAT(IBETA)**MACHEP /= 1.0, except that
!         MACHEP is bounded below by  -(IT+3)
!   NEGEPS  - the largest negative integer such that
!         1.0-FLOAT(IBETA)**NEGEPS /= 1.0, except that
!         NEGEPS is bounded below by  -(IT+3)
!   IEXP    - the number of bits (decimal places if IBETA = 10)
!         reserved for the representation of the exponent
!         (including the bias or sign) of a floating-point
!         number
!   MINEXP  - the largest in magnitude negative integer such that
!         FLOAT(IBETA)**MINEXP is positive and normalized
!   MAXEXP  - the smallest positive power of  BETA  that overflows
!   EPS     - FLOAT(IBETA)**MACHEP.
!   EPSNEG  - FLOAT(IBETA)**NEGEPS.
!   XMIN    - the smallest non-vanishing normalized floating-point
!         power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
!   XMAX    - the largest finite floating-point number.  In
!         particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
!         Note - on some machines  XMAX  will be only the
!         second, or perhaps third, largest number, being
!         too small by 1 or 2 units in the last digit of
!         the significand.
!
  implicit none

  integer ( kind = 4 ) i,ibeta,iexp,irnd,it,itemp,iz,j,k,machep,maxexp, &
      minexp,mx,negep,ngrd,nxres
  real ( kind = 8 ) &
     a,b,beta,betain,betah,conv,eps,epsneg,one,t,temp,tempa, &
     temp1,two,xmax,xmin,y,z,zero

  conv(i) = dble(i)
  one = conv(1)
  two = one + one
  zero = one - one
!
!  Determine IBETA, BETA ala Malcolm.
!
  a = one

  do

    a = a + a
    temp = a + one
    temp1 = temp - a

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  b = one

  do

    b = b + b
    temp = a+b
    itemp = int(temp-a)

    if ( itemp /= 0 ) then
      exit
    end if

  end do

  ibeta = itemp
  beta = conv(ibeta)
!
!  Determine IT, IRND.
!
  it = 0
  b = one

  do

    it = it + 1
    b = b * beta
    temp = b+one
    temp1 = temp-b

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  irnd = 0
  betah = beta / two
  temp = a + betah
  if ( temp - a /= zero ) then
    irnd = 1
  end if
  tempa = a + beta
  temp = tempa + betah
  if ( ( irnd == 0 ) .and. ( temp - tempa /= zero ) ) then
    irnd = 2
  end if
!
!  Determine NEGEP, EPSNEG.
!
  negep = it + 3
  betain = one / beta
  a = one
  do i = 1, negep
    a = a * betain
  end do

  b = a

  do
    temp = one - a
    if ( temp - one /= zero ) then
      exit
    end if
    a = a * beta
    negep = negep - 1
  end do

  negep = -negep
  epsneg = a
!
!  Determine MACHEP, EPS.
!
  machep = -it - 3
  a = b

  do

    temp = one+a
    if ( temp - one /= zero ) then
      exit
    end if
    a = a * beta
    machep = machep + 1
  end do

  eps = a
!
!  Determine NGRD.
!
  ngrd = 0
  temp = one + eps
  if ( ( irnd == 0 ) .and. ( temp * one - one /= zero ) ) then
    ngrd = 1
  end if
!
!  Determine IEXP, MINEXP, XMIN.
!
!  Loop to determine largest I and K = 2**I such that
!     (1/BETA) ** (2**(I))
!  does not underflow.
!  Exit from loop is signaled by an underflow.
!
  i = 0
  k = 1
  z = betain
  t = one + eps
  nxres = 0

  do

    y = z
    z = y * y
!
!  Check for underflow here.
!
    a = z * one
    temp = z * t

    if ( ( a + a == zero ) .or. ( abs ( z ) >= y ) ) then
      exit
    end if

    temp1 = temp * betain

    if ( temp1 * beta == z ) then
      exit
    end if

    i = i + 1
    k = k + k

  end do

  if ( ibeta /= 10 ) then

    iexp = i + 1
    mx = k + k

  else
!
!  This segment is for decimal machines only.
!
    iexp = 2
    iz = ibeta

    do

      if (k < iz) then
        exit
      end if

      iz = iz * ibeta
      iexp = iexp + 1

    end do

    mx = iz + iz - 1

  end if
!
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!
450 continue

  xmin = y
  y = y * betain
!
!  Check for underflow here.
!
  a = y * one
  temp = y * t
  if (((a+a) == zero) .or. (abs(y) >= xmin)) go to 460
  k = k + 1
  temp1 = temp * betain

  if ((temp1*beta /= y) .or. (temp == y)) then
    go to 450
  else
    nxres = 3
    xmin = y
  end if

460 continue

  minexp = -k
!
!  Determine MAXEXP, XMAX.
!
  if ( ( mx <= k + k - 3 ) .and. ( ibeta /= 10 ) ) then
    mx = mx + mx
    iexp = iexp + 1
  end if

  maxexp = mx + minexp
!
!  Adjust IRND to reflect partial underflow.
!
  irnd = irnd + nxres
!
!  Adjust for IEEE-style machines.
!
  if (irnd >= 2) maxexp = maxexp - 2
!
!  Adjust for machines with implicit leading bit in binary
!  significand, and machines with radix point at extreme
!  right of significand.
!
  i = maxexp + minexp
  if ((ibeta == 2) .and. (i == 0)) maxexp = maxexp - 1
  if (i > 20) maxexp = maxexp - 1
  if (a /= y) maxexp = maxexp - 2
  xmax = one - epsneg
  if ( xmax * one /= xmax ) then
    xmax = one - beta * epsneg
  end if
  xmax = xmax / ( beta * beta * beta * xmin )
  i = maxexp + minexp + 3

  do j = 1, i
    if (ibeta == 2) xmax = xmax + xmax
    if (ibeta /= 2) xmax = xmax * beta
  end do

  return
end
function psi ( xx )

!*****************************************************************************80
!
!! PSI evaluates the Psi function.
!
!  Discussion:
!
!    This function evaluates the logarithmic derivative of the gamma function,
!      psi(x) = d/dx (gamma(x)) / gamma(x) = d/dx (ln gamma(x))
!    for real x, where either
!     -xmax1 < x < -xmin (x not a negative integer), or
!     xmin < x.
!
!    The calling sequence for this function is
!      Y = PSI(X)
!    The main computation uses rational Chebyshev approximations
!    published in Math. Comp. 27, 123-127 (1973) by Cody, Strecok and
!    Thacher.  This transportable program is patterned after the
!    machine-dependent FUNPACK program PSI(X), but cannot match that
!    version for efficiency or accuracy.  This version uses rational
!    approximations that are theoretically accurate to 20 significant
!    decimal digits.  The accuracy achieved depends on the arithmetic
!    system, the compiler, the intrinsic functions, and proper selection
!    of the machine-dependent constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.
!
!   XINF   = largest positive machine number
!   XMAX1  = beta ** (p-1), where beta is the radix for the
!    floating-point system, and p is the number of base-beta
!    digits in the floating-point significand.  This is an
!    upper bound on non-integral floating-point numbers, and
!    the negative of the lower bound on acceptable negative
!    arguments for PSI.  If rounding is necessary, round this
!    value down.
!   XMIN1  = the smallest in magnitude acceptable argument.  We
!    recommend XMIN1 = MAX(1/XINF,xmin) rounded up, where
!    xmin is the smallest positive floating-point number.
!   XSMALL = absolute argument below which  PI*COTAN(PI*X)  may be
!    represented by 1/X.  We recommend XSMALL < sqrt(3 eps)/pi,
!    where eps is the smallest positive number such that
!    1+eps > 1.
!   XLARGE = argument beyond which PSI(X) may be represented by
!    LOG(X).  The solution to the equation
!       x*ln(x) = beta ** p
!    is a safe value.
!
! Error Returns
!
!  The program returns XINF for  X < -XMAX1, for X zero or a negative
!    integer ( kind = 4 ), or when X lies in (-XMIN1, 0), and returns -XINF
!    when X lies in (0, XMIN1).
!
  implicit none

  real ( kind = 8 ) psi
  integer ( kind = 4 ) i,n,nq
  real ( kind = 8 ) &
     aug,conv,den,four,fourth,half,one,p1,p2,piov4,q1,q2, &
     sgn,three,xlarge,upper,w,x,xinf,xmax1,xmin1,xsmall,x01, &
     x01d,x02,xx,z,zero
  dimension p1(9),p2(7),q1(8),q2(6)
!
!  Mathematical constants.  PIOV4 = pi / 4
!
  data zero,fourth,half,one/0.0d0,0.25d0,0.5d0,1.0d0/
  data three,four/3.0d0,4.0d0/,piov4/7.8539816339744830962d-01/
!
!  Machine-dependent constants
!
  data xinf/1.79d+308/, xmin1/2.23d-308/, xmax1/4.50d+15/, &
       xsmall/5.80d-09/, xlarge/2.71d+14/
!
!  Zero of psi(x)
!
  data x01/187.0d0/,x01d/128.0d0/,x02/6.9464496836234126266d-04/
!
!  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
!
  data p1/4.5104681245762934160d-03,5.4932855833000385356d+00, &
      3.7646693175929276856d+02,7.9525490849151998065d+03, &
      7.1451595818951933210d+04,3.0655976301987365674d+05, &
      6.3606997788964458797d+05,5.8041312783537569993d+05, &
      1.6585695029761022321d+05/
  data q1/9.6141654774222358525d+01,2.6287715790581193330d+03, &
      2.9862497022250277920d+04,1.6206566091533671639d+05, &
      4.3487880712768329037d+05,5.4256384537269993733d+05, &
      2.4242185002017985252d+05,6.4155223783576225996d-08/
!
!  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
!     for  x > 3.0
!
  data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01, &
      -1.9784554148719218667d+01,-8.8100958828312219821d+00, &
      -1.4479614616899842986d+00,-7.3689600332394549911d-02, &
      -6.5135387732718171306d-21/
  data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02, &
       2.4736979003315290057d+02, 1.0742543875702278326d+02, &
       1.7463965060678569906d+01, 8.8427520398873480342d-01/

  conv(i) = dble(i)
  x = xx
  w = abs(x)
  aug = zero
!
!  Check for valid arguments, then branch to appropriate algorithm
!
  if ((-x >= xmax1) .or. (w < xmin1)) then

    psi = xinf
    if (x > zero) psi = -xinf
    return

  else if (x >= half) then

    go to 200
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
  else if (w <= xsmall) then

    aug = -one / x
    go to 150

  end if
!
!  Argument reduction for cot
!
  if (x < zero) then
    sgn = piov4
  else
    sgn = -piov4
  end if

  w = w - aint(w)
  nq = int(w * four)
  w = four * (w - conv(nq) * fourth)
!
!  W is now related to the fractional part of  4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
  n = nq / 2
  if ((n+n) /= nq) w = one - w
  z = piov4 * w
  if (mod(n,2) /= 0) sgn = - sgn
!
!  Determine the final value for  -pi * cotan(pi*x)
!
  n = (nq + 1) / 2

  if (mod(n,2) == 0) then
!
!  Check for singularity
!
    if (z == zero) then
      psi = xinf
      if (x > zero) psi = -xinf
      return
    end if

    aug = sgn * (four / tan(z))

  else

    aug = sgn * (four * tan(z))

  end if

150 continue

  x = one - x

200 continue
!
!  0.5 <= X <= 3.0
!
  if (x <= three) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = (den + q1(i)) * x
      upper = (upper + p1(i+1)) * x
    end do
    den = (upper + p1(9)) / (den + q1(8))
    x = (x-x01/x01d) - x02
    psi = den * x + aug
    return

  end if
!
!  3.0 < X
!
  if (x < xlarge) then

    w = one / (x * x)
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = (den + q2(i)) * w
      upper = (upper + p2(i+1)) * w
    end do
    aug = (upper + p2(7)) / (den + q2(6)) - half / x + aug

  end if

  psi = aug + log(x)

  return
end
function ren ( k )

!*****************************************************************************80
!
!! REN is a random number generator.
!
!  Discussion:
!
!    This function is a random number generator - based on 
!
!    This function is intended for use on computers with
!    fixed point wordlength of at least 29 bits.  It is
!    best if the floating-point significand has at most 29 bits.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Algorithm 266
!    Pike and Hill (modified by Hansson), 
!    Communications of the ACM,
!    Vol. 8, No. 10, October 1965.
!
  implicit none

  integer ( kind = 4 ) iy,j,k
  real ( kind = 8 ) conv,c1,c2,c3,one
  real ( kind = 8 ) ren

  data iy/100001/
  data one,c1,c2,c3/1.0d0,2796203.0d0,1.0d-6,1.0d-12/
!
!  Statement functions for conversion between integer and float
!
  conv(j) = dble(j)

  j = k
  iy = iy * 125
  iy = iy - (iy/2796203) * 2796203
  ren = conv(iy) / c1 * (one + c2 + c3)

  return
end
subroutine ribesl ( x, alpha, nb, ize, b, ncalc )

!*****************************************************************************80
!
!! RIBESL evaluates a sequence of Bessel I functions.
!
!  Discussion:
!
!    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA,
!    with or without exponential scaling.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
! X     - Working precision non-negative real argument for which
!     I's or exponentially scaled I's (I*EXP(-X))
!     are to be calculated.  If I's are to be calculated,
!     X must be less than EXPARG (see below).
!
! ALPHA - Working precision fractional part of order for which
!     I's or exponentially scaled I's (I*EXP(-X)) are
!     to be calculated.  0 <= ALPHA < 1.0.
!
! NB    - Integer number of functions to be calculated, NB > 0.
!     The first function calculated is of order ALPHA, and the
!     last is of order (NB - 1 + ALPHA).
!
! IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
!     and 2 if exponentially scaled I's are to be calculated.
!
! B     - Working precision output vector of length NB.  If the routine
!     terminates normally (NCALC=NB), the vector B contains the
!   functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
!     corresponding exponentially scaled functions.
!
! NCALC - Integer output variable indicating possible errors.
!     Before using the vector B, the user should check that
!     NCALC=NB, i.e., all orders have been calculated to
!     the desired accuracy.  See error returns below.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   it     = Number of bits in the mantissa of a working precision
!    variable
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   NSIG   = Decimal significance desired.  Should be set to
!    INT(LOG10(2)*it+1).  Setting NSIG lower will result
!    in decreased accuracy while setting NSIG higher will
!    increase CPU time without increasing accuracy.  The
!    truncation error is limited to a relative error of
!    T=.5*10**(-NSIG).
!   ENTEN  = 10.0 ** K, where K is the largest integer such that
!    ENTEN is machine-representable in working precision
!   ENSIG  = 10.0 ** NSIG
!   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
!    K >= NSIG/4
!   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
!   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
!    in mind that if ABS(X)=N, then at least N iterations
!    of the backward recursion will be executed.  The value
!    of 10.0 ** 4 is used on every machine.
!   EXPARG = Largest working precision argument that the library
!    EXP routine can handle and upper limit on the
!    magnitude of X when IZE=1; approximately
!    LOG(beta**maxexp)
!
! Error returns
!
!  In case of an error,  NCALC /= NB, and not all I's are
!  calculated to the desired accuracy.
!
!  NCALC < 0:  An argument is out of range. For example,
!     NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= EXPARG.
!     In this case, the B-vector is not calculated, and NCALC is
!     set to MIN0(NB,0)-1 so that NCALC /= NB.
!
!  NB > NCALC > 0: Not all requested function values could
!     be calculated accurately.  This usually occurs because NB is
!     much larger than ABS(X).  In this case, B(N) is calculated
!     to the desired accuracy for N <= NCALC, but precision
!     is lost for NCALC < N <= NB.  If B(N) does not vanish
!     for N > NCALC (because it is too small to be represented),
!     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!     significant figures of B(N) can be trusted.
!
!  Acknowledgement:
!
!    This program is based on a program written by David J.
!    Sookne (2) that computes values of the Bessel functions J or
!    I of real argument and integer order.  Modifications include
!    the restriction of the computation to the I Bessel function
!    of non-negative real argument, the extension of the computation
!    to arbitrary positive order, the inclusion of optional
!    exponential scaling, and the elimination of most underflow.
!    An earlier version was published in (3).
!
!  Reference: 
!
!    "A Note on Backward Recurrence Algorithms," Olver,
!      F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
!      pp 941-947.
!
!     "Bessel Functions of Real Argument and Integer Order,"
!      Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
!      125-132.
!
!     "ALGORITHM 597, Sequence of Modified Bessel Functions
!      of the First Kind," Cody, W. J., Trans. Math. Soft.,
!      1983, pp. 242-245.
!
  implicit none

  integer ( kind = 4 ) ize,k,l,magx,n,nb,nbmx,ncalc,nend,nsig,nstart
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b,const,conv,em,empal,emp2al,en,enmten,ensig, &
   enten,exparg,func,half,halfx,one,p,plast,pold,psave,psavel, &
   rtnsig,sum,tempa,tempb,tempc,test,tover,two,x,xlarge,zero
  dimension b(nb)
!
!  Mathematical constants
!
  data one,two,zero,half,const/1.0d0,2.0d0,0.0d0,0.5d0,1.585d0/
!
!  Machine-dependent parameters
!
  data nsig,xlarge,exparg /16,1.0d4,709.0d0/
  data enten,ensig,rtnsig/1.0d308,1.0d16,1.0d-4/
  data enmten/8.9d-308/
!
!  Statement functions for conversion
!
  conv(n) = dble(n)
!
!  Check for X, NB, OR IZE out of range.
!
  if ( &
    ( nb > 0 ) .and. &
    (x >= zero) .and. &
    (alpha >= zero) .and. &
    (alpha < one) .and. &
    (((ize == 1) .and. (x <= exparg)) .or. &
    ((ize == 2) .and. (x <= xlarge)))) then
!
!  Use 2-term ascending series for small X
!
    ncalc = nb
    magx = int(x)

    if (x >= rtnsig) then
!
!  Initialize the forward sweep, the P-sequence of Olver
!
      nbmx = nb - magx
      n = magx+1
      en = conv(n+n) + (alpha+alpha)
      plast = one
      p = en / x
!
! Calculate general significance test
!
      test = ensig + ensig

      if (2*magx > 5*nsig) then
        test = sqrt(test*p)
      else
        test = test / const**magx
      end if

      if (nbmx >= 3) then
!
! Calculate P-sequence until N = NB-1.  Check for possible overflow.
!
        tover = enten / ensig
        nstart = magx+2
        nend = nb - 1

        do k = nstart, nend
          n = k
          en = en + two
          pold = plast
          plast = p
          p = en * plast/x + pold
          if (p > tover) then
!
!  To avoid overflow, divide P-sequence by TOVER.  Calculate
!  P-sequence until ABS(P) > 1.
!
            tover = enten
            p = p / tover
            plast = plast / tover
            psave = p
            psavel = plast
            nstart = n + 1

   60       continue

            n = n + 1
            en = en + two
            pold = plast
            plast = p
            p = en * plast/x + pold
            if (p <= one) go to 60

            tempb = en / x
!
!  Calculate backward test, and find NCALC, the highest N
!  such that the test is passed.
!
            test = pold*plast / ensig
            test = test*(half-half/(tempb*tempb))
            p = plast * tover
            n = n - 1
            en = en - two
            nend = min0(nb,n)
            do l = nstart, nend
              ncalc = l
              pold = psavel
              psavel = psave
              psave = en * psavel/x + pold
              if (psave*psavel > test) go to 90
            end do

            ncalc = nend + 1

   90       continue

            ncalc = ncalc - 1
            go to 120
          end if

        end do

        n = nend
        en = conv(n+n) + (alpha+alpha)
!
!  Calculate special significance test for NBMX > 2.
!
        test = max(test,sqrt(plast*ensig)*sqrt(p+p))

      end if
!
!  Calculate P-sequence until significance test passed.
!
110   continue

      n = n + 1
      en = en + two
      pold = plast
      plast = p
      p = en * plast/x + pold
      if (p < test) go to 110
!
!  Initialize the backward recursion and the normalization sum.
!
120   continue

      n = n + 1
      en = en + two
      tempb = zero
      tempa = one / p
      em = conv(n) - one
      empal = em + alpha
      emp2al = (em - one) + (alpha + alpha)
      sum = tempa * empal * emp2al / em
      nend = n - nb

      if (nend < 0) then
!
!  N < NB, so store B(N) and set higher orders to zero.
!
        b(n) = tempa
        nend = -nend
        do l = 1, nend
          b(n+l) = zero
        end do

      else

        if (nend > 0) then
!
!  Recur backward via difference equation, calculating (but
!  not storing) B(N), until N = NB.
!
          do l = 1, nend
            n = n - 1
            en = en - two
            tempc = tempb
            tempb = tempa
            tempa = (en*tempb) / x + tempc
            em = em - one
            emp2al = emp2al - one
            if (n == 1) then
              exit
            end if
            if (n == 2) emp2al = one
            empal = empal - one
            sum = (sum + tempa*empal) * emp2al / em
          end do

        end if
!
!  Store B(NB)
!
        b(n) = tempa

        if (nb <= 1) then
          sum = (sum + sum) + tempa
          go to 230
        end if
!
!  Calculate and Store B(NB-1)
!
        n = n - 1
        en = en - two
        b(n)  = (en*tempa) / x + tempb
        if (n == 1) go to 220
        em = em - one
        emp2al = emp2al - one
        if (n == 2) emp2al = one
        empal = empal - one
        sum = (sum + b(n)*empal) * emp2al / em

      end if

      nend = n - 2

      if (nend > 0) then
!
!  Calculate via difference equation and store B(N), until N = 2.
!
        do l = 1, nend
          n = n - 1
          en = en - two
          b(n) = (en*b(n+1)) / x +b(n+2)
          em = em - one
          emp2al = emp2al - one
          if (n == 2) emp2al = one
          empal = empal - one
          sum = (sum + b(n)*empal) * emp2al / em
        end do

      end if
!
!  Calculate B(1)
!
      b(1) = two*empal*b(2) / x + b(3)

220   continue

      sum = (sum + sum) + b(1)
!
!  Normalize.  Divide all B(N) by sum.
!
230   continue

      if (alpha /= zero) then
         sum = sum * gamma ( one + alpha ) * (x*half)**(-alpha)
      end if

      if (ize == 1) sum = sum * exp(-x)
      tempa = enmten
      if (sum > one) tempa = tempa * sum

      do n = 1, nb
        if (b(n) < tempa) b(n) = zero
        b(n) = b(n) / sum
      end do

      return
!
!  Two-term ascending series for small X.
!
    else

      tempa = one
      empal = one + alpha
      halfx = zero
      if (x > enmten) halfx = half * x
      if (alpha /= zero) tempa = halfx**alpha / gamma ( empal )
      if (ize == 2) tempa = tempa * exp(-x)
      tempb = zero
      if ((x+one) > one) tempb = halfx * halfx
      b(1) = tempa + tempa*tempb / empal
      if ((x /= zero) .and. (b(1) == zero)) ncalc = 0
 
      if (nb > 1) then

        if (x == zero) then

          do n = 2, nb
            b(n) = zero
          end do

        else
!
!  Calculate higher-order functions.
!
          tempc = halfx
          tover = (enmten + enmten) / x
          if (tempb /= zero) tover = enmten / tempb

          do n = 2, nb
            tempa = tempa / empal
            empal = empal + one
            tempa = tempa * tempc
            if (tempa <= tover*empal) tempa = zero
            b(n) = tempa + tempa*tempb / empal
            if ((b(n) == zero) .and. (ncalc > n)) then
              ncalc = n-1
            end if
          end do

        end if

      end if

    end if

  else

    ncalc = min0(nb,0)-1

  end if

  return
end
subroutine rjbesl ( x, alpha, nb, b, ncalc )

!*****************************************************************************80
!
!! RJBESL evaluates a sequence of Bessel J functions.
!
!  Discussion:
!
!    This routine calculates Bessel functions J sub(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA.
!
!    This program is based on a program written by David Sookne
!    that computes values of the Bessel functions J or I of real
!    argument and integer order.  Modifications include the restriction
!    of the computation to the J Bessel function of non-negative real
!    argument, the extension of the computation to arbitrary positive
!    order, and the elimination of most underflow.
!
!  Modified:
!
!    15 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference: 
!
!    F Olver, David Sookne,
!    A Note on Backward Recurrence Algorithms," 
!    Math. Comp.,
!    Volume 26, 1972, pages 941-947.
!
!    David Sookne,
!    Bessel Functions of Real Argument and Integer Order,
!    NBS Journal of Res. B,
!    Volume 77B, 1973, pages 125-132.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, non-negative real argument for which
!    J's are to be calculated.
!
!    Input, real ( kind = 8 ) ALPHA, fractional part of order for which
!    J's or exponentially scaled J'r (J*exp(X)) are
!    to be calculated.  0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, number of functions to be calculated, 
!    NB > 0.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Output, real ( kind = 8 ) B(NB).  If RJBESL
!    terminates normally (NCALC=NB), the vector B contains the
!    functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
!    corresponding exponentially scaled functions.
!
!    Output, integer ( kind = 4 ) NCALC, indicates possible errors.
!    Before using the vector B, the user should check that
!    NCALC=NB, i.e., all orders have been calculated to
!    the desired accuracy.  See Error Returns below.
!
!  Internal Parameters:
!
!    IT = Number of bits in the mantissa of a working precision variable
!
!    NSIG   = Decimal significance desired.  Should be set to
!    INT(LOG10(2)*it+1).  Setting NSIG lower will result
!    in decreased accuracy while setting NSIG higher will
!    increase CPU time without increasing accuracy.  The
!    truncation error is limited to a relative error of
!    T=.5*10**(-NSIG).
!
!    Then the following machine-dependent constants must be declared
!    in DATA statements.  IEEE values are provided as a default.
!
!    ENTEN  = 10.0 ** K, where K is the largest integer such that
!    ENTEN is machine-representable in working precision.
!
!    ENSIG  = 10.0 ** NSIG
!
!    RTNSIG = 10.0 ** (-K) for the smallest integer K such that K >= NSIG/4
!
!    ENMTEN = Smallest ABS(X) such that X/4 does not underflow
!
!    XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
!    then at least N iterations of the backward recursion
!    will be executed.  The value of 10.0 ** 4 is used on
!    every machine.
!
!  Error returns:
!
!    In case of an error,  NCALC /= NB, and not all J's are
!    calculated to the desired accuracy.
!
!    NCALC < 0:  An argument is out of range. For example,
!    NBES <= 0, ALPHA < 0 or > 1, or X is too large.
!    In this case, B(1) is set to zero, the remainder of the
!    B-vector is not calculated, and NCALC is set to
!    MIN(NB,0)-1 so that NCALC /= NB.
!
!    NB > NCALC > 0: Not all requested function values could
!    be calculated accurately.  This usually occurs because NB is
!    much larger than ABS(X).  In this case, B(N) is calculated
!    to the desired accuracy for N <= NCALC, but precision
!    is lost for NCALC < N <= NB.  If B(N) does not vanish
!    for N > NCALC (because it is too small to be represented),
!    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!    significant figures of B(N) can be trusted.
!
  implicit none

  integer ( kind = 4 ) nb

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpem
  real ( kind = 8 ) alp2em
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) capp
  real ( kind = 8 ) capq
  real ( kind = 8 ) eighth
  real ( kind = 8 ) em
  real ( kind = 8 ) en
  real ( kind = 8 ) enmten
  real ( kind = 8 ) ensig
  real ( kind = 8 ) enten
  real ( kind = 8 ) fact(25)
  real ( kind = 8 ) four
  real ( kind = 8 ) gnu
  real ( kind = 8 ) half
  real ( kind = 8 ) halfx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) magx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbmx
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) nend
  integer ( kind = 4 ) nstart
  real ( kind = 8 ) one
  real ( kind = 8 ) one30
  real ( kind = 8 ) p
  real ( kind = 8 ) pi2
  real ( kind = 8 ) plast
  real ( kind = 8 ) pold
  real ( kind = 8 ) psave
  real ( kind = 8 ) psavel
  real ( kind = 8 ) rtnsig
  real ( kind = 8 ) s
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) tempa
  real ( kind = 8 ) tempb
  real ( kind = 8 ) tempc
  real ( kind = 8 ) test
  real ( kind = 8 ) three
  real ( kind = 8 ) three5
  real ( kind = 8 ) tover
  real ( kind = 8 ) two
  real ( kind = 8 ) twofiv
  real ( kind = 8 ) twopi1
  real ( kind = 8 ) twopi2
  real ( kind = 8 ) x
  real ( kind = 8 ) xc
  real ( kind = 8 ) xin
  real ( kind = 8 ) xk
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xm
  real ( kind = 8 ) vcos
  real ( kind = 8 ) vsin
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
!
!  Mathematical constants
!
!  PI2    - 2 / PI
!  TWOPI1 - first few significant digits of 2 * PI
!  TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
!  TWOPI1 + TWOPI2 = 2 * PI to extra precision.
!
  data pi2 / 0.636619772367581343075535d0 /
  data twopi1 / 6.28125d0 /
  data twopi2 / 1.935307179586476925286767d-3 /
  data zero / 0.0d0 /
  data eighth / 0.125d0 /
  data half / 0.5d0 /
  data one / 1.0d0 /
  data two / 2.0d0 /
  data three / 3.0d0 /
  data four / 4.0d0 /
  data twofiv / 25.0d0 /
  data one30 / 130.0d0 /
  data three5 / 35.0d0 /
!
!  Machine-dependent parameters
!
  data enten / 1.0D+308 /
  data ensig / 1.0D+16 /
  data rtnsig / 1.0D-04 /
  data enmten / 8.90D-308 /
  data xlarge / 1.0D+04 /
!
!  Factorial(N)
!
  data fact / &
    1.0d0, &
    1.0d0, &
    2.0d0, &
    6.0d0, &
    24.0d0, &
    1.2d2, &
    7.2d2, &
    5.04d3, &
    4.032d4, &
    3.6288d5, &
    3.6288d6, &
    3.99168d7, &
    4.790016d8, &
    6.2270208d9, &
    8.71782912d10, &
    1.307674368d12, &
    2.0922789888d13, &
    3.55687428096d14, &
    6.402373705728d15, &
    1.21645100408832d17, &
    2.43290200817664d18, &
    5.109094217170944d19, &
    1.12400072777760768d21, &
    2.585201673888497664d22, &
    6.2044840173323943936d23 /

  jump = .false.
!
!  Check for out of range arguments.
!
  magx = int ( x )

  if ( &
    ( 0 < nb ) .and. &
    ( zero <= x ) .and. &
    ( x <= xlarge ) .and. &
    ( zero <= alpha ) .and. &
    ( alpha < one ) ) then
!
!  Initialize result array to zero.
!
    ncalc = nb
    b(1:nb) = zero
!
!  Branch to use 2-term ascending series for small X and asymptotic
!  form for large X when NB is not too large.
!
    if ( x < rtnsig ) then
!
!  Two-term ascending series for small X.
!
      tempa = one
      alpem = one + alpha

      if ( enmten < x ) then
        halfx = half * x
      else
        halfx = zero
      end if

      if ( alpha /= zero ) then
        tempa = halfx ** alpha / ( alpha * gamma ( alpha ) )
      end if

      if ( one < ( x + one ) ) then
        tempb = - halfx * halfx
      else
        tempb = zero
      end if

      b(1) = tempa + tempa * tempb / alpem

      if ( ( x /= zero ) .and. ( b(1) == zero ) ) then
        ncalc = 0
      end if

      if ( nb /= 1 ) then

        if ( x <= zero ) then

          do n = 2, nb
            b(n) = zero
          end do

        else
!
!  Calculate higher order functions.
!
          tempc = halfx

          if ( tempb /= zero ) then
            tover = enmten / tempb
          else
            tover = ( enmten + enmten ) / x
          end if

          do n = 2, nb

            tempa = tempa / alpem
            alpem = alpem + one

            tempa = tempa * tempc
            if ( tempa <= tover * alpem ) then
              tempa = zero
            end if

            b(n) = tempa + tempa * tempb / alpem

            if ( ( b(n) == zero ) .and. ( n < ncalc ) ) then
              ncalc = n - 1
            end if

          end do

        end if

      end if

    else if ( ( twofiv < x ) .and. ( nb <= magx + 1 ) ) then
!
!  Asymptotic series for 21.0 < X.
!
      xc = sqrt ( pi2 / x )
      xin = ( eighth / x ) ** 2

      if ( x < three5 ) then
        m = 11
      else if ( x < one30 ) then
        m = 8
      else
        m = 4
      end if

      xm = four * real ( m, kind = 8 )
!
!  Argument reduction for SIN and COS routines.
!
      t = aint ( x / ( twopi1 + twopi2 ) + half )
      z = ( ( x - t * twopi1 ) - t * twopi2 ) - ( alpha + half ) / pi2
      vsin = sin ( z )
      vcos = cos ( z )
      gnu = alpha + alpha

      do i = 1, 2

        s = ( ( xm - one ) - gnu ) * ( ( xm - one ) + gnu ) * xin * half
        t = ( gnu - ( xm - three ) ) * ( gnu + ( xm - three ) )
        capp = s * t / fact(2*m+1)
        t1 = ( gnu - ( xm + one ) ) * ( gnu + ( xm + one ) )
        capq = s * t1 / fact(2*m+2)
        xk = xm
        k = m + m
        t1 = t

        do j = 2, m
          xk = xk - four
          s = ( ( xk - one ) - gnu ) * ( ( xk - one ) + gnu )
          t = ( gnu - ( xk - three ) ) * ( gnu + ( xk - three ) )
          capp = ( capp + one / fact(k-1) ) * s * t * xin
          capq = ( capq + one / fact(k) ) * s * t1 * xin
          k = k - 2
          t1 = t
        end do

        capp = capp + one
        capq = ( capq + one ) * ( gnu * gnu - one ) * ( eighth / x )
        b(i) = xc * ( capp * vcos - capq * vsin )

        if ( nb == 1 ) then
          return
        end if

        t = vsin
        vsin = - vcos
        vcos = t
        gnu = gnu + two

      end do
!
!  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1
!
      if ( 2 < nb ) then
        gnu = alpha + alpha + two
        do j = 3, nb
          b(j) = gnu * b(j-1) / x - b(j-2)
          gnu = gnu + two
        end do
      end if
!
!  Use recurrence to generate results.  First initialize the
!  calculation of P*S.
!
    else

      nbmx = nb - magx
      n = magx + 1
      en = real ( n + n, kind = 8 ) + ( alpha + alpha )
      plast = one
      p = en / x
!
!  Calculate general significance test.
!
      test = ensig + ensig

      if ( 3 <= nbmx ) then
!
!  Calculate P*S until N = NB-1.  Check for possible overflow.
!
        tover = enten / ensig
        nstart = magx + 2
        nend = nb - 1
        en = real ( nstart + nstart, kind = 8 ) - two + ( alpha + alpha )

        do k = nstart, nend

          n = k
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x - pold

          if ( tover < p ) then
!
!  To avoid overflow, divide P*S by TOVER.  Calculate P*S until 1 < ABS(P).
!
            tover = enten
            p = p / tover
            plast = plast / tover
            psave = p
            psavel = plast
            nstart = n + 1

            do

              n = n + 1
              en = en + two
              pold = plast
              plast = p
              p = en * plast / x - pold
              if ( one < p ) then
                exit
              end if

            end do

            tempb = en / x
!
!  Calculate backward test and find NCALC, the highest N such that
!  the test is passed.
!
            test = pold * plast * ( half - half / ( tempb * tempb ) )
            test = test / ensig
            p = plast * tover
            n = n - 1
            en = en - two
            nend = min ( nb, n )

            do l = nstart, nend
              pold = psavel
              psavel = psave
              psave = en * psavel / x - pold
              if ( test < psave * psavel ) then
                ncalc = l - 1
                jump = .true.
                exit
              end if
            end do

            if ( jump ) then
              exit
            end if

            ncalc = nend
            jump = .true.
            exit

          end if

        end do

        if ( .not. jump ) then

          n = nend
          en = real ( n + n, kind = 8 ) + ( alpha + alpha )
!
!  Calculate special significance test for 2 < NBMX.
!
          test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

        end if

      end if
!
!  Calculate P*S until significance test passes.
!
      if ( .not. jump ) then

        do

          n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x - pold

          if ( test <= p ) then
            exit
          end if

        end do

      end if
!
!  Initialize the backward recursion and the normalization sum.
!
      n = n + 1
      en = en + two
      tempb = zero
      tempa = one / p
      m = 2 * n - 4 * ( n / 2 )
      sum = zero
      em = real ( n / 2, kind = 8 )
      alpem = ( em - one ) + alpha
      alp2em = ( em + em ) + alpha
      if ( m /= 0 ) then
        sum = tempa * alpem * alp2em / em
      end if
      nend = n - nb

      if ( 0 < nend ) then
!
!  Recur backward via difference equation, calculating (but not
!  storing) B(N), until N = NB.
!
        do l = 1, nend

          n = n - 1
          en = en - two
          tempc = tempb
          tempb = tempa
          tempa = ( en * tempb ) / x - tempc
          m = 2 - m

          if ( m /= 0 ) then
            em = em - one
            alp2em = ( em + em ) + alpha
            if ( n == 1 ) then
              exit
            end if
            alpem = ( em - one ) + alpha
            if ( alpem == zero ) then
              alpem = one
            end if
            sum = ( sum + tempa * alp2em ) * alpem / em
          end if

        end do

      end if
!
!  Store B(NB).
!
      b(n) = tempa

      if ( 0 <= nend ) then

        if ( nb <= 1 ) then

          alp2em = alpha
          if ( ( alpha + one ) == one ) then
            alp2em = one
          end if
          sum = sum + b(1) * alp2em

          if ( ( alpha + one ) /= one ) then
            sum = sum * gamma ( alpha ) * ( x * half ) ** ( - alpha )
          end if

          tempa = enmten

          if ( one < sum ) then
            tempa = tempa * sum
          end if

          do n = 1, nb
            if ( abs ( b(n) ) < tempa ) then
              b(n) = zero
            end if
            b(n) = b(n) / sum
          end do

          return

        else
!
!  Calculate and store B(NB-1).
!
          n = n - 1
          en = en - two
          b(n) = ( en * tempa ) / x - tempb

          if ( n == 1 ) then

            em = em - one
            alp2em = ( em + em ) + alpha
            if ( alp2em == zero ) then
              alp2em = one
            end if
            sum = sum + b(1) * alp2em
!
!  Normalize.  Divide all B(N) by sum.
!
            if ( ( alpha + one ) /= one ) then
              sum = sum * gamma ( alpha ) * ( x * half ) ** ( - alpha )
            end if

            tempa = enmten

            if ( one < sum ) then
              tempa = tempa * sum
            end if
  
            do n = 1, nb
              if ( abs ( b(n) ) < tempa ) then
                b(n) = zero
              end if
              b(n) = b(n) / sum
            end do

            return

          end if

          m = 2 - m

          if ( m /= 0 ) then
            em = em - one
            alp2em = ( em + em ) + alpha
            alpem = ( em - one ) + alpha
            if ( alpem == zero ) then
              alpem = one
            end if
            sum = ( sum + b(n) * alp2em ) * alpem / em
          end if

        end if

      end if

      nend = n - 2

      if ( nend /= 0 ) then
!
!  Calculate via difference equation and store B(N), until N = 2.
!
        do l = 1, nend

          n = n - 1
          en = en - two
          b(n) = ( en * b(n+1) ) / x - b(n+2)
          m = 2 - m

          if ( m /= 0 ) then
            em = em - one
            alp2em = ( em + em ) + alpha
            alpem = ( em - one ) + alpha
            if ( alpem == zero ) then
              alpem = one
            end if
            sum = ( sum + b(n) * alp2em ) * alpem / em
          end if

        end do

      end if
!
!  Calculate B(1).
!
      b(1) = two * ( alpha + one ) * b(2) / x - b(3)

      em = em - one
      alp2em = ( em + em ) + alpha
      if ( alp2em == zero ) then
        alp2em = one
      end if
      sum = sum + b(1) * alp2em
!
!  Normalize.  Divide all B(N) by sum.
!
      if ( ( alpha + one ) /= one ) then
        sum = sum * gamma ( alpha ) * ( x * half ) ** ( - alpha )
      end if

      tempa = enmten

      if ( one < sum ) then
        tempa = tempa * sum
      end if

      do n = 1, nb
        if ( abs ( b(n) ) < tempa ) then
          b(n) = zero
        end if
        b(n) = b(n) / sum
      end do

    end if
!
!  Error return.
!
  else

    b(1) = zero
    ncalc = min ( nb, 0 ) - 1

  end if

  return
end
subroutine rkbesl ( x, alpha, nb, ize, bk, ncalc )

!*****************************************************************************80
!
!! RKBESL evaluates a sequence of Bessel K functions.
!
!  Discussion:
!
!    This routine calculates modified Bessel functions
!    of the second kind, K SUB(N+ALPHA) (X), for non-negative
!    argument X, and non-negative order N+ALPHA, with or without
!    exponential scaling.
!
!    This program is based on a program written by J. B. Campbell
!    that computes values of the Bessel functions K of real
!    argument and real order.  Modifications include the addition
!    of non-scaled functions, parameterization of machine
!    dependencies, and the use of more accurate approximations
!    for SINH and SIN.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference: 
!
!    J B Campbell,
!    On Temme's Algorithm for the Modified Bessel functions of the Third Kind,
!    ACM TOMS,
!    Volume 6, Number 4, December 1980, pages 581-586.
!
!    J B Campbell,
!    A FORTRAN IV Subroutine for the Modified Bessel functions of the Third 
!    Kind of Real Order and Real Argument,
!    Report NRC/ERB-925,
!    National Research Council, Canada.
!
!  Parameters:
!
! X     - Working precision non-negative real argument for which
!     K's or exponentially scaled K's (K*EXP(X))
!     are to be calculated.  If K's are to be calculated,
!     X must not be greater than XMAX (see below).
!
! ALPHA - Working precision fractional part of order for which
!     K's or exponentially scaled K's (K*EXP(X)) are
!     to be calculated.  0 <= ALPHA < 1.0.
!
! NB    - Integer number of functions to be calculated, NB > 0.
!     The first function calculated is of order ALPHA, and the
!     last is of order (NB - 1 + ALPHA).
!
! IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
!     and 2 if exponentially scaled K's are to be calculated.
!
! BK    - Working precision output vector of length NB.  If the
!     routine terminates normally (NCALC=NB), the vector BK
!     contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
!     or the corresponding exponentially scaled functions.
!     If (0 < NCALC < NB), BK(I) contains correct function
!     values for I <= NCALC, and contains the ratios
!     K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
! NCALC - Integer output variable indicating possible errors.
!     Before using the vector BK, the user should check that
!     NCALC=NB, i.e., all orders have been calculated to
!     the desired accuracy.  See error returns below.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   EPS    = The smallest positive floating-point number such that
!    1.0+EPS > 1.0
!   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution
!    to equation:
!       W(X) * (1-1/8X+9/128X**2) = beta**minexp
!    where  W(X) = EXP(-X)*SQRT(PI/2X)
!   SQXMIN = Square root of beta**minexp
!   XINF   = Largest positive machine number; approximately
!    beta**maxexp
!   XMIN   = Smallest positive machine number; approximately
!    beta**minexp
!
! Error returns
!
!  In case of an error, NCALC /= NB, and not all K's are
!  calculated to the desired accuracy.
!
!  NCALC < -1:  An argument is out of range. For example,
!   NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >=
!   XMAX.  In this case, the B-vector is not calculated,
!   and NCALC is set to MIN0(NB,0)-2  so that NCALC /= NB.
!  NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
!   K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.  In this case,
!   the B-vector is not calculated.  Note that again
!   NCALC /= NB.
!
!  0 < NCALC < NB: Not all requested function values could
!   be calculated accurately.  BK(I) contains correct function
!   values for I <= NCALC, and contains the ratios
!   K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
  implicit none

  integer ( kind = 4 ) i,iend,itemp,ize,j,k,m,mplus1,nb,ncalc
  real ( kind = 8 ) &
      a,alpha,blpha,bk,bk1,bk2,c,d,dm,d1,d2,d3,enu,eps,estf,estm, &
      ex,four,f0,f1,f2,half,one,p,p0,q,q0,r,ratio,s,sqxmin,t,tinyx, &
      two,twonu,twox,t1,t2,wminf,x,xinf,xmax,xmin,x2by4,zero
  dimension bk(nb),p(8),q(7),r(5),s(4),t(6),estm(6),estf(7)
!
!  Mathematical constants
!    A = LOG(2.D0) - Euler's constant
!    D = SQRT(2.D0/PI)
!
  data half,one,two,zero/0.5d0,1.0d0,2.0d0,0.0d0/
  data four,tinyx/4.0d0,1.0d-10/
  data a/ 0.11593151565841244881d0/,d/0.797884560802865364d0/
!
!  Machine dependent parameters
!
  data eps/2.22d-16/,sqxmin/1.49d-154/,xinf/1.79d+308/
  data xmin/2.23d-308/,xmax/705.342d0/
!
!  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
!                 + Euler's constant
!     Coefficients converted from hex to decimal and modified
!     by W. J. Cody, 2/26/82
!  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
!  T    - Approximation for SINH(Y)/Y
!
  data p/ 0.805629875690432845d00,    0.204045500205365151d02, &
      0.157705605106676174d03,    0.536671116469207504d03, &
      0.900382759291288778d03,    0.730923886650660393d03, &
      0.229299301509425145d03,    0.822467033424113231d00/
  data q/ 0.294601986247850434d02,    0.277577868510221208d03, &
      0.120670325591027438d04,    0.276291444159791519d04, &
      0.344374050506564618d04,    0.221063190113378647d04, &
      0.572267338359892221d03/
  data r/-0.48672575865218401848d+0,  0.13079485869097804016d+2, &
     -0.10196490580880537526d+3,  0.34765409106507813131d+3, &
      0.34958981245219347820d-3/
  data s/-0.25579105509976461286d+2,  0.21257260432226544008d+3, &
     -0.61069018684944109624d+3,  0.42269668805777760407d+3/
  data t/ 0.16125990452916363814d-9, 0.25051878502858255354d-7, &
      0.27557319615147964774d-5, 0.19841269840928373686d-3, &
      0.83333333333334751799d-2, 0.16666666666666666446d+0/
  data estm/5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2, &
        9.3715d0/
  data estf/4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0, &
        8.45096d1, 2.0d1/

  ex = x
  enu = alpha
  ncalc = min ( nb, 0 ) - 2

  if ( &
    ( 0 < nb ) .and. &
    ((enu >= zero) .and. (enu < one)) .and. &
    ((ize >= 1) .and. (ize <= 2)) .and. &
    ((ize /= 1) .or. (ex <= xmax)) .and. &
    (ex > zero)) then

    k = 0
    if (enu < sqxmin) enu = zero

    if (enu > half) then
      k = 1
      enu = enu - one
    end if

    twonu = enu+enu
    iend = nb+k-1
    c = enu*enu
    d3 = -c

    if (ex <= one) then
!
!  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
!  Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
!
      d1 = zero
      d2 = p(1)
      t1 = one
      t2 = q(1)

      do i = 2,7,2
        d1 = c*d1+p(i)
        d2 = c*d2+p(i+1)
        t1 = c*t1+q(i)
        t2 = c*t2+q(i+1)
      end do

      d1 = enu*d1
      t1 = enu*t1
      f1 = log(ex)
      f0 = a+enu*(p(8)-enu*(d1+d2)/(t1+t2))-f1
      q0 = exp(-enu*(a-enu*(p(8)+enu*(d1-d2)/(t1-t2))-f1))
      f1 = enu*f0
      p0 = exp(f1)
!
!  Calculation of F0 =
!
      d1 = r(5)
      t1 = one

      do i = 1,4
        d1 = c*d1+r(i)
        t1 = c*t1+s(i)
      end do

      if (abs(f1) <= half) then

        f1 = f1*f1
        d2 = zero
        do i = 1,6
          d2 = f1*d2+t(i)
        end do
        d2 = f0+f0*f1*d2

      else

        d2 = sinh(f1)/enu

      end if

      f0 = d2-enu*d1/(t1*p0)

      if (ex <= tinyx) then
!
!  X<=1.0E-10
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
!
        bk(1) = f0+ex*f0
        if (ize == 1) bk(1) = bk(1)-ex*bk(1)
        ratio = p0/f0
        c = ex*xinf

        if (k /= 0) then
!
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
!  ALPHA >= 1/2
!
          ncalc = -1

          if (bk(1) >= c/ratio) then
            return
          end if

          bk(1) = ratio*bk(1)/ex
          twonu = twonu+two
          ratio = twonu

        end if

        ncalc = 1

        if (nb == 1) then
          return
        end if
!
!  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
!
        ncalc = -1

        do i = 2,nb

          if (ratio >= c) then
            return
          end if

          bk(i) = ratio/ex
          twonu = twonu+two
          ratio = twonu

        end do

        ncalc = 1
        j = ncalc+1

        do i = j,nb
          if (bk(ncalc) >= xinf/bk(i)) then
            return
          end if
          bk(i) = bk(ncalc)*bk(i)
          ncalc = i
        end do

        return

      else
!
!  1.0E-10 < X <= 1.0
!
        c = one
        x2by4 = ex*ex/four
        p0 = half*p0
        q0 = half*q0
        d1 = -one
        d2 = zero
        bk1 = zero
        bk2 = zero
        f1 = f0
        f2 = p0

  100   continue

        d1 = d1+two
        d2 = d2+one
        d3 = d1+d3
        c = x2by4*c/d2
        f0 = (d2*f0+p0+q0)/d3
        p0 = p0/(d2-enu)
        q0 = q0/(d2+enu)
        t1 = c*f0
        t2 = c*(p0-d2*f0)
        bk1 = bk1+t1
        bk2 = bk2+t2

        if ((abs(t1/(f1+bk1)) > eps) .or. &
           (abs(t2/(f2+bk2)) > eps))  go to 100

        bk1 = f1+bk1
        bk2 = two*(f2+bk2)/ex

        if (ize == 2) then
          d1 = exp(ex)
          bk1 = bk1*d1
          bk2 = bk2*d1
        end if

        wminf = estf(1)*ex+estf(2)

      end if

    else if (eps*ex > one) then
!
!  X > ONE/EPS
!
      ncalc = nb
      bk1 = one / (d*sqrt(ex))
      do i = 1, nb
        bk(i) = bk1
      end do

      return

    else
!
!  X > 1.0
!
      twox = ex+ex
      blpha = zero
      ratio = zero

      if (ex <= four) then
!
!  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
!
        d2 = aint(estm(1)/ex+estm(2))
        m = int(d2)
        d1 = d2+d2
        d2 = d2-half
        d2 = d2*d2
        do i = 2,m
          d1 = d1-two
          d2 = d2-d1
          ratio = (d3+d2)/(twox+d1-ratio)
        end do
!
!  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
!  recurrence and K(ALPHA,X) from the wronskian
!
        d2 = aint(estm(3)*ex+estm(4))
        m = int(d2)
        c = abs(enu)
        d3 = c+c
        d1 = d3-one
        f1 = xmin
        f0 = (two*(c+d2)/ex+half*ex/(c+d2+one))*xmin
        do i = 3,m
          d2 = d2-one
          f2 = (d3+d2+d2)*f0
          blpha = (one+d1/d2)*(f2+blpha)
          f2 = f2/ex+f1
          f1 = f0
          f0 = f2
        end do
        f1 = (d3+two)*f0/ex+f1
        d1 = zero
        t1 = one
        do i = 1,7
          d1 = c*d1+p(i)
          t1 = c*t1+q(i)
        end do
        p0 = exp(c*(a+c*(p(8)-c*d1/t1)-log(ex)))/ex
        f2 = (c+half-ratio)*f1/ex
        bk1 = p0+(d3*f0-f2+f0+blpha)/(f2+f1+f0)*p0
        if (ize == 1) bk1 = bk1*exp(-ex)
        wminf = estf(3)*ex+estf(4)

      else
!
!  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
!  recurrence, for  X > 4.0
!
        dm = aint(estm(5)/ex+estm(6))
        m = int(dm)
        d2 = dm-half
        d2 = d2*d2
        d1 = dm+dm

        do i = 2,m
          dm = dm-one
          d1 = d1-two
          d2 = d2-d1
          ratio = (d3+d2)/(twox+d1-ratio)
          blpha = (ratio+ratio*blpha)/dm
        end do

        bk1 = one/((d+d*blpha)*sqrt(ex))
        if (ize == 1) bk1 = bk1*exp(-ex)
        wminf = estf(5)*(ex-abs(ex-estf(7)))+estf(6)

      end if
!
!  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
!  K(ALPHA+1,X)/K(ALPHA,X)
!
      bk2 = bk1+bk1*(enu+half-ratio)/ex

    end if
!
!  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
!  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
!
    ncalc = nb
    bk(1) = bk1

    if (iend == 0) then
      return
    end if

    j = 2-k
    if (j > 0) bk(j) = bk2
 
    if (iend == 1) then
      return
    end if

    m = min(int(wminf-enu),iend)

    do i = 2,m

      t1 = bk1
      bk1 = bk2
      twonu = twonu+two

      if (ex < one) then
        if (bk1 >= (xinf/twonu)*ex) then
          exit
        end if
      else
        if (bk1/ex >= xinf/twonu) then
          exit
        end if
      end if

      bk2 = twonu/ex*bk1+t1
      itemp = i
      j = j+1
      if (j > 0) bk(j) = bk2

    end do

    m = itemp

    if (m == iend) then
      return
    end if

    ratio = bk2/bk1
    mplus1 = m+1
    ncalc = -1

    do i = mplus1,iend

      twonu = twonu+two
      ratio = twonu/ex+one/ratio
      j = j+1

      if (j > 1) then
        bk(j) = ratio
      else
        if (bk2 >= xinf/ratio) then
          return
        end if
        bk2 = ratio*bk2
      end if

    end do

    ncalc = max(mplus1-k,1)
    if (ncalc == 1) bk(1) = bk2

    if (nb == 1) then
      return
    end if

    j = ncalc+1
    do i = j,nb
      if (bk(ncalc) >= xinf/bk(i)) then
        return
      end if
      bk(i) = bk(ncalc)*bk(i)
      ncalc = i
    end do

  end if

  return
end
subroutine rybesl ( x, alpha, nb, by, ncalc )

!*****************************************************************************80
!
!! RYBESL evaluates a sequence of Bessel Y functions.
!
!  Discussion:
!
!    This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
! X     - Working precision positive real argument for which
!     Y's are to be calculated.
!
! ALPHA - Working precision fractional part of order for which
!     Y's are to be calculated.  0 <= ALPHA < 1.0.
!
! NB    - Integer number of functions to be calculated, NB > 0.
!     The first function calculated is of order ALPHA, and the
!     last is of order (NB - 1 + ALPHA).
!
! BY    - Working precision output vector of length NB.  If the
!     routine terminates normally (NCALC=NB), the vector BY
!     contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
!     If (0 < NCALC < NB), BY(I) contains correct function
!     values for I <= NCALC, and contains the ratios
!     Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
!
! NCALC - Integer output variable indicating possible errors.
!     Before using the vector BY, the user should check that
!     NCALC=NB, i.e., all orders have been calculated to
!     the desired accuracy.  See error returns below.
!
! Explanation of machine-dependent constants.  Let
!
!   beta   = Radix for the floating-point system
!   p  = Number of significant base-beta digits in the
!    significand of a floating-point number
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   EPS    = beta ** (-p)
!   DEL    = Machine number below which sin(x)/x = 1; approximately
!    SQRT(EPS).
!   XMIN   = Smallest acceptable argument for RBESY; approximately
!    max(2*beta**minexp,2/XINF), rounded up
!   XINF   = Largest positive machine number; approximately
!    beta**maxexp
!   THRESH = Lower bound for use of the asymptotic form; approximately
!    AINT(-LOG10(EPS/2.0))+1.0
!   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
!    and cosine functions have lost about half of their
!    precision at that point.
!
! Error returns
!
!  In case of an error, NCALC /= NB, and not all Y's are
!  calculated to the desired accuracy.
!
!  NCALC <= -1:  An argument is out of range. For example,
!   NB <= 0, or ABS(X) >= XLARGE.  In this case,
!   BY(1) = 0.0, the remainder of the BY-vector is not
!   calculated, and NCALC is set to MIN0(NB,0)-1  so that
!   NCALC /= NB.
!  1 < NCALC < NB: Not all requested function values could
!   be calculated accurately.  BY(I) contains correct function
!   values for I <= NCALC, and and the remaining NB-NCALC
!   array elements contain 0.0.
!
! Acknowledgement
!
!  This program draws heavily on Temme's Algol program for Y(a,x)
!  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
!  scheme is used for  x < THRESH, and Campbell's scheme is used
!  in the asymptotic region.  Segments of code from both sources
!  have been translated into Fortran 77, merged, and heavily modified.
!  Modifications include parameterization of machine dependencies,
!  use of a new approximation for ln(gamma(x)), and built-in
!  protection against overflow and destructive underflow.
!
!  Reference: 
!
!    "Bessel functions J_nu(x) and Y_nu(x) of real
!      order and real argument," Campbell, J. B.,
!      Comp. Phy. Comm. 18, 1979, pp. 133-142.
!
!     "On the numerical evaluation of the ordinary
!      Bessel function of the second kind," Temme,
!      N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
!
  implicit none

  integer ( kind = 4 ) i,k,na,nb,ncalc
  real ( kind = 8 ) &
    alfa,alpha,aye,b,by,c,ch,cosmu,d,del,den,ddiv,div,dmu,d1,d2, &
    e,eight,en,enu,en1,eps,even,ex,f,fivpi,g,gamma,h,half,odd, &
    onbpi,one,one5,p,pa,pa1,pi,piby2,pim5,q,qa,qa1,q0,r,s,sinmu, &
    sq2bpi,ten9,term,three,thresh,two,twobyx,x,xinf,xlarge,xmin, &
    xna,x2,ya,ya1,zero
  dimension by(nb),ch(21)
!
!  Mathematical constants
!    FIVPI = 5*PI
!    PIM5 = 5*PI - 15
!    ONBPI = 1/PI
!    PIBY2 = PI/2
!    SQ2BPI = SQUARE ROOT OF 2/PI
!
  data zero,half,one,two,three/0.0d0,0.5d0,1.0d0,2.0d0,3.0d0/
  data eight,one5,ten9/8.0d0,15.0d0,1.9d1/
  data fivpi,piby2/1.5707963267948966192d1,1.5707963267948966192d0/
  data pi,sq2bpi/3.1415926535897932385d0,7.9788456080286535588d-1/
  data pim5,onbpi/7.0796326794896619231d-1,3.1830988618379067154d-1/
!
!  Machine-dependent constants
!
  data del,xmin,xinf,eps/1.0d-8,4.46d-308,1.79d308,1.11d-16/
  data thresh,xlarge/16.0d0,1.0d8/
!
!  Coefficients for Chebyshev polynomial expansion of
!     1/gamma(1-x), abs(x) <= .5
!
  data ch/-0.67735241822398840964d-23,-0.61455180116049879894d-22, &
       0.29017595056104745456d-20, 0.13639417919073099464d-18, &
       0.23826220476859635824d-17,-0.90642907957550702534d-17, &
      -0.14943667065169001769d-14,-0.33919078305362211264d-13, &
      -0.17023776642512729175d-12, 0.91609750938768647911d-11, &
       0.24230957900482704055d-09, 0.17451364971382984243d-08, &
      -0.33126119768180852711d-07,-0.86592079961391259661d-06, &
      -0.49717367041957398581d-05, 0.76309597585908126618d-04, &
       0.12719271366545622927d-02, 0.17063050710955562222d-02, &
      -0.76852840844786673690d-01,-0.28387654227602353814d+00, &
       0.92187029365045265648d+00/

  ex = x
  enu = alpha

  if ( &
    (nb > 0) .and. &
    (x >= xmin) .and. &
    (ex < xlarge) .and. &
    (enu >= zero) .and. &
    (enu < one))  then

    xna = aint(enu+half)
    na = int(xna)
    if (na == 1) enu = enu - xna

    if (enu == -half) then

      p = sq2bpi/sqrt(ex)
      ya = p * sin(ex)
      ya1 = -p * cos(ex)

    else if (ex < three) then
!
!  Use Temme's scheme for small X
!
      b = ex * half
      d = -log(b)
      f = enu * d
      e = b**(-enu)

      if (abs(enu) < del) then
        c = onbpi
      else
        c = enu / sin(enu*pi)
      end if
!
!  Computation of sinh(f)/f
!
      if (abs(f) < one) then
        x2 = f*f
        en = ten9
        s = one
        do i = 1, 9
          s = s*x2/en/(en-one)+one
          en = en - two
        end do
      else
        s = (e - one/e) * half / f
      end if
!
!  Computation of 1/gamma(1-a) using Chebyshev polynomials
!
      x2 = enu*enu*eight
      aye = ch(1)
      even = zero
      alfa = ch(2)
      odd = zero

      do i = 3, 19, 2
        even = -(aye+aye+even)
        aye = -even*x2 - aye + ch(i)
        odd = -(alfa+alfa+odd)
        alfa = -odd*x2 - alfa + ch(i+1)
      end do

      even = (even*half+aye)*x2 - aye + ch(21)
      odd = (odd+alfa)*two
      gamma = odd*enu + even
!
!  End of computation of 1/gamma(1-a)
!
      g = e * gamma
      e = (e + one/e) * half
      f = two*c*(odd*e+even*s*d)
      e = enu*enu
      p = g*c
      q = onbpi / g
      c = enu*piby2

      if (abs(c) < del) then
        r = one
      else
        r = sin(c)/c
      end if

      r = pi*c*r*r
      c = one
      d = - b*b
      h = zero
      ya = f + r*q
      ya1 = p
      en = zero

100   continue

      en = en + one

      if (abs(g/(one+abs(ya))) + abs(h/(one+abs(ya1))) > eps) then
        f = (f*en+p+q)/(en*en-e)
        c = c * d/en
        p = p/(en-enu)
        q = q/(en+enu)
        g = c*(f+r*q)
        h = c*p - en*g
        ya = ya + g
        ya1 = ya1+h
        go to 100
      end if

      ya = -ya
      ya1 = -ya1/b

    else if (ex < thresh) then
!
!  Use Temme's scheme for moderate X
!
      c = (half-enu)*(half+enu)
      b = ex + ex
      e = (ex*onbpi*cos(enu*pi)/eps)
      e = e*e
      p = one
      q = -ex
      r = one + ex*ex
      s = r
      en = two

      do while (r*en*en < e)
        en1 = en+one
        d = (en-one+c/en)/s
        p = (en+en-p*d)/en1
        q = (-b+q*d)/en1
        s = p*p + q*q
        r = r*s
        en = en1
      end do

      f = p/s
      p = f
      g = -q/s
      q = g

220   continue

      en = en - one

      if (en > zero) then
        r = en1*(two-p)-two
        s = b + en1*q
        d = (en-one+c/en)/(r*r+s*s)
        p = d*r
        q = d*s
        e = f + one
        f = p*e - g*q
        g = q*e + p*g
        en1 = en
        go to 220
      end if

      f = one + f
      d = f*f + g*g
      pa = f/d
      qa = -g/d
      d = enu + half -p
      q = q + ex
      pa1 = (pa*q-qa*d)/ex
      qa1 = (qa*q+pa*d)/ex
      b = ex - piby2*(enu+half)
      c = cos(b)
      s = sin(b)
      d = sq2bpi/sqrt(ex)
      ya = d*(pa*s+qa*c)
      ya1 = d*(qa1*s-pa1*c)

    else
!
!  Use Campbell's asymptotic scheme.
!
      na = 0
      d1 = aint(ex/fivpi)
      i = int(d1)
      dmu = ((ex-one5*d1)-d1*pim5)-(alpha+half)*piby2

      if (i-2*(i/2) == 0) then
        cosmu = cos(dmu)
        sinmu = sin(dmu)
         else
        cosmu = -cos(dmu)
        sinmu = -sin(dmu)
      end if

      ddiv = eight * ex
      dmu = alpha
      den = sqrt(ex)

      do k = 1, 2

        p = cosmu
        cosmu = sinmu
        sinmu = -p
        d1 = (two*dmu-one)*(two*dmu+one)
        d2 = zero
        div = ddiv
        p = zero
        q = zero
        q0 = d1/div
        term = q0

        do i = 2, 20
          d2 = d2 + eight
          d1 = d1 - d2
          div = div + ddiv
          term = -term*d1/div
          p = p + term
          d2 = d2 + eight
          d1 = d1 - d2
          div = div + ddiv
          term = term*d1/div
          q = q + term

          if (abs(term) <= eps) then
            exit
          end if

        end do

        p = p + one
        q = q + q0
        if (k == 1) then
          ya = sq2bpi * (p*cosmu-q*sinmu) / den
        else
          ya1 = sq2bpi * (p*cosmu-q*sinmu) / den
        end if

        dmu = dmu + one

      end do

    end if

    if (na == 1) then

      h = two*(enu+one)/ex

      if (h > one) then
        if (abs(ya1) > xinf/h) then
          h = zero
          ya = zero
        end if
      end if

      h = h*ya1 - ya
      ya = ya1
      ya1 = h

    end if
!
!  Now have first one or two Y's
!
    by(1) = ya
    by(2) = ya1

    if (ya1 == zero) then

      ncalc = 1

    else

      aye = one + alpha
      twobyx = two/ex
      ncalc = 2

      do i = 3, nb

        if (twobyx < one) then
          if (abs(by(i-1))*twobyx >= xinf/aye) then
            exit
          end if
        else
          if (abs(by(i-1)) >= xinf/aye/twobyx ) then
            exit
          end if
        end if

        by(i) = twobyx*aye*by(i-1) - by(i-2)
        aye = aye + one
        ncalc = ncalc + 1

      end do

    end if

    do i = ncalc + 1, nb
      by(i) = zero
    end do

  else

    by(1) = zero
    ncalc = min ( nb, 0 ) - 1

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
