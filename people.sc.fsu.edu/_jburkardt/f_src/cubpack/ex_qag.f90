MODULE Integrand
  PRIVATE
  PUBLIC :: F
  CONTAINS
    FUNCTION F(NUMFUN,X) RESULT(Value)
    USE Precision_Model
    INTEGER, INTENT(IN) :: NUMFUN
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
    REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
    Value(1) = (abs(x(1)-1.0_stnd/3.0_stnd))**(0.8_stnd)
   
    RETURN 
    END FUNCTION F
END MODULE Integrand

PROGRAM Example_QAG

!*****************************************************************************80
!
!! EXAMPLE_QAQ is equivalent to NAG's test for QAG from QUADPACK.
!
!  Discussion:
!
!    This program uses CUBPACK to estimate the integral of a function over
!    an interval.
!
!  Local Parameters:
!
!    Local integer N, the spatial dimension.
!
  USE Precision_Model
  USE CUI                      ! Cubpack User Interface
  USE Integrand

  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 1

  REAL ( kind = stnd ) :: AbsErr
  REAL ( kind = stnd ) :: epsrel
  INTEGER, PARAMETER :: Finite_interval = 1
  REAL ( kind = stnd ) :: IntegralValue
  INTEGER :: j
  INTEGER :: key
  INTEGER :: NEval
  INTEGER :: RgType
  REAL ( kind = stnd ), DIMENSION(1:n,0:n) :: Vertices

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'EX_QAG:'
  write ( *, '(a)' ) '  Use QAQ to estimate the integral of a function'
  write ( *, '(a)' ) '  over a 1D interval.'

  epsrel = 0.1_stnd

  do j = 1, 10

    epsrel = epsrel * 0.1_stnd

    do key = 1, 6

      RgType = Finite_interval
      Vertices(1,:) = (/ 0 , 1 /)

      CALL CUBATR ( n, F, Vertices, RgType, IntegralValue, AbsErr, &
        MAXPTS = 50000, EpsRel = epsrel, Key = key, NEval = NEval, JOB = 1)

      write ( *, '(a)' ) ''
      write(unit=*,fmt="( "" Results for key = "",I3 )") KEY
      write(unit=*,fmt="( "" Results for epsrel = "",es9.2 )") EPSREL
      write(unit=*,fmt="( "" Integral approximation = "",es15.8 )") IntegralValue
      write(unit=*,fmt="( "" Estimate of absolute error = "",es9.2 )") ABSERR
      write(unit=*,fmt="( "" Number of function evaluattions = "",i5 )") NEVAL

    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'EX_QAG:'
  write ( *, '(a)' ) '  Normal end of execution.'

  STOP 
END PROGRAM Example_QAG
                                        
