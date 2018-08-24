program main

!*****************************************************************************80
!
!! MAIN is the main program for NSPCG_PRB.
!
!  Discussion:
!
!    NSPCG_PRB tests the NSPCG library.
!
!  Modified:
!
!    25 June 2007
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSPCG_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NSPCG library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSPCG_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NSPCG with MIC1 and CG.
!
!  Discussion:
!
!    This routine is meant to be a simple demonstration of using NSPCG.
!
!    So, in particular, the linear system is one of the simplest, being
!    just the tridiagonal matrix A = tridiag [ -1, 2, -1], with solution
!    vector X = [ 1, 2, 3, ..., N ] and right hand side 
!    B = [ 0, 0, ..., 0, N+1].
!
!    The storage mode chosen is the "primary storage mode", which
!    corresponds to NSTORE = 1.  The default storage mode is NSTORE = 2.
!    In order to request NSTORE = 1, we must reset IPARM(12) to 1.
!
!    Since we are using storage mode 1, the preconditioner that we
!    use must have a last digit of 1 (so it's MIC1, and we cannot
!    use MIC2, MIC3, MIC4 or MIC5!).
!
!    CG is the name of the solver we are going to use.
!    MIC1 is the name of the preconditioner we are going to use.
!
  implicit none

  integer ( kind = 4 ), parameter :: itmax = 100
  integer ( kind = 4 ), parameter :: maxnz = 3
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ), parameter :: inw = 20
  integer ( kind = 4 ), parameter :: nw = 3*n + 2*itmax + n + n + 10

  external cg
  real ( kind = 8 ) coef(n,maxnz)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) iparm(30)
  integer ( kind = 4 ) iwksp(inw)
  integer ( kind = 4 ) jcoef(n,maxnz)
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) maxnzz
  integer ( kind = 4 ) mdim
  external mic1
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) p(n)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(n)
  real ( kind = 8 ) wksp(nw)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate a simple use of NSPCG.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use "primary storage" for the matrix'
  write ( *, '(a)' ) '  (storage mode 1).'
!
!  Because A is a tridiagonal matrix, we can use the "primary storage format",
!  in which A is an N by 3 array of values, and JCOEF is an N by 3 array
!  of column indices.
!
   coef(1,1:3) =   (/ 2.0D+00, -1.0D+00, 0.0D+00 /)
   jcoef(1,1:3) =  (/ 1, 2, 0 /)

   do i = 2, n-1
     coef(i,1:3) = (/ -1.0D+00, 2.0D+00, -1.0D+00 /)
     jcoef(i,1:3) = (/ i-1, i, i+1 /)
   end do

   coef(n,1:3) = (/ -1.0D+00, 2.0D+00, 0.0D+00 /)
   jcoef(n,1:3) = (/ n-1, n, 0 /)
!
!  Set an estimate for the solution.
!
  u(1:n) = 0.0D+00
!
!  Set the exact solution.
!  Obviously, this is not normally known.
!  We supply it because NSPCG can use it to tell you how many
!  iterations it would take to reach the exact solution, it
!  it is known.
!
  do i = 1, n
    ubar(i) = real ( i, kind = 8 )
  end do
!
!  Set the right hand side.
!
  rhs(1:n-1) = 0.0D+00
  rhs(n) = real ( n + 1, kind = 8 )
!
!  Initialize the parameters to their default values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize parameters to default values.'

  call dfault ( iparm, rparm )
!
!  Adjust some parameters.
!
!  Right now, we've set ITMAX = 100, which is the default.  But
!  if we change it, we'll want to let NSPCG know too, by resetting
!  IPARM(2).
!
!  IPARM(12) indicates the storage mode we have chosen.  We have
!  chosen storage mode 1 (which is why, for instance, we are calling
!  MIC1), so we need to let NSPCG know this, since it's not the
!  default choice.)
!
  iparm(2) = itmax
  iparm(12) = 1
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call NSPCG to solve the system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The preconditioner is MIC1;'
  write ( *, '(a)' ) '  The accelerator is CG.'
!
!  Warning:
!
!    Some of the arguments to NSPCG are used as both input and output variables.
!    Thus, it is ILLEGAL to pass MAXNZ directly, since it's a parameter,
!    and NSPCG is going to change it on output.  
!
  mdim = 3
  ndim = n
  maxnzz = maxnz

  call nspcg ( mic1, cg, ndim, mdim, n, maxnzz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) u(klo:khi)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests NSPCG with MIC1 and CG.
!
!  Discussion:
!
!    This routine solves the same linear system as the previous one,
!    but uses a different storage mode, called "nonsymmetric
!    coordinate storage", which is NSTORE = 5.
!
!    Since the default storage mode is NSTORE = 2.
!    In order to request NSTORE = 5, we must reset IPARM(12) to 1.
!
!    Since we are using storage mode 5, the preconditioner that we
!    use must have a last digit of 5.  There is no "MIC5" preconditioner,
!    so we went with a Jacobi preconditioner called "JAC5".  Our
!    choice of accelerator and preconditioner determines the necessary
!    sizes NW and INW of the work arrays.
!
!    Our choice of storage mode also affects the shape and size of the 
!    COEF and JCOEF arrays.
!
!    CG is the name of the solver we are going to use.
!
!    JAC5 is the name of the preconditioner we are going to use.
!
!    Note that JCOEF is declared to have first dimension MAXNZ, but
!    actually only uses a first dimension of NZ (whose value we don't
!    know yet.)  We will have to be careful to give both of these numbers
!    to NPSCG when it is time to call it.
!
  implicit none

  integer ( kind = 4 ), parameter :: itmax = 100
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ), parameter :: maxnz = 3*n
  integer ( kind = 4 ), parameter :: inw = 2*n+1
  integer ( kind = 4 ), parameter :: nw = 3*n + 2*itmax + n + n + 10

  external cg
  real ( kind = 8 ) coef(maxnz)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) iparm(30)
  integer ( kind = 4 ) iwksp(inw)
  external jac5
  integer ( kind = 4 ) jcoef(maxnz,2)
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) maxnzz
  integer ( kind = 4 ) mdim
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) p(n)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(n)
  real ( kind = 8 ) wksp(nw)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Demonstrate a simple use of NSPCG.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use "nonsymmetric coordinate storage" '
  write ( *, '(a)' ) '  for the matrix (storage mode 5).'
!
!  Store A using nonsymmetric coordinate format.
!  NZ will keep track of the next free space in COEF and JCOEF,
!  and at the end of the loop will equal the total number of
!  nonzero entries in the matrix.
!
  nz = 0
   
  do i = 1, n

    if ( 1 < i ) then
      nz = nz + 1
      coef(nz) = -1.0D+00
      jcoef(nz,1) = i
      jcoef(nz,2) = i-1
    end if

    nz = nz + 1
    coef(nz) = 2.0D+00
    jcoef(nz,1) = i
    jcoef(nz,2) = i

    if ( i < n ) then
      nz = nz + 1
      coef(nz) = -1.0D+00
      jcoef(nz,1) = i
      jcoef(nz,2) = i+1
    end if

  end do
!
!  Set an estimate for the solution.
!
  u(1:n) = 0.0D+00
!
!  Set the exact solution.
!  Obviously, this is not normally known.
!  We supply it because NSPCG can use it to tell you how many
!  iterations it would take to reach the exact solution, it
!  it is known.
!
  do i = 1, n
    ubar(i) = real ( i, kind = 8 )
  end do
!
!  Set the right hand side.
!
  rhs(1:n-1) = 0.0D+00
  rhs(n) = real ( n + 1, kind = 8 )
!
!  Initialize the parameters to their default values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize parameters to default values.'

  call dfault ( iparm, rparm )
!
!  Adjust some parameters.
!
!  Right now, we've set ITMAX = 100, which is the default.  But
!  if we change it, we'll want to let NSPCG know too, by resetting
!  IPARM(2).
!
!  IPARM(12) indicates the storage mode we have chosen.  We have
!  chosen storage mode 1 (which is why, for instance, we are calling
!  MIC1), so we need to let NSPCG know this, since it's not the
!  default choice.)
!
  iparm(2) = itmax
  iparm(12) = 5
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call NSPCG to solve the system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The preconditioner is JAC5;'
  write ( *, '(a)' ) '  The accelerator is CG.'
!
!  Note that NDIM is the leading dimension of JCOEF.
!  JCOEF is a two dimensional array, and the first dimension must
!  be at least NZ (the number of nonzeros) but could be more if
!  you overestimate the space you need.
!
  ndim = maxnz
  mdim = 1
  maxnzz = nz

  call nspcg ( jac5, cg, ndim, mdim, n, maxnzz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) u(klo:khi)
  end do

  return
end
