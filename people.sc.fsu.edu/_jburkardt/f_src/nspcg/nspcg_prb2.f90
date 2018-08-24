program main

!*******************************************************************************
!
!! MAIN is the main program for NSPCG_PRB2.
!
!  Discussion:
!
!    NSPCG_PRB2 runs the NSPCG tests.
!
!  Modified:
!
!    08 August 2006
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSPCG_PRB2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tests for NSPCG.'

  call test01
  call test02
  call test03
  call test04
  call test05

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSPCG_PRB2'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*******************************************************************************
!
!! TEST01 tests NSPCG with MIC2 and CG.
!
!  Modified:
!
!    08 August 2006
!
  implicit none

  integer, parameter :: inw = 300
  integer, parameter :: mdim = 4
  integer, parameter :: nw = 600
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10

  integer, parameter :: n = nx * ny
  integer, parameter :: ndim = n

  external cg
  real ( kind = 8 ) coef(ndim,mdim)
  real ( kind = 8 ) h
  integer i
  integer ier
  integer ip(1)
  integer iparm(30)
  integer isym
  integer iwksp(inw)
  integer j
  integer jcoef(5)
  integer k
  integer maxnz
  external mic2
  integer p(1)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(1)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) wksp(nw)
  real ( kind = 8 ) xx
  real ( kind = 8 ) yy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Set up and solve Laplace''s equation'
  write ( *, '(a)' ) '  inside the unit square, using a mesh'
  write ( *, '(a,i8,a,i8)' ) '  of NX = ', nx, ' by NY = ', ny
  write ( *, '(a)' ) '  strictly interior nodes.'

  maxnz = 3
  h = 1.0D+00 / 11.0D+00
!
!  Generate Laplace's equation.
!
  isym = 0
  call matgen ( n, nx, ny, ndim, maxnz, jcoef, coef, rhs, h, isym )
!
!  Initialize the parameters to their default values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize NSPCG parameters to default values.'

  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(2) = 50
  iparm(3) = 3

  rparm(1) = 1.0D-04
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call NSPCG to solve the system.'
  write ( *, '(a)' ) '  We use the CG conjugate gradient module.'

  call nspcg ( mic2, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       I       J      X         Y        U(exact)      U(computed)'
  write ( *, '(a)' ) ' '

  k = 0

  do j = 1, ny

    yy = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      xx = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )
      k = k + 1
      uexact = 1.0D+00 + xx * yy

      write ( *, '(2x,i6,2x,i6,2x,f8.4,2x,f8.4,2x,g14.6,2x,g14.6)' ) &
        i, j, xx, yy, uexact, u(k)
  
    end do
  end do

  return
end
subroutine matgen ( n, nx, ny, ndim, maxnz, jcoef, coef, rhs, h, isym )

!*******************************************************************************
!
!! MATGEN generates the data for a discretization of a PDE.
!
!  Discussion:
!
!    MATGEN generates the itpack matrices jcoef, coef, and rhs
!    for a 5-point central difference approximation of the
!    self-adjoint pde
!
!      a * uxx + c * uyy + f * u = g
!
!    where the coefficients a, c, f, and g are functions of (x,y)
!    (supplied by the user) and the domain is the unit square
!    (0,1) x (0,1).  Dirichlet boundary conditions are imposed
!    upon the boundary in the form of the function ub(x,y) also
!    supplied by the user.
!
!  Parameters:
!
!    Output, integer N, the number of linear equations.
!
!    Input, integer NX, NY, the number of nodes in the X and Y directions.
!
!    Input, integer NDIM, the row dimension of JCOEF and COEF.
!
!    Output, integer MAXNZ, the maximum number of nonzeros per row.
!
!    Output, integer JCOEF(5), defines the NSPCG data structure.
!
!    Output, real ( kind = 8 ) COEF(NDIM,5), the equation coefficients.
!
!    Output, real ( kind = 8 ) RHS(N), the right-hand-side values.
!
!    Input, real ( kind = 8 ) H, the mesh spacing.
!
!    Input, integer ISYM, 0 if a symmetric matrix is to be generated.
!
  implicit none

  integer n
  integer ndim

  real ( kind = 8 ) a
  real ( kind = 8 ) ae
  real ( kind = 8 ) aw
  real ( kind = 8 ) c
  real ( kind = 8 ) cc
  real ( kind = 8 ) cn
  real ( kind = 8 ) coef(ndim,5)
  real ( kind = 8 ) cs
  real ( kind = 8 ) f
  real ( kind = 8 ) fp
  real ( kind = 8 ) g
  real ( kind = 8 ) gp
  real ( kind = 8 ) h
  integer i
  integer isym
  integer j
  integer jcoef(5)
  integer maxnz
  integer neq
  integer nx
  integer ny
  real ( kind = 8 ) rhs(n)
  logical symm
  real ( kind = 8 ) ub
  real ( kind = 8 ) x
  real ( kind = 8 ) xx
  real ( kind = 8 ) y
  real ( kind = 8 ) yy
!
!  Statement functions
!
  a(x,y) = 1.0D+00
  c(x,y) = 2.0D+00
  f(x,y) = 0.0D+00
  g(x,y) = 0.0D+00
  ub(x,y) = 1.0D+00 + x * y

  symm = ( isym == 0 )
  maxnz = 5

  if ( symm ) then
    maxnz = 3
  else
    maxnz = 5
  end if
!
!  Loop on equations, from left to right and from down to up.
!
  neq = 0

  do j = 1, ny

    yy = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      xx = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )

      neq = neq + 1
      ae = a ( xx + 0.5D+00 * h, yy )
      aw = a ( xx - 0.5D+00 * h, yy )
      cn = c ( xx, yy + 0.5D+00 * h )
      cs = c ( xx, yy - 0.5D+00 * h )
      fp = f ( xx, yy )
      gp = g ( xx, yy )

      cc = ae + cn + aw + cs - fp * h**2
!
!  Center point
!
      coef(neq,1) = cc
      rhs(neq) = - h**2 * gp
!
!  East point.
!
      if ( i /= nx ) then
        coef(neq,2) = -ae
      else
        coef(neq,2) = 0.0D+00
        rhs(neq) = rhs(neq) + ae * ub ( 1.0D+00, yy )
      end if
!
!  North point.
!
      if ( j /= ny ) then
        coef(neq,3) = -cn
      else
        coef(neq,3) = 0.0D+00
        rhs(neq) = rhs(neq) + cn * ub ( xx, 1.0D+00 )
      end if
!
!  West point.
!
      if ( i /= 1 ) then
        if ( .not. symm ) then
          coef(neq,4) = -aw
        end if
      else
        if ( .not. symm ) then
          coef(neq,4) = 0.0D+00
        end if
        rhs(neq) = rhs(neq) + aw * ub ( 0.0D+00, yy )
      end if
!
!  South point.
!
      if ( j /= 1 ) then
        if (.not. symm ) then
          coef(neq,5) = -cs
        end if
      else
        if ( .not. symm ) then
          coef(neq,5) = 0.0D+00
        end if
        rhs(neq) = rhs(neq) + cs * ub ( xx, 0.0D+00 )
      end if

    end do
  end do
!
!  Define data structure 2.
!
  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx

  return
end
subroutine test02

!*******************************************************************************
!
!! TEST02 ...
!
  implicit none

  integer, parameter :: mdim = 3
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10

  integer, parameter :: n = nx * ny

  integer, parameter :: ndim = n
  integer, parameter :: nw = 10 * n
  integer, parameter :: inw = 3 * n

  external bic2
  external cg
  real ( kind = 8 ) coef(ndim,mdim)
  integer i
  integer ier
  integer ip(1)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(3)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer p(1)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(1)
  real ( kind = 8 ) wksp(nw)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'

  maxnz = 3
!
!  Set COEF.
!
  coef(1:n,1) = 6.0D+00
  coef(1:n,2) = -1.0D+00
  coef(1:n,3) = -2.0D+00
!
!  Set the right hand side.
!
  rhs(1:n) = 0.0D+00
!
!  Make modifications for boundary values.
!
  k = 0

  do j = 1, ny

    y = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      x = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0D+00
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0D+00 * ( 1.0D+00 + x )
        coef(k,3) = 0.0D+00
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0D+00
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0D+00 + y
        coef(k,2) = 0.0D+00
      end if

    end do
  end do

  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(4) = 6
  iparm(18) = 1
  iparm(19) = nx

  rparm(1) = 1.0D-04
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Solve the system.
!
  call nspcg ( bic2, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test03

!*******************************************************************************
!
!! TEST03 ...
!
  implicit none

  integer, parameter :: mdim = 5
  integer, parameter :: inw = 500
  integer, parameter :: nw = 800
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10

  integer, parameter :: n = nx * ny
  integer, parameter :: ndim = n

  external cg
  real ( kind = 8 ) coef(ndim,mdim)
  integer i
  integer ier
  integer ip(n)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(ndim,5)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer p(n)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  external rs6
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(1)
  real ( kind = 8 ) wksp(nw)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
!
!  generate coef, jcoef, and rhs.
!
  maxnz = 5

  do i = 1, n

    coef(i,1) = 6.0D+00
    coef(i,2) = -1.0D+00
    coef(i,3) = -2.0D+00
    coef(i,4) = -1.0D+00
    coef(i,5) = -2.0D+00

    jcoef(i,1) = i
    jcoef(i,2) = i+1
    jcoef(i,3) = i+nx
    jcoef(i,4) = i-1
    jcoef(i,5) = i-nx

  end do

  rhs(1:n) = 0.0D+00

  k = 0

  do j = 1, ny

    y = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      x = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0D+00
        coef(k,5) = 0.0D+00
        jcoef(k,5) = 0
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0D+00 * ( 1.0D+00 + x )
        coef(k,3) = 0.0D+00
        jcoef(k,3) = 0
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0D+00
        coef(k,4) = 0.0D+00
        jcoef(k,4) = 0
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0D+00 + y
        coef(k,2) = 0.0D+00
        jcoef(k,2) = 0
      end if

    end do
  end do
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(12) = 1
  iparm(14) = 1
  iparm(23) = 0

  rparm(1) = 1.0D-04
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Call REDBLK to check for property A.
!
  call redblk ( ndim, n, maxnz, coef, jcoef, p, ip, 1, iwksp, ier )
!
!  Solve the system.
!
  call nspcg ( rs6, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, &
    u, ubar, rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test04

!*******************************************************************************
!
!! TEST04 ...
!
  implicit none

  integer, parameter :: inw = 400
  integer, parameter :: mdim = 5
  integer, parameter :: nw = 1000
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10
  integer, parameter :: nz = 1

  integer, parameter :: n = nx * ny * nz
  integer, parameter :: ndim = n

  real ( kind = 8 ) coef(ndim,mdim)
  integer i
  integer ier
  integer ip(n)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(5)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer nxp
  integer nyp
  integer nzp
  integer p(n)
  integer patt(2)
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  external sor
  external sor7
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(1)
  real ( kind = 8 ) wksp(nw)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
!
!  generate coef, jcoef, and rhs.
!
  maxnz = 3

  do i = 1, n
    coef(i,1) = 6.0D+00
    coef(i,2) = -1.0D+00
    coef(i,3) = -2.0D+00
  end do

  rhs(1:n) = 0.0D+00

  k = 0
  do j = 1, ny

    y = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      x = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0D+00
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0D+00 * ( 1.0D+00 + x ) 
        coef(k,3) = 0.0D+00
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0D+00
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0D+00 + y
        coef(k,2) = 0.0D+00
      end if

    end do
  end do

  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(14) = 1
  iparm(19) = nx

  rparm(1) = 1.0D-04
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0
!
!  Define a color pattern.
!
  nxp = 1
  nyp = 2
  nzp = 1

  patt(1) = 1
  patt(2) = 2

  call color ( nxp, nyp, nzp, nx, ny, nz, patt, p )
!
!  Solve the system.
!
  call nspcg ( sor7, sor, ndim, mdim, n, maxnz, coef, jcoef, p, ip, &
    u, ubar, rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test05

!*******************************************************************************
!
!! TEST05 tests SORW.
!
  implicit none

  integer, parameter :: nw = 800
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10

  integer, parameter :: n = nx * ny

  real ( kind = 8 ) coef(1)
  integer i
  integer ier
  integer iparm(30)
  integer j
  integer jcoef(1)
  integer jwfac(1)
  integer k
  integer khi
  integer klo
  external mult
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) rparm(30)
  external sorpass
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ubar(1)
  real ( kind = 8 ) wfac(1)
  real ( kind = 8 ) wksp(nw)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SORW implements the SOR iteration.'
!
!  Initialize the right hand side to zero, but adjust it for boundary values.
!
  rhs(1:n) = 0.0D+00

  k = 0
  do j = 1, ny

    y = real ( j, kind = 8 ) / real ( ny + 1, kind = 8 )

    do i = 1, nx

      x = real ( i, kind = 8 ) / real ( nx + 1, kind = 8 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0D+00
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0D+00 * ( 1.0D+00 + x )
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0D+00
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0D+00 + y
      end if

    end do
  end do
!
!  Initialize the parameters to their default values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize parameters to default values.'

  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  rparm(1) = 1.0D-04
!
!  Set an initial value for the solution.
!
  u(1:n) = 0.0D+00
!
!  Solve the system.
!
  call sorw ( mult, sorpass, coef, jcoef, wfac, jwfac, n, u, ubar, rhs, &
    wksp, nw, iparm, rparm, ier )
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
subroutine mult ( coef, jcoef, wfac, jwfac, n, x, y )

!*******************************************************************************
!
!! MULT ...
!
  implicit none

  integer n

  real ( kind = 8 ) coef
  integer i
  integer jcoef
  integer jwfac
  integer nx
  real ( kind = 8 ) wfac
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  nx = 10
!
!  Compute the product as if first superdiagonal and first
!  subdiagonal were full.
!
  y(1) = 6.0D+00 * x(1) - x(2) - 2.0D+00 * x(nx+1)
  do i = 2, nx
    y(i) = 6.0D+00 * x(i) - x(i+1) - x(i-1) - 2.0D+00 * x(i+nx)
  end do

  do i = nx+1, n-nx
    y(i) = 6.0D+00 * x(i) - x(i+1) - x(i-1) - 2.0D+00 * ( x(i+nx) + x(i-nx) )
  end do

  do i = n-nx+1, n-1
    y(i) = 6.0D+00 * x(i) - x(i+1) - x(i-1) - 2.0D+00 * x(i-nx)
  end do
  y(n) = 6.0D+00 * x(n) - x(n-1) - 2.0D+00 * x(n-nx)
!
!  Make corrections to y vector for zeros in first superdiagonal
!  and first subdiagonal.
!
  do i = nx, n-nx,nx
    y(i) = y(i) + x(i+1)
  end do

  do i = nx+1, n-nx+1, nx
    y(i) = y(i) + x(i-1)
  end do

  return
end
subroutine sorpass ( coef, jcoef, wfac, jwfac, n, u, rhs, unew )

!*******************************************************************************
!
!! SORPASS does an SOR iteration.
!
!  Discussion:
!
!    unew = inv ( (1/w)*d+l) * (((1/w)-1)*d*un+rhs-u*un)
!
!  Parameters:
!
!    n      order of system
!
!    u      current solution vector
!
!    rhs    right hand side
!
!    unew   updated solution vector
!
  implicit none

  real ( kind = 8 ) alphab
  real ( kind = 8 ) betab
  real ( kind = 8 ) coef(1)
  real ( kind = 8 ) con
  real ( kind = 8 ) fff
  integer i
  integer ibgn
  integer iend
  integer j
  integer jcoef(1)
  integer jwfac(1)
  integer n
  integer nx
  real ( kind = 8 ) omega
  logical omgadp
  real ( kind = 8 ) rhs(1)
  real ( kind = 8 ) specr
  real ( kind = 8 ) u(1)
  real ( kind = 8 ) unew(1)
  real ( kind = 8 ) wfac(1)

  common / itcom5 / omgadp
  common / itcom55 /  omega, alphab, betab, fff, specr
!
!  temp=((1/w)-1)*d*un+rhs-u*un
!     unew is used for temp.
!
  nx = 10
  con = 6.0D+00 * ( 1.0D+00 / omega - 1.0D+00 )

  do i = 1, n
    unew(i) = con * u(i) + rhs(i)
  end do

  do i = 1, n-1
    unew(i) = unew(i) + u(i+1)
  end do

  do i = 1, n-nx
    unew(i) = unew(i) + 2.0D+00 * u(i+nx)
  end do

  do i = nx, n-nx, nx
    unew(i) = unew(i) - u(i+1)
  end do
!
!  unew = inv((1/w)*d+l) * temp
!
  con = omega / 6.0D+00

  do j = 1, nx

    ibgn = ( j - 1 ) * nx + 1
    iend = j * nx
    unew(ibgn) = con * unew(ibgn)

    do i = ibgn+1, iend
      unew(i) = con * ( unew(i) + unew(i-1) )
    end do

    if ( j /= nx ) then
      do i = ibgn, iend
        unew(i+nx) = unew(i+nx) + 2.0D+00 * unew(i)
      end do
    end if

  end do

  return
end
