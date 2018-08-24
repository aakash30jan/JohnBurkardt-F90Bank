program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM2D_POISSON_RECTANGLE_LINEAR.
!
!  Discussion:
!
!    This program solves
!
!      - d2U(X,Y)/dx2 - d2U(X,Y)/dy2 = F(X,Y)
!
!    in a rectangular region in the plane.
!
!    Along the boundary of the region, Dirichlet conditions
!    are imposed:
!
!      U(X,Y) = G(X,Y)
!
!    The code uses continuous piecewise linear basis functions on
!    triangles determined by a uniform grid of NX by NY points.
!
!    u    =      sin ( pi * x ) * sin ( pi * y ) + x
!
!    dudx = pi * cos ( pi * x ) * sin ( pi * y ) + 1
!    dudy = pi * sin ( pi * x ) * cos ( pi * y )
!
!    d2udx2 = - pi * pi * sin ( pi * x ) * sin ( pi * y )
!    d2udy2 = - pi * pi * sin ( pi * x ) * sin ( pi * y )
!
!    rhs  = 2 * pi * pi * sin ( pi * x ) * sin ( pi * y )
!
!    THINGS YOU CAN EASILY CHANGE:
!
!    1) Change NX or NY, the number of nodes in the X and Y directions.
!    2) Change XL, XR, YB, YT, the left, right, bottom and top limits of the 
!       rectangle.
!    3) Change the exact solution in the EXACT routine, but make sure you also
!       modify the formula for RHS in the assembly portion of the program!
!
!    HARDER TO CHANGE:
!
!    4) Change from "linear" to "quadratic" triangles;
!    5) Change the region from a rectangle to a general triangulated region;
!    6) Store the matrix as a sparse matrix so you can solve bigger systems.
!    7) Handle Neumann boundary conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 17
  integer ( kind = 4 ), parameter :: ny = 17

  real ( kind = 8 ), allocatable, dimension (:,:) :: a
  real ( kind = 8 ) area
  real ( kind = 8 ), allocatable, dimension (:) :: b
  real ( kind = 8 ) dqidx
  real ( kind = 8 ) dqidy
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  integer ( kind = 4 ) e
  integer ( kind = 4 ), allocatable, dimension(:,:) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) nq1
  integer ( kind = 4 ) nq2
  integer ( kind = 4 ) nti1
  integer ( kind = 4 ) nti2
  integer ( kind = 4 ) nti3
  integer ( kind = 4 ) ntj1
  integer ( kind = 4 ) ntj2
  integer ( kind = 4 ) ntj3
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) q1
  integer ( kind = 4 ) q2
  real ( kind = 8 ) qi
  real ( kind = 8 ) qj
  integer ( kind = 4 ) ti1
  integer ( kind = 4 ) ti2
  integer ( kind = 4 ) ti3
  integer ( kind = 4 ) tj1
  integer ( kind = 4 ) tj2
  integer ( kind = 4 ) tj3
  real ( kind = 8 ) rhs
  real ( kind = 8 ) u
  real ( kind = 8 ) wq
  real ( kind = 8 ), parameter :: xl = 0.0D+00
  real ( kind = 8 ) xq
  real ( kind = 8 ), parameter :: xr = 1.0D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xy
  real ( kind = 8 ), parameter :: yb = 0.0D+00
  real ( kind = 8 ) yq
  real ( kind = 8 ), parameter :: yt = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_POISSON_RECTANGLE_LINEAR'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution of the Poisson equation:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  - Uxx - Uyy = F(x,y) inside the region,'
  write ( *, '(a)' ) '       U(x,y) = G(x,y) on the boundary of the region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The region is a rectangle, defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6,a,g14.6)' ) xl,' = XL<= X <= XR = ', xr
  write ( *, '(g14.6,a,g14.6)' ) yb,' = YB<= Y <= YT = ', yt
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The finite element method is used, with piecewise'
  write ( *, '(a)' ) '  linear basis functions on 3 node triangular'
  write ( *, '(a)' ) '  elements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The corner nodes of the triangles are generated by an'
  write ( *, '(a)' ) '  underlying grid whose dimensions are'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NX =                       ', nx
  write ( *, '(a,i8)' ) '  NY =                       ', ny
!
!  NODE COORDINATES
!
!  Numbering of nodes is suggested by the following 5x10 example:
!
!    J=5 | K=41  K=42 ... K=50
!    ... |
!    J=2 | K=11  K=12 ... K=20
!    J=1 | K= 1  K= 2     K=10
!        +--------------------
!          I= 1  I= 2 ... I=10
!
  node_num = nx * ny

  write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num

  allocate ( xy(1:2,1:node_num) )

  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1

      xy(1,k) = ( real ( nx - i,     kind = 8 ) * xl   &
                + real (      i - 1, kind = 8 ) * xr ) &
                / real ( nx     - 1, kind = 8 )

      xy(2,k) = ( real ( ny - j,     kind = 8 ) * yb   &
                + real (      j - 1, kind = 8 ) * yt ) &
                / real ( ny     - 1, kind = 8 )

    end do
  end do
!
!  ELEMENT array
!
!  Organize the nodes into a grid of 3-node triangles.
!  Here is part of the diagram for a 5x10 example:
!
!    |  \ |  \ |  \ |
!    |   \|   \|   \|
!   21---22---23---24--
!    |\ 8 |\10 |\12 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |  \ |
!    |  7\|  9\| 11\|   \|
!   11---12---13---14---15---16---17---18---19---20
!    |\ 2 |\ 4 |\ 6 |\  8|                   |\ 18|
!    | \  | \  | \  | \  |                   | \  |
!    |  \ |  \ |  \ |  \ |      ...          |  \ |
!    |  1\|  3\|  5\| 7 \|                   |17 \|
!    1----2----3----4----5----6----7----8----9---10
!
  element_num = 2 * ( nx - 1 ) * ( ny - 1 )

  write ( *, '(a,i8)' ) '  Number of elements =       ', element_num

  allocate ( element_node(1:3,1:element_num) )

  k = 0

  do j = 1, ny - 1
    do i = 1, nx - 1

      k = k + 1
      element_node(1,k) = i     + ( j - 1 ) * nx
      element_node(2,k) = i + 1 + ( j - 1 ) * nx
      element_node(3,k) = i     +   j       * nx

      k = k + 1
      element_node(1,k) = i + 1 +   j       * nx
      element_node(2,k) = i     +   j       * nx
      element_node(3,k) = i + 1 + ( j - 1 ) * nx

    end do
  end do
!
!  ASSEMBLE THE SYSTEM
!
!  Assemble the coefficient matrix A and the right-hand side B of the
!  finite element equations, ignoring boundary conditions.
!
  allocate ( b(1:node_num) )
  allocate ( a(1:node_num,1:node_num) )

  b(1:node_num) = 0.0D+00
  a(1:node_num,1:node_num) = 0.0D+00

  do e = 1, element_num

    i1 = element_node(1,e)
    i2 = element_node(2,e)
    i3 = element_node(3,e)

    area = 0.5D+00 * &
      ( xy(1,i1) * ( xy(2,i2) - xy(2,i3) ) &
      + xy(1,i2) * ( xy(2,i3) - xy(2,i1) ) &
      + xy(1,i3) * ( xy(2,i1) - xy(2,i2) ) )
!
!  Consider each quadrature point.
!  Here, we use the midside nodes as quadrature points.
!
    do q1 = 1, 3

      q2 = mod ( q1, 3 ) + 1

      nq1 = element_node(q1,e)
      nq2 = element_node(q2,e)

      xq = 0.5D+00 * ( xy(1,nq1) + xy(1,nq2) )
      yq = 0.5D+00 * ( xy(2,nq1) + xy(2,nq2) )
      wq = 1.0D+00 / 3.0D+00
!
!  Consider each test function in the element.
!
      do ti1 = 1, 3

        ti2 = mod ( ti1,     3 ) + 1
        ti3 = mod ( ti1 + 1, 3 ) + 1

        nti1 = element_node(ti1,e)
        nti2 = element_node(ti2,e)
        nti3 = element_node(ti3,e)

        qi = 0.5D+00 * ( &
            ( xy(1,nti3) - xy(1,nti2) ) * ( yq - xy(2,nti2) ) &
          - ( xy(2,nti3) - xy(2,nti2) ) * ( xq - xy(1,nti2) ) ) / area
        dqidx = - 0.5D+00 * ( xy(2,nti3) - xy(2,nti2) ) / area
        dqidy =   0.5D+00 * ( xy(1,nti3) - xy(1,nti2) ) / area

        rhs = 2.0D+00 * r8_pi * r8_pi * sin ( r8_pi * xq ) * sin ( r8_pi * yq )

        b(nti1) = b(nti1) + area * wq * rhs * qi
!
!  Consider each basis function in the element.
!
        do tj1 = 1, 3

          tj2 = mod ( tj1,     3 ) + 1
          tj3 = mod ( tj1 + 1, 3 ) + 1

          ntj1 = element_node(tj1,e)
          ntj2 = element_node(tj2,e)
          ntj3 = element_node(tj3,e)

          qj = 0.5D+00 * ( &
              ( xy(1,ntj3) - xy(1,ntj2) ) * ( yq - xy(2,ntj2) ) &
            - ( xy(2,ntj3) - xy(2,ntj2) ) * ( xq - xy(1,ntj2) ) ) / area
          dqjdx = - 0.5D+00 * ( xy(2,ntj3) - xy(2,ntj2) ) / area
          dqjdy =   0.5D+00 * ( xy(1,ntj3) - xy(1,ntj2) ) / area

          a(nti1,ntj1) = a(nti1,ntj1) &
            + area * wq * ( dqidx * dqjdx + dqidy * dqjdy )

        end do

      end do

    end do

  end do
!
!  BOUNDARY CONDITIONS
!
!  If the K-th variable is at a boundary node, replace the K-th finite
!  element equation by a boundary condition that sets the variable to U(K).
!
  k = 0

  do j = 1, ny

    do i = 1, nx

      k = k + 1

      if ( i == 1 .or. &
           i == nx .or. &
           j == 1 .or. &
           j == ny ) then

        call exact ( xy(1,k), xy(2,k), u, dudx, dudy )

        a(k,1:node_num) = 0.0D+00
        a(k,k)          = 1.0D+00
        b(k)            = u

      end if
    end do
  end do
!
!  SOLVE the linear system A * X = B.
!
!  The solution X is actually returned in the space occupied by B.
!
  call r8ge_fs ( node_num, a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_POISSON_RECTANGLE_LINEAR - Fatal error!'
    write ( *, '(a)' ) '  R8GE_FS returned an error condition.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The linear system was not solved, and the'
    write ( *, '(a)' ) '  algorithm cannot proceed.'
    stop
  end if
!
!  COMPARE computed and exact solutions.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     K     I     J          X           Y        U               ' // &
    'U               Error'
  write ( *, '(a)' ) '                                                ' // &
    ' exact           computed '
  write ( *, '(a)' ) ' '

  k = 0

  do j = 1, ny
    do i = 1, nx

      k = k + 1

      call exact ( xy(1,k), xy(2,k), u, dudx, dudy )
      write ( *, &
        '(2x,i4,2x,i4,2x,i4,2x,f10.2,2x,f10.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        k, i, j, xy(1,k), xy(2,k), u, b(k), abs ( u - b(k) )

    end do
    write ( *, '(a)' ) ' '
  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( element_node )
  deallocate ( xy )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_POISSON_RECTANGLE_LINEAR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine r8ge_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8GE_FS factors and solves a R8GE system.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_FS does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    R8GE_FS uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 1999
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
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the linear system.
!    On output, B is the solution of the linear system.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) temp

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

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      temp = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = temp

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol) = b(jcol) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine exact ( x, y, u, dudx, dudy )

!*****************************************************************************80
!
!! EXACT calculates the exact solution and its first derivatives.
!
!  Discussion:
!
!    The function specified here depends on the problem being
!    solved.  The user must be sure to change both EXACT and RHS
!    or the program will have inconsistent data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point
!    in the region, at which the exact solution is to be evaluated.
!
!    Output, real ( kind = 8 ) U, DUDX, DUDY, the value of
!    the exact solution U and its derivatives dUdX
!    and dUdY at the point (X,Y).
!
  implicit none

  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  u    =         sin ( r8_pi * x ) * sin ( r8_pi * y ) + x
  dudx = r8_pi * cos ( r8_pi * x ) * sin ( r8_pi * y ) + 1.0D+00
  dudy = r8_pi * sin ( r8_pi * x ) * cos ( r8_pi * y )

  return
end
