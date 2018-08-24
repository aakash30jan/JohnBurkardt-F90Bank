program main

!*****************************************************************************80
!
!! NSASM_TEST_SMALL shows how to call NSASM with the "small" data.
!
!  Discussion:
!
!    The flow region is the unit square [0,1]x[0,1] on which
!    we impose a 5 x 5 grid of nodes indexed as follows:
!
!    1.0
!     ^       7  19   8  23   9
!     |      21  20  22  24  25
!     Y       4  10   5  15   6
!     |      12  11  13  16  17
!     |       1  14   2   18  3
!     |
!     0-------------- X -->  1.0
!
!    The node coordinates are stored in the file 'small_nodes.txt'.
!    The pressure nodes, numbers 1 through 9, form a 3x3 subgrid.
!
!    The nodes are used to form 8 triangular elements:  
!
!     7--19---8--23---9
!     |     / |     / |
!    21  20  22  24  25
!     | /     | /     |
!     4--10---5--15---6
!     |     / |     / |
!    12  11  13  16  17
!     | /     | /     |
!     1--14---2--18---3
!
!    The 8 elements are quadratic six node triangles, with a peculiar
!    ordering for the nodes, and one which is not consistently
!    counterclockwise.
!
!      1: 1, 4, 5, 10, 11, 12.
!      2: 1, 2, 5, 13, 11, 14
!      3: 2, 5, 6, 15, 16, 13
!      4: 2, 3, 6, 17, 16, 18
!      5: 4, 7, 8, 19, 20, 21
!      6: 4, 5, 8, 22, 20, 10
!      7: 5, 8, 9, 23, 24, 22
!      8: 5, 6, 9, 25, 24, 15
!
!    The node indices forming each element are stored in 'small_elements.txt'.
!
!    All boundary nodes are constrained to have zero horizontal velocity.
!
!    Boundary nodes 3, 17, 6, 25 and 9 are constrained to have a vertical
!    velocity of 1; the other boundary nodes are constrained to have a
!    vertical velocity of 0.
!
!    Node 1 is constrained to have a pressure of 1.
!
!    The node index, the variable being set (0=horizontal velocity, 1=vertical,
!    2=pressure), and the value of the variable being constrained are stored 
!    as triples in 'small_constraints.txt'.
!
!    The number of degrees of freedom "NDOF", is:
!      2 * NP (horizontal and vertical velocity at every node)
!      + NP0 (pressure at every interior node)
!      + NE (constraints)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson,
!    Implementation of Finite Element-Based Navier-Stokes Solver,
!    April 2002
!
!  Local Parameters:
!
!    Local, string E_FILE, the name of the constraint file, which contains 
!    3 rows and NE columns, defining constraints on the data, including 
!    Dirichlet boundary values in particular.
!    Item #1 is a node, #2 is a variable index (0 = horizontal velocity,
!    1 = vertical velocity, 2 = pressure) and #3 is an associated value.
!
!    Local, sparse real K(*), the stiffness matrix, in compressed column format.
!
!    Local, real L(2*NP+NP0+NE), the residual vector.
!
!    Local, integer NE, the number of constraints.
!
!    Local, integer NP, the number of nodes.
!
!    Local, integer NP0, the number of pressure nodes (vertices of triangles.)
!
!    Local, integer NT, the number of triangles.
!
!    Local, real NU, the viscosity.
!
!    Local, string P_FILE, the name of the file that contains 2 rows and NP 
!    columns of (X,Y) node coordinates.
!
!    Local, string T_FILE, the name of the element file, which contains 6 rows 
!    and NT columns, with each column containing the (1-based) indices of nodes 
!    forming the triangles, in a particular order.
!
  implicit none

  character * ( 255 ) e_file
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: k_col(:)
  integer ( kind = 4 ) k_nz
  integer ( kind = 4 ), allocatable :: k_row(:)
  real ( kind = 8 ), allocatable :: k_val(:)
  real ( kind = 8 ), allocatable :: l(:)
  real ( kind = 8 ) l_norm
  integer ( kind = 4 ) ndof
  integer ( kind = 4 ) np0
  real ( kind = 8 ) nu
  character * ( 255 ) p_file
  real ( kind = 8 ) r8vec_amax
  character * ( 255 ) t_file
  real ( kind = 8 ), allocatable :: u(:)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NSASM_SMALL_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Call NSASM_INTERFACE with the "small" test data'
!
!  Define the input to NSASM_INTERFACE:
!
!    P_FILE contains 2 rows and NP columns of (X,Y) node coordinates.
!    T_FILE contains 6 rows and NT columns of node indices forming triangles. 
!    E_FILE contains 3 rows and NE columns, defining constraints.
!    NP0 is the number of pressure nodes.
!    NU is the fluid viscosity.
!
  p_file = 'small_nodes.txt'
  t_file = 'small_elements.txt'
  e_file = 'small_constraints.txt'
  np0 = 9
  nu = 100.0D+00
!
!  Call NSASM_INTERFACE to compute the system matrix and right hand side.
!
  call nsasm_interface ( p_file, t_file, e_file, np0, nu, ndof, k_nz, ...
    k_row, k_col, k_val, l, u )
!
!  Print some of L.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of degrees of freedom, = ', ndof

  l_norm = maxval ( abs ( l(1:ndof) ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L-Infinity norm of L = ', l_norm

  call r8vec_print_part ( ndof, l, 15, '  L vector (partial):' )
!
!  Get the number of nonzero entries in K.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Matrix nonzeros K_NZ = ', k_nz
!
!  Print (some of) the sparse triplet data.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Matrix sparse triplet representation:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       ROW       COL     VAL'
  write ( *, '(a)' ) ''

  do i = 1, k_nz

    if ( i < 10 .or. k_nz - 10 < i ) then
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) k_row(i), k_col(i), k_val(i)
    end if

    if ( i == 10 ) then
      write ( *, '(a)' ) '..(skipping some entries)...'
    end if

  end do
!
!  Print (some of) the matrix.
!
  call r8sp_print_some ( ndof, ndof, k_nz, k_col, k_row, k_val, 1, 1, 10, &
    10, '  Initial part of K as a matrix' )
!
!  Print some of U.
!
  call r8vec_print_part ( ndof, u, 15, '  U vector (partial):' )
!
!  Free memory.
!
  deallocate ( k_col )
  deallocate ( k_row )
  deallocate ( k_val )
  deallocate ( l )
  deallocate ( u )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NSASM_TEST_SMALL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end

