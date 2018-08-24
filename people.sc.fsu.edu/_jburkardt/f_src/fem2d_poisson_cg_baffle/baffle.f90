subroutine dirichlet_condition ( node_num, node_xy, node_bc )

!*****************************************************************************80
!
!! DIRICHLET_CONDITION sets the value of a Dirichlet boundary condition.
!
!  Discussion:
!
!    This routine specifies that the solution is zero on the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM),
!    the coordinates of the points.
!
!    Output, real ( kind = 8 ) NODE_BC(NODE_NUM), the value of the
!    Dirichlet boundary conditions at the points.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) node_bc(node_num)
  real ( kind = 8 ) node_xy(2,node_num)

  node_bc(1:node_num) = 0.0D+00

  return
end
subroutine h_coef ( node_num, node_xy, node_h )

!*****************************************************************************80
!
!! H_COEF evaluates the coefficient H(X,Y) of DEL U in the Poisson equation.
!
!  Discussion:
!
!    The equation is
!
!      - Del H(X,Y) DEL U + K(X,Y) * U = RHS(X,Y)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM),
!    the coordinates of the points.
!
!    Output, real ( kind = 8 ) NODE_H(NODE_NUM),
!    the value of the coefficient of DEL U.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) node_h(node_num)
  real ( kind = 8 ) node_xy(2,node_num)

  node_h(1:node_num) = 1.0D+00

  return
end
subroutine k_coef ( node_num, node_xy, node_k )

!*****************************************************************************80
!
!! K_COEF evaluates the coefficient K(X,Y) of U in the Poisson equation.
!
!  Discussion:
!
!    The equation is
!
!      - Del H(X,Y) DEL U + K(X,Y) * U = RHS(X,Y)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM),
!    the coordinates of the points.
!
!    Output, real ( kind = 8 ) NODE_K(NODE_NUM),
!    the value of the coefficient of U.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) node_k(node_num)
  real ( kind = 8 ) node_xy(2,node_num)

  node_k(1:node_num) = 0.0D+00

  return
end
subroutine rhs ( node_num, node_xy, node_rhs )

!*****************************************************************************80
!
!! RHS gives the right-hand side of the differential equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM),
!    the coordinates of the points.
!
!    Output, real ( kind = 8 ) NODE_RHS(NODE_NUM), the value of the right
!    hand side of the differential equation at (X,Y).
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) j
  real ( kind = 8 ) node_rhs(node_num)
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) z

  do j = 1, node_num
    z = 4.0D+00 - sqrt ( ( node_xy(1,j) - 2.0 )**2 & 
                       + ( node_xy(2,j) - 2.0 )**2 )
    node_rhs(j) = max ( z, 0.0D+00 )
  end do

  return
end
