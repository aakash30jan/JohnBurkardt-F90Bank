      program main

!*****************************************************************************80
!
!! MAIN is the main program for TSG_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2014
!
!  Author:
!
!    Miro Stoyanov
!

!
!  Include the interace module.
!  At this point, I am not sure how to properly allow more than one instance.
!
  use tsgFGrid 

  integer i
  integer num_dimensions
  integer num_points
  integer num_outputs
  double precision, pointer :: points(:,:)
      
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSG_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TSG library.'
        
  call tsgFMakeLocalPolynomialGrid ( 2, 1, 4, 1, 1 )
      
  call tsgFGetPoints ( points )
      
  num_points = tsgFGetNumPoints ( )
  num_dimensions = tsgFGetNumDimensions ( )
  num_outputs = tsgFGetNumOutputs ( )
    
  write ( *, '(a,i6)' ) '  Number of Points ', num_points
  write ( *, '(a,i6)' ) '  Number of Dimensions ', num_dimensions
  write ( *, '(a,i6)' ) '  Number of Outputs ', num_outputs
      
  do i = 1, num_points
    write ( *, * ) i, points(i,1:num_dimensions)
  end do

! call tsgDeleteGrid()
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSG_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
