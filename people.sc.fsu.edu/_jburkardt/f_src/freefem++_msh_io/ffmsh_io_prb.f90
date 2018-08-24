program main

!*****************************************************************************80
!
!! MAIN is the main program for FFMSH_IO_PRB.
!
!  Discussion:
!
!    FFMSH_IO_PRB tests the FFMSH_IO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FFMSH_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FFMSH_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FFMSH_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 gets the example data and prints it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: e_l(:)
  integer ( kind = 4 ) e_num
  integer ( kind = 4 ), allocatable :: e_v(:,:)
  integer ( kind = 4 ), allocatable :: t_l(:)
  integer ( kind = 4 ) t_num
  integer ( kind = 4 ), allocatable :: t_v(:,:)
  integer ( kind = 4 ), allocatable :: v_l(:)
  integer ( kind = 4 ) v_num
  real ( kind = 8 ), allocatable :: v_xy(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Get the example 2D data and print it.'
!
!  Get the sizes.
!
  call ffmsh_2d_size_example ( v_num, e_num, t_num )
!
!  Print the sizes.
!
  call ffmsh_2d_size_print ( '  Example Sizes:', v_num, e_num, t_num ) 
!
!  Allocate memory.
!
  allocate ( v_xy(2,v_num) )
  allocate ( v_l(v_num) )
  allocate ( e_v(2,e_num) )
  allocate ( e_l(e_num) )
  allocate ( t_v(3,t_num) )
  allocate ( t_l(t_num) )
!
!  Get the data.
!
  call ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, &
    t_v, t_l )
!
!  Print the data.
!
  call ffmsh_2d_data_print ( '  Example data:', v_num, e_num, t_num, v_xy, &
    v_l, e_v, e_l, t_v, t_l )
!
!  Free memory.
!
  deallocate ( e_l )
  deallocate ( e_v )
  deallocate ( t_l )
  deallocate ( t_v )
  deallocate ( v_l )
  deallocate ( v_xy )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 writes the example data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: e_l(:)
  integer ( kind = 4 ) e_num
  integer ( kind = 4 ), allocatable :: e_v(:,:)
  character ( len = 255 ) filename
  integer ( kind = 4 ), allocatable :: t_l(:)
  integer ( kind = 4 ) t_num
  integer ( kind = 4 ), allocatable :: t_v(:,:)
  integer ( kind = 4 ), allocatable :: v_l(:)
  integer ( kind = 4 ) v_num
  real ( kind = 8 ), allocatable :: v_xy(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Write the example 2D data to a file.'
!
!  Get example sizes.
!
  call ffmsh_2d_size_example ( v_num, e_num, t_num )
!
!  Allocate memory.
!
  allocate ( v_xy(2,v_num) )
  allocate ( v_l(v_num) )
  allocate ( e_v(2,e_num) )
  allocate ( e_l(e_num) )
  allocate ( t_v(3,t_num) )
  allocate ( t_l(t_num) )
!
!  Get example data.
!
  call ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, &
    t_v, t_l )
!
!  Write the data to a file.
!
  filename = 'output.msh'

  call ffmsh_2d_write ( filename, v_num, e_num, t_num, v_xy, &
    v_l, e_v, e_l, t_v, t_l )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The data was written to "' // trim ( filename ) // '"'
!
!  Free memory.
!
  deallocate ( e_l )
  deallocate ( e_v )
  deallocate ( t_l )
  deallocate ( t_v )
  deallocate ( v_l )
  deallocate ( v_xy )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 gets the example data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: e_l(:)
  integer ( kind = 4 ) e_num
  integer ( kind = 4 ), allocatable :: e_v(:,:)
  character ( len = 255 ) ffmsh_filename
  integer ( kind = 4 ), allocatable :: t_l(:)
  integer ( kind = 4 ) t_num
  integer ( kind = 4 ), allocatable :: t_v(:,:)
  integer ( kind = 4 ), allocatable :: v_l(:)
  integer ( kind = 4 ) v_num
  real ( kind = 8 ), allocatable :: v_xy(:,:)

  ffmsh_filename = 'input.msh'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Read the example 2D data from a file.'
!
!  Read sizes.
!
  call ffmsh_2d_size_read ( ffmsh_filename, v_num, e_num, t_num )
!
!  Print sizes.
!
  call ffmsh_2d_size_print ( ffmsh_filename, v_num, e_num, t_num ) 
!
!  Allocate memory.
!
  allocate ( v_xy(2,v_num) )
  allocate ( v_l(v_num) )
  allocate ( e_v(2,e_num) )
  allocate ( e_l(e_num) )
  allocate ( t_v(3,t_num) )
  allocate ( t_l(t_num) )
!
!  Read data.
!
  call ffmsh_2d_data_read ( ffmsh_filename, v_num, e_num, t_num, v_xy, &
    v_l, e_v, e_l, t_v, t_l )
!
!  Print data.
!
  call ffmsh_2d_data_print ( ffmsh_filename, v_num, e_num, t_num, v_xy, &
    v_l, e_v, e_l, t_v, t_l )
!
!  Free memory.
!
  deallocate ( e_l )
  deallocate ( e_v )
  deallocate ( t_l )
  deallocate ( t_v )
  deallocate ( v_l )
  deallocate ( v_xy )

  return
end
