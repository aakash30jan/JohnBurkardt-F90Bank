program main

!*****************************************************************************80
!
!! MAIN is the main program for HB_READ_PRB.
!
!  Discussion:
!
!    HB_READ_PRB tests the HB_READ library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) filename

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'HB_READ_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HB_READ library.'
!
!  Quick read real single precision matrices.
!
  filename = 'g05_rua.txt'
  call test01 ( filename )
!
!  Quick read real double precision matrices.
!
  filename = 'g05_rua.txt'
  call test02 ( filename )
!
!  Quick read complex single precision matrices.
!
  filename = 'cg05_cua.txt'
  call test03 ( filename )
!
!  Quick read complex single precision matrices.
!
  filename = 'cg05_cua.txt'
  call test04 ( filename )
!
!  Full read real single precision matrices.
!
  filename = 'g05_rua.txt'
  call test05 ( filename )

  filename = 'g10_rua.txt'
  call test05 ( filename )

  filename = 'g20_rua.txt'
  call test05 ( filename )

  filename = 'cavity_small_rua.txt'
  call test05 ( filename )

  filename = 'cavity_big_rua.txt'
  call test05 ( filename )
!
!  Full read real double precision matrices.
!
  filename = 'g05_rua.txt'
  call test06 ( filename )
!
!  Full read complex single precision matrices.
!
  filename = 'cg20_cua.txt'
  call test07 ( filename )
!
!  Full read complex double precision matrices.
!
  filename = 'cg20_cua.txt'
  call test08 ( filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'HB_READ_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine test01 ( filename )

!*****************************************************************************80
!
!! TEST01 tests R4_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30
  integer ( kind = 4 ), parameter :: nnzero_max = 200

  integer ( kind = 4 ) colptr(n_max+1)
  character ( len = * ) filename
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero_max)
  real ( kind = 4 ) values(nnzero_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R4_HB_QUICK_READ reads a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using real single precision arithmetic.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  We assume NNZERO is no greater than ', nnzero_max
  write ( *, '(a,i6)' ) '  We assume N is no greater than ', n_max
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the information.
!
  call r4_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, colptr )
!
!  #3: Print (some of) the information.
!
  call r4_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
    colptr )
!
!  #4: Close the file.
!
  close ( unit = input )

  return
end
subroutine test02 ( filename )

!*****************************************************************************80
!
!! TEST02 tests R8_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30
  integer ( kind = 4 ), parameter :: nnzero_max = 200

  integer ( kind = 4 ) colptr(n_max+1)
  character ( len = * ) filename
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero_max)
  real ( kind = 4 ) values(nnzero_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8_HB_QUICK_READ reads a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using real double precision arithmetic.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  We assume NNZERO is no greater than ', nnzero_max
  write ( *, '(a,i6)' ) '  We assume N is no greater than ', n_max
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST02 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the information.
!
  call r4_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, colptr )
!
!  #3: Print (some of) the information.
!
  call r4_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
    colptr )
!
!  #4: Close the file.
!
  close ( unit = input )

  return
end
subroutine test03 ( filename )

!*****************************************************************************80
!
!! TEST03 tests C4_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30
  integer ( kind = 4 ), parameter :: nnzero_max = 200

  integer ( kind = 4 ) colptr(n_max+1)
  character ( len = * ) filename
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero_max)
  complex ( kind = 4 ) values(nnzero_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  C4_HB_QUICK_READ reads a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using complex single precision arithmetic.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  We assume NNZERO is no greater than ', nnzero_max
  write ( *, '(a,i6)' ) '  We assume N is no greater than ', n_max
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST03 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the information.
!
  call c4_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, colptr )
!
!  #3: Print (some of) the information.
!
  call c4_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
    colptr )
!
!  #4: Close the file.
!
  close ( unit = input )

  return
end
subroutine test04 ( filename )

!*****************************************************************************80
!
!! TEST04 tests C8_HB_QUICK_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30
  integer ( kind = 4 ), parameter :: nnzero_max = 200

  integer ( kind = 4 ) colptr(n_max+1)
  character ( len = * ) filename
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rowind(nnzero_max)
  complex ( kind = 8 ) values(nnzero_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  C8_HB_QUICK_READ reads a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using complex single precision arithmetic.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  We assume NNZERO is no greater than ', nnzero_max
  write ( *, '(a,i6)' ) '  We assume N is no greater than ', n_max
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST04 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the information.
!
  call c8_hb_quick_read ( input, nrow, ncol, nnzero, values, rowind, colptr )
!
!  #3: Print (some of) the information.
!
  call c8_hb_quick_print ( filename, nrow, ncol, nnzero, values, rowind, &
    colptr )
!
!  #4: Close the file.
!
  close ( unit = input )

  return
end
subroutine test05 ( filename )

!*****************************************************************************80
!
!! TEST05 tests R4_HB_HEADER_READ and R4_HB_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), allocatable :: colptr(:)
  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  real ( kind = 4 ), allocatable :: rhsval(:,:)
  integer ( kind = 4 ), allocatable :: rowind(:)
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 4 ), allocatable :: values(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  R4_HB_HEADER_READ and R4_HB_DATA_READ'
  write ( *, '(a)' ) '  read a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using real single precision arithmetic.'
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the header.
!
  call r4_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #3: Print the header.
!
  call r4_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #4: Allocate space.
!
  allocate ( colptr(1:ncol+1) )

  if ( 0 < nrhs ) then
    allocate ( rhsval(nrow,nrhs) )
  end if

  if ( 0 < nnzero ) then
    allocate ( rowind(1:nnzero) )
    allocate ( values(1:nnzero) )
  end if
!
!  #5: Read the structure and data.
!
  call r4_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
    rhsval )
!
!  #6: Print (some of) the structure.
!
  call r4_hb_structure_print ( ncol, nnzero, colptr, rowind )
!
!  #7: Print (some of) the data.
!
  call r4_hb_values_print ( ncol, colptr, nnzero, values )
!
!  #8: Print (some of) the right hand sides.
!
  if ( 0 < nrhs ) then
    call r4mat_print_some ( nrow, nrhs, rhsval, 1, 1, 10, 5, &
      '  10x5 portion of right hand sides:' )
  end if
!
!  #9: Close the file.
!
  close ( unit = input )

  return
end
subroutine test06 ( filename )

!*****************************************************************************80
!
!! TEST06 tests R8_HB_HEADER_READ and R8_HB_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), allocatable :: colptr(:)
  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  real ( kind = 8 ), allocatable :: rhsval(:,:)
  integer ( kind = 4 ), allocatable :: rowind(:)
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  real ( kind = 8 ), allocatable :: values(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  R8_HB_HEADER_READ and R8_HB_DATA_READ'
  write ( *, '(a)' ) '  read a sparse matrix from a Harwell-Boeing file'
  write ( *, '(a)' ) '  using real double precision arithmetic.'
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST06 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the header.
!
  call r8_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #3: Print the header.
!
  call r8_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #4: Allocate space.
!
  allocate ( colptr(1:ncol+1) )

  if ( 0 < nrhs ) then
    allocate ( rhsval(nrow,nrhs) )
  end if

  if ( 0 < nnzero ) then
    allocate ( rowind(1:nnzero) )
    allocate ( values(1:nnzero) )
  end if
!
!  #5: Read the structure and data.
!
  call r8_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
    rhsval )
!
!  #6: Print (some of) the structure.
!
  call r8_hb_structure_print ( ncol, nnzero, colptr, rowind )
!
!  #7: Print (some of) the data.
!
  call r8_hb_values_print ( ncol, colptr, nnzero, values )
!
!  #8: Print (some of) the right hand sides.
!
  if ( 0 < nrhs ) then
    call r8mat_print_some ( nrow, nrhs, rhsval, 1, 1, 10, 5, &
      '  10x5 portion of right hand sides:' )
  end if
!
!  #9: Close the file.
!
  close ( unit = input )

  return
end
subroutine test07 ( filename )

!*****************************************************************************80
!
!! TEST07 tests C4_HB_HEADER_READ and C4_HB_DATA_READ.
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
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), allocatable :: colptr(:)
  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  complex ( kind = 4 ), allocatable :: rhsval(:,:)
  integer ( kind = 4 ), allocatable :: rowind(:)
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 4 ), allocatable :: values(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  C4_HB_HEADER_READ and C4_HB_DATA_READ'
  write ( *, '(a)' ) '  read a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using complex single precision arithmetic.'
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST07 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the header.
!
  call c4_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #3: Print the header.
!
  call c4_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #4: Allocate space.
!
  allocate ( colptr(1:ncol+1) )

  if ( 0 < nrhs ) then
    allocate ( rhsval(nrow,nrhs) )
  end if

  if ( 0 < nnzero ) then
    allocate ( rowind(1:nnzero) )
    allocate ( values(1:nnzero) )
  end if
!
!  #5: Read the structure and data.
!
  call c4_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
    rhsval )
!
!  #6: Print (some of) the structure.
!
  call c4_hb_structure_print ( ncol, nnzero, colptr, rowind )
!
!  #7: Print (some of) the data.
!
  call c4_hb_values_print ( ncol, colptr, nnzero, values )
!
!  #8: Print (some of) the right hand sides.
!
  if ( 0 < nrhs ) then
    call c4mat_print_some ( nrow, nrhs, rhsval, 1, 1, 10, 5, &
      '  10x5 portion of right hand sides:' )
  end if
!
!  #9: Close the file.
!
  close ( unit = input )

  return
end
subroutine test08 ( filename )

!*****************************************************************************80
!
!! TEST08 tests C8_HB_HEADER_READ and C8_HB_DATA_READ.
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
!    Input, character ( len = * ) FILENAME, the name of the HB file to read.
!
  implicit none

  integer ( kind = 4 ), allocatable :: colptr(:)
  character ( len = * ) filename
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrhsix
  integer ( kind = 4 ) nrow
  character ( len = 16 ) ptrfmt
  integer ( kind = 4 ) rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  complex ( kind = 8 ), allocatable :: rhsval(:,:)
  integer ( kind = 4 ), allocatable :: rowind(:)
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
  complex ( kind = 8 ), allocatable :: values(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  C8_HB_HEADER_READ and C8_HB_DATA_READ'
  write ( *, '(a)' ) '  read a sparse matrix from a Harwell-Boeing file,'
  write ( *, '(a)' ) '  using complex double precision arithmetic.'
!
!  #1: Open the file.
!
  call get_unit ( input )
 
  open ( unit = input, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TEST08 - Warning!'
    write ( *, '(a)' ) '  Could not open file "' // trim ( filename ) // '".'
    return
  end if
!
!  #2: Read the header.
!
  call c8_hb_header_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #3: Print the header.
!
  call c8_hb_header_print ( filename, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt )
!
!  #4: Allocate space.
!
  allocate ( colptr(1:ncol+1) )

  if ( 0 < nrhs ) then
    allocate ( rhsval(nrow,nrhs) )
  end if

  if ( 0 < nnzero ) then
    allocate ( rowind(1:nnzero) )
    allocate ( values(1:nnzero) )
  end if
!
!  #5: Read the structure and data.
!
  call c8_hb_data_read ( input, nrow, ncol, nnzero, rhstyp, nrhs, nrhsix, &
    valcrd, rhscrd, ptrfmt, indfmt, valfmt, rhsfmt, colptr, rowind, values, &
    rhsval )
!
!  #6: Print (some of) the structure.
!
  call c8_hb_structure_print ( ncol, nnzero, colptr, rowind )
!
!  #7: Print (some of) the data.
!
  call c8_hb_values_print ( ncol, colptr, nnzero, values )
!
!  #8: Print (some of) the right hand sides.
!
  if ( 0 < nrhs ) then
    call c8mat_print_some ( nrow, nrhs, rhsval, 1, 1, 10, 5, &
      '  10x5 portion of right hand sides:' )
  end if
!
!  #9: Close the file.
!
  close ( unit = input )

  return
end
