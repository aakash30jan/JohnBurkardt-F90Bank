program main

!*****************************************************************************80
!
!! MAIN is the main program for IHS_DATASET.
!
!  Discussion:
!
!    IHS_DATASET generates an IHS dataset and writes it to a file.
!
!    This program is meant to be used interactively.  It's also
!    possible to prepare a simple input file beforehand and use it
!    in batch mode.
!
!    The program requests input values from the user:
!
!    * DIM_NUM, the spatial dimension,
!    * N, the number of points to generate,
!    * D, the duplication factor,
!    * SEED, a seed for the random number generator.
!
!    The program generates the data, writes it to the file
!
!      ihs_M_N.txt
!
!    where "M" and "N" are the numeric values specified by the user,
!    and then asks the user for more input.   To indicate that no further
!    computations are desired, it is enough to input a nonsensical
!    value, such as -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) dim_num
  character ( len = 255 ) :: file_out_name
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: i4_ihs
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r8_ihs

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IHS_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate an IHS dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is meant to be used interactively.'
  write ( *, '(a)' ) '  It is also possible to prepare a simple input'
  write ( *, '(a)' ) '  file beforehand and use it in batch mode.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * DIM_NUM, the spatial dimension,'
  write ( *, '(a)' ) '  * N, the number of points to generate,'
  write ( *, '(a)' ) '  * D, the duplication factor,'
  write ( *, '(a)' ) '  * SEED, a seed for the random number generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program generates the data, writes it to the file'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ihs_DIM_NUM_N.txt'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where "DIM_NUM" and "N" are the values specified'
  write ( *, '(a)' ) '  by the user, and then asks the user for more input.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To indicate that no further computations are'
  write ( *, '(a)' ) '  desired, it is enough to input a nonsensical value,'
  write ( *, '(a)' ) '  such as -1.'

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter DIM_NUM, the spatial dimension:'
    write ( *, '(a)' ) '  (2 is a small typical value).'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) dim_num

    if ( ios /= 0 ) then
      exit
    end if

    if ( dim_num <= 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of points to generate:'
    write ( *, '(a)' ) '  (10 is a small typical value).'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) n

    if ( ios /= 0 ) then
      exit
    end if

    if ( n <= 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter D, the duplication factor.'
    write ( *, '(a)' ) '  This must be at least 1, but not too large.'
    write ( *, '(a)' ) '  (5 is a typical value.)'
    write ( *, '(a)' ) '  (a negative or 0 value terminates execution.)'

    read ( *, *, iostat = ios ) d

    if ( ios /= 0 ) then
      exit
    end if

    if ( d <= 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SEED, the seed for the UNIFORM,'
    write ( *, '(a)' ) '  a portable random number generator:'
    write ( *, '(a)' ) '  (123456789 is a fun value.)'
    write ( *, '(a)' ) '  (0 means you want a seed to be chosen for you.)'
    write ( *, '(a)' ) '  (a negative value terminates execution.)'

    read ( *, *, iostat = ios ) seed

    if ( ios /= 0 ) then
      exit
    end if

    if ( seed < 0 ) then
      exit
    end if

    if ( seed == 0 ) then
      call get_seed ( seed )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The chosen SEED = ', seed
    end if
!
!  Compute integer IHS values.
!
    allocate ( i4_ihs(1:dim_num,1:n) )

    seed_init = seed

    call ihs ( dim_num, n, d, seed, i4_ihs )
!
!  Convert from integer in [1,n] to real in [0,1].
!
    allocate ( r8_ihs(1:dim_num,1:n) )

    r8_ihs(1:dim_num,1:n) = &
      real ( 2 * i4_ihs(1:dim_num,1:n) - 1, kind = 8 ) / &
      real ( 2 * n, kind = 8 )
!
!  Write real data to a file.
!
    write ( file_out_name, '(a,i2.2,a,i5.5,a,i2.2,a,i10.10,a)' ) &
      'ihs_', dim_num, '_', n, '_', d, '_', seed_init, '.txt'

    call r8mat_write ( file_out_name, dim_num, n, r8_ihs )

    deallocate ( i4_ihs )
    deallocate ( r8_ihs )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The IHS data was written to the file "' &
       // trim ( file_out_name ) // '".'

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IHS_DATASET'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
