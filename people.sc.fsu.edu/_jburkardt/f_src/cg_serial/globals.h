      include 'npbparams.h'

!---------------------------------------------------------------------
!  Note: please observe that in the routine conj_grad three
!  implementations of the sparse matrix-vector multiply have
!  been supplied.  The default matrix-vector multiply is not
!  loop unrolled.  The alternate implementations are unrolled
!  to a depth of 2 and unrolled to a depth of 8.  Please
!  experiment with these to find the fastest for your particular
!  architecture.  If reporting timing results, any of these three may
!  be used without penalty.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!  Class specific parameters:
!  It appears here for reference only.
!  These are their values, however, this info is imported in the npbparams.h
!  include file, which is written by the sys/setparams.c program.
!---------------------------------------------------------------------

!----------
!  Class S:
!----------
!!       parameter( na=1400, &
!!                  nonzer=7, &
!!                  shift=10., &
!!                  niter=15, &
!!                  rcond=1.0d-1 )
!----------
!  Class W:
!----------
!!       parameter( na=7000, &
!!                  nonzer=8, &
!!                  shift=12., &
!!                  niter=15,&
!!                  rcond=1.0d-1 )
!----------
!  Class A:
!----------
!!       parameter( na=14000,&
!!                  nonzer=11,&
!!                  shift=20.,&
!!                  niter=15,&
!!                  rcond=1.0d-1 )
!----------
!  Class B:
!----------
!!       parameter( na=75000, &
!!                  nonzer=13, &
!!                  shift=60., &
!!                  niter=75,&
!!                  rcond=1.0d-1 )
!----------
!  Class C:
!----------
!!       parameter( na=150000,&
!!                  nonzer=15, &
!!                  shift=110., &
!!                  niter=75,&
!!                  rcond=1.0d-1 )
!----------
!  Class D:
!----------
!!       parameter( na=1500000, &
!!                  nonzer=21, &
!!                  shift=500., &
!!                  niter=100,&
!!                  rcond=1.0d-1 )
!----------
!  Class E:
!----------
!!       parameter( na=9000000,&
!!                  nonzer=26, &
!!                  shift=1500.,&
!!                  niter=100,&
!!                  rcond=1.0d-1 )


      integer    nz, naz
      parameter( nz = na*(nonzer+1)*(nonzer+1) )
      parameter( naz = na*(nonzer+1) )


      common / partit_size  / 	 naa, nzz, &
                              	 firstrow, &
                                 lastrow, &
                                 firstcol, &
                                 lastcol
      integer                 	 naa, nzz, &
                              	 firstrow, &
                                 lastrow, &
                                 firstcol, &
                                 lastcol

      common /urando/         	 amult, tran
      double precision           amult, tran

      external         timer_read
      double precision timer_read

      integer T_init, T_bench, T_conj_grad, T_last
      parameter (T_init=1, T_bench=2, T_conj_grad=3, T_last=3)
      logical timeron
      common /timers/ timeron
