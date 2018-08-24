program main

!*****************************************************************************80
!
!! R83V_PRB tests the R83V library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the R83V library.'

  call i4_log_10_test ( )

  call r8_uniform_01_test ( )

  call r83_indicator_test ( )
  call r83_print_test ( )
  call r83_print_some_test ( )

  call r83v_cg_test ( )
  call r83v_copy_test ( )
  call r83v_cr_fa_test ( )
  call r83v_cr_sl_test ( )
  call r83v_cr_sls_test ( )
  call r83v_dif2_test ( )
  call r83v_fs_test ( )
  call r83v_gs_sl_test ( )
  call r83v_indicator_test ( )
  call r83v_jac_sl_test ( )
  call r83v_mtv_test ( )
  call r83v_mv_test ( )
  call r83v_print_test ( )
  call r83v_print_some_test ( )
  call r83v_random_test ( )
  call r83v_res_test ( )
  call r83v_to_r8ge_test ( )
  call r83v_to_r8vec_test ( )
  call r83v_transpose_test ( )
  call r83v_zeros_test ( )

  call r8ge_indicator_test ( )
  call r8ge_mtv_test ( )
  call r8ge_mv_test ( )
  call r8ge_print_test ( )
  call r8ge_print_some_test ( )
  call r8ge_to_r83v_test ( )

  call r8vec_indicator1_test ( )
  call r8vec_norm_test ( )
  call r8vec_norm_affine_test ( )
  call r8vec_print_test ( )
  call r8vec_to_r83v_test ( )
  call r8vec_uniform_01_test ( )

  call r8vec2_print_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end

