!*****************************************************************************************
!>
!  Step trace: van der Pol, Adams + plain fixed-point iteration.
!  Fortran counterpart of `reference/trace_vdp_adams_fp.c`.

    program trace_vdp_adams_fp
    use fvode_kinds, only: wp
    use vdp_adams_driver, only: run_vdp
    implicit none
    call run_vdp(0, 0, 1.0_wp, 'reference/out/fvode_vdp_adams_fp.csv')
    end program trace_vdp_adams_fp
!*****************************************************************************************
