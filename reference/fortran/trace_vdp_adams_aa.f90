!*****************************************************************************************
!>
!  Step trace: van der Pol, Adams + Anderson acceleration (m = 2, up to 5 iterations).
!  Fortran counterpart of `reference/trace_vdp_adams_aa.c`.

    program trace_vdp_adams_aa
    use fvode_kinds, only: wp
    use vdp_adams_driver, only: run_vdp
    implicit none
    call run_vdp(2, 5, 1.0_wp, 'reference/out/fvode_vdp_adams_aa.csv')
    end program trace_vdp_adams_aa
!*****************************************************************************************
