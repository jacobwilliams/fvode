!*****************************************************************************************
!>
!  Step trace: van der Pol, Adams + damped Anderson acceleration (m = 3, up to 5 iterations, beta = 0.75).
!  Fortran counterpart of `reference/trace_vdp_adams_aa_damped.c`.

    program trace_vdp_adams_aa_damped
    use fvode_kinds, only: wp
    use vdp_adams_driver, only: run_vdp
    implicit none
    call run_vdp(3, 5, 0.75_wp, 'reference/out/fvode_vdp_adams_aa_damped.csv')
    end program trace_vdp_adams_aa_damped
!*****************************************************************************************
