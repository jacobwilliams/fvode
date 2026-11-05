!*****************************************************************************************
!>
!  Step trace: `CVodeResizeHistory` after every step (same size, CV_ADAMS,
!  fixed point). Fortran counterpart of `reference/trace_resize_adams_fp_same.c`.

    program trace_resize_adams_fp_same

    use fvode
    use resize_driver

    implicit none

    call run_resize(1, CV_ADAMS, .false., 'reference/out/fvode_resize_adams_fp_same.csv')

    end program trace_resize_adams_fp_same
!*****************************************************************************************
