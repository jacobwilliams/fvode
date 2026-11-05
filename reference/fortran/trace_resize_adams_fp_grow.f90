!*****************************************************************************************
!>
!  Step trace: `CVodeResizeHistory` after every step (grow by one element each step, CV_ADAMS,
!  fixed point). Fortran counterpart of `reference/trace_resize_adams_fp_grow.c`.

    program trace_resize_adams_fp_grow

    use fvode
    use resize_driver

    implicit none

    call run_resize(2, CV_ADAMS, .false., 'reference/out/fvode_resize_adams_fp_grow.csv')

    end program trace_resize_adams_fp_grow
!*****************************************************************************************
