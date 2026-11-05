!*****************************************************************************************
!>
!  Step trace: `CVodeResizeHistory` after every step (grow by one element each step, CV_BDF,
!  Newton + band). Fortran counterpart of `reference/trace_resize_bdf_newton_grow.c`.

    program trace_resize_bdf_newton_grow

    use fvode
    use resize_driver

    implicit none

    call run_resize(2, CV_BDF, .true., 'reference/out/fvode_resize_bdf_newton_grow.csv')

    end program trace_resize_bdf_newton_grow
!*****************************************************************************************
