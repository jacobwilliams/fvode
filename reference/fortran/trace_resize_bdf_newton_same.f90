!*****************************************************************************************
!>
!  Step trace: `CVodeResizeHistory` after every step (same size, CV_BDF,
!  Newton + band). Fortran counterpart of `reference/trace_resize_bdf_newton_same.c`.

    program trace_resize_bdf_newton_same

    use fvode
    use resize_driver

    implicit none

    call run_resize(1, CV_BDF, .true., 'reference/out/fvode_resize_bdf_newton_same.csv')

    end program trace_resize_bdf_newton_same
!*****************************************************************************************
