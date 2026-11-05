!*****************************************************************************************
!>
!  FVODE: a Modern Fortran port of CVODE (SUNDIALS v7.9.0).
!
!  This is the umbrella module: `use fvode` gives access to the whole public API.

    module fvode

    use fvode_kinds
    use fvode_parameters

    implicit none

    public

    character(len=*),parameter :: fvode_version    = '0.0.1'  !! FVODE version
    character(len=*),parameter :: sundials_version = '7.9.0'  !! upstream SUNDIALS version ported

    end module fvode
!*****************************************************************************************
