!*****************************************************************************************
!>
!  Numeric kind parameters for FVODE.
!
!  These replace the SUNDIALS `sunrealtype`, `sunindextype` and `long int`
!  typedefs. The real kind is selected at compile time, mirroring the SUNDIALS
!  `SUNDIALS_SINGLE_PRECISION` / `SUNDIALS_DOUBLE_PRECISION` /
!  `SUNDIALS_EXTENDED_PRECISION` options:
!
!  * default        : `real64`
!  * `-DREAL32`     : `real32`
!  * `-DREAL128`    : `real128`
!
!  Similarly `-DINT32` selects a 32-bit `sunindextype` (default is 64-bit).

    module fvode_kinds

    use iso_fortran_env, only: real32, real64, real128, int32, int64

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! working real kind (`sunrealtype`)
#elif REAL128
    integer,parameter,public :: wp = real128  !! working real kind (`sunrealtype`)
#else
    integer,parameter,public :: wp = real64   !! working real kind (`sunrealtype`)
#endif

#ifdef INT32
    integer,parameter,public :: ip = int32    !! index kind (`sunindextype`)
#else
    integer,parameter,public :: ip = int64    !! index kind (`sunindextype`)
#endif

    integer,parameter,public :: li = int64    !! kind for counters (C `long int`)

    end module fvode_kinds
!*****************************************************************************************
