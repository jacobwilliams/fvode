!*****************************************************************************************
!>
!  Test driver: runs all FVODE unit test suites (test-drive).

    program test_main

    use iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use test_parameters, only: collect_parameters
    use test_nvector, only: collect_nvector
    use test_core, only: collect_core
    use test_linsol, only: collect_linsol
    use test_krylov, only: collect_krylov

    implicit none

    integer :: stat, is  !! number of failed tests, and suite counter
    type(testsuite_type),allocatable :: testsuites(:)  !! all test suites
    character(len=*),parameter :: fmt = '("#", *(1x, a))'  !! format for the suite name

    stat = 0

    testsuites = [ new_testsuite("parameters", collect_parameters), &
                   new_testsuite("nvector",    collect_nvector),   &
                   new_testsuite("core",       collect_core),      &
                   new_testsuite("linsol",     collect_linsol),    &
                   new_testsuite("krylov",     collect_krylov)     ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

    end program test_main
!*****************************************************************************************
