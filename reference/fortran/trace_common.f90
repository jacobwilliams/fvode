!*****************************************************************************************
!>
!  Step-trace output for the FVODE trace drivers: the same CSV columns as
!  `reference/trace_common.h`, so `reference/compare.py` can compare them.

    module trace_common

    use fvode_kinds, only: wp
    use fvode_core, only: cvode_solver

    implicit none

    private

    public :: trace_header, trace_row

    contains
!*****************************************************************************************

    subroutine trace_header(iunit, neq)
    !! write the CSV header line
    integer,intent(in) :: iunit  !! output unit
    integer,intent(in) :: neq  !! number of equations
    integer :: i  !! component counter
    write(iunit, '(a)', advance='no') 'nst,t,hu,qu,hcur,qcur,nfe,nni,nncf,netf,nsetups,njev,nfels,flag'
    do i = 1, neq
        write(iunit, '(a,i0)', advance='no') ',y', i
    end do
    write(iunit, '(a)') ''
    end subroutine trace_header

    subroutine trace_row(iunit, cv, flag, t, y)
    !! write the CSV row for the current step
    integer,intent(in) :: iunit  !! output unit
    class(cvode_solver),intent(in) :: cv  !! the CVODE solver
    integer,intent(in) :: flag  !! return flag of the last solve call
    real(wp),intent(in) :: t  !! time returned by the solver
    real(wp),dimension(:),intent(in) :: y  !! solution returned by the solver
    integer :: i  !! component counter
    write(iunit, '(i0,",",es25.17e3,",",es25.17e3,",",i0,",",es25.17e3,",",i0,8(",",i0))', &
          advance='no') &
        cv%get_num_steps(), t, cv%get_last_step(), cv%get_last_order(), &
        cv%get_current_step(), cv%get_current_order(), cv%get_num_rhs_evals(), &
        cv%get_num_nonlin_solv_iters(), cv%get_num_nonlin_solv_conv_fails(), &
        cv%get_num_err_test_fails(), cv%get_num_lin_solv_setups(), cv%get_num_jac_evals(), &
        cv%get_num_lin_rhs_evals(), flag
    do i = 1, size(y)
        write(iunit, '(",",es25.17e3)', advance='no') y(i)
    end do
    write(iunit, '(a)') ''
    end subroutine trace_row

    end module trace_common
!*****************************************************************************************
