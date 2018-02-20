module my_fblas

use iso_fortran_env
use iso_c_binding
implicit none

integer, parameter :: DP = REAL64

contains

subroutine daxpy(n, da, dx, incx, dy, incy)
    integer,  intent(in)    :: n, incx, incy
    real(DP), intent(in)    :: da, dx(*)
    real(DP), intent(inout) :: dy(*)

    integer :: ix, iy, i

    ix = 1
    iy = 1
    do i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx
        iy = iy + incy
    end do

end subroutine daxpy

subroutine impersonate_c_daxpy(n, da, dx, incx, dy, incy) bind(c)
    integer(C_INT), intent(in), value :: n, incx, incy
    real(C_DOUBLE), intent(in), value :: da
    real(C_DOUBLE), intent(in)        :: dx(*)
    real(C_DOUBLE), intent(inout)     :: dy(*)

    call daxpy(n, da, dx, incx, dy, incy)

end subroutine

end module my_fblas
