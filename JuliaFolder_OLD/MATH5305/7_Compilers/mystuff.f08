module mystuff

    use iso_fortran_env
    implicit none

    integer, parameter :: DP = REAL64

contains

    function quadform(A, v) result(ans)
    real(DP), intent(in) :: A(:,:), v(:)
    real(DP) :: ans
    real(DP) :: inner_sum
    integer :: i, j, n

    n = size(A, 1)
    if ( size(A, 2) /= n ) then
        stop "Array must be square"
    end if
    if ( size(v) /= n ) then
        stop "Vector size does not match matrix"
    end if
    ans = 0.0
    do j = 1, n
        inner_sum = 0.0
        do i = 1, n
            inner_sum = inner_sum + v(i) * a(i,j)
        end do
        ans = ans + inner_sum * v(j)
    end do

end function quadform

end module mystuff
