program example1

    use iso_fortran_env
    implicit none

    integer, parameter :: DP = REAL64

    integer, parameter :: M = 1000
    real(DP) :: A(M, M), v(M), ans1, ans2

    print "(i0, '-by-', i0, ' quadratic form.')", M, M
    call random_number(A)
    call random_number(v)
    ans1 = dot_product(v, matmul(A,v))
    print "('Fortran intrinsics: ans = ', e20.14)", ans1
    ans2 = quadform(A, v)
    print "('quadform function:  ans = ', e20.14)", ans2

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

end program example1
