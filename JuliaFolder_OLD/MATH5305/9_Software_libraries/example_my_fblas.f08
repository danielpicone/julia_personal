program test_daxpy

    use my_fblas

    real(DP), dimension(4) :: x = [ 1.0, 0.0, -3.0, 2.0 ], &
                              y = [ -2.0, 1.0, 5.0, 0.0 ]
    real(DP) :: a = 1.5

    print *, "a = ", a
    call print_vector(x, "x = ", "(f8.4)")
    call print_vector(y, "y = ", "(f8.4)")
    call daxpy(size(x), a, x, 1, y, 1)
    call print_vector(y, "ax + y = ", "(f8.4)")

contains

    subroutine print_vector(x, label, fmt_str)
        real(DP),         intent(in) :: x(:)
        character(len=*), intent(in) :: label, fmt_str

        integer :: j
        
        write(*, "(a10)", advance='no') label
        do j = 1, size(x)
            write(*, fmt_str, advance='no') x(j)
        end do
        print *

    end subroutine print_vector

end program test_daxpy
