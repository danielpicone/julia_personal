program maxeig
    !
    ! Usage: ./maxeig N alpha
    !
    use iso_fortran_env
    implicit none

    integer, parameter :: DP = REAL64, MAXITS = 10000

    real(DP), allocatable :: A(:,:), v(:), w(:)
    real(DP) :: alpha, h, lambda, start, finish, resid
    integer :: length, stat, N, i, j, iter
    character(len=10) :: arg
 
    call get_command_argument(1, arg, length, stat)
    if ( length > len(arg) .or. stat /= 0 ) then
        print *, "Error reading command-line argument 1"
        stop 1
    else
        read(arg,*,iostat=stat) N
        if ( stat /= 0 ) then
            print *, "Failed to parse command-line argument 1"
            stop 2
        end if
    end if
    call get_command_argument(2, arg, length, stat)
    if ( length > len(arg) .or. stat /= 0 ) then
        print *, "Error reading command-line argument 2"
        stop 3
    else
        read(arg,*,iostat=stat) alpha
        if ( stat /= 0 ) then
            print *, "Failed to parse command-line argument 2"
            stop 4
        end if
    end if
    print "('Started ', i0, '-by-', i0, ' matrix with alpha = ', f4.2)", &
           N, N, alpha

    call cpu_time(start)
    h = 1.0_DP / N
    allocate(A(N,N), v(N), w(N))
    do j = 1, N
        do i = 1, N
            A(i,j) = ( abs(i-j) + 1 )**(alpha-1) * h**alpha / gamma(alpha) 
        end do
    end do

    call power_method()
    call matvec(w, A, v)
    w = w - lambda * v 
    resid = norm2(w)
    deallocate(A, v, w)
    call cpu_time(finish)

    print "('Finished ', i0, '-by-', i0, ' matrix with alpha = ', f4.2)", &
           N, N, alpha
    print "(t4, 'Max eigenvalue = ', f0.14)", lambda
    print "(t4, 'Used ', i0, ' iterations')", iter
    print "(t4, '|| Av - lambda v ||_2 = ', e8.2)", resid
    print "(t4, 'Elapsed time = ', f0.4, ' seconds')", finish-start

contains

    subroutine power_method()
        real(DP) :: prev_lambda, rtol

        rtol = 10 * epsilon(lambda)
        call random_number(v)                
        call matvec(w, A, v) ! w = Av
        prev_lambda = dot_product(w,v) / dot_product(v,v)
        v = w
        do iter = 1, MAXITS
            call matvec(w, A, v)
            lambda = dot_product(w,v) / dot_product(v,v)
            if ( abs(lambda-prev_lambda) < rtol * lambda ) exit
            v = w / norm2(w)
            prev_lambda = lambda
        end do

    end subroutine power_method

    subroutine matvec(y, M, x)
        ! y = M x
        real(DP), intent(out) :: y(:)
        real(DP), intent(in) :: M(:,:)
        real(DP), intent(in) :: x(:)

        integer :: nrows, ncols, i, j

        nrows = size(M, 1)
        ncols = size(M, 2)
        y = 0
        do j = 1, ncols
            do i = 1, nrows
                y(i) = y(i) + M(i,j)*x(j)
            end do
        end do

    end subroutine matvec

end program maxeig
