program example2

    use mystuff, only: DP, quadform
    implicit none

    integer, parameter :: M = 1000
    real(DP) :: A(M, M), v(M), ans1, ans2

    print "(i0, '-by-', i0, ' quadratic form.')", M, M
    call random_number(A)
    call random_number(v)
    ans1 = dot_product(v, matmul(A,v))
    print "('Fortran intrinsics: ans = ', e20.14)", ans1
    ans2 = quadform(A, v)
    print "('quadform function:  ans = ', e20.14)", ans2

end program example2
