module Elliptic

  use iso_c_binding, only: C_DOUBLE, C_INT
  !$ use omp_lib
  implicit none

  integer, parameter :: NPROCS = 2

contains

subroutine ser_Richardson(U, r, nu, omega, Dx, Dy, maxP, maxQ, b) &
    bind(c, name="ser_Richardson")
    real(C_DOUBLE), intent(inout)     :: U(0:maxP,0:maxQ), r(maxP-1,maxQ-1)
    real(C_DOUBLE), intent(in), value :: omega, Dx, Dy
    integer(C_INT), intent(in), value :: nu, maxP, maxQ
    real(C_DOUBLE), intent(in)        :: b(maxP-1,maxQ-1)

    real(C_DOUBLE) :: rdx2, rdy2, resid
    integer :: k, p, q

    rdx2 = 1 / Dx**2
    rdy2 = 1 / Dy**2
    
    do k = 1, nu
        do q = 1, maxQ-1
            do p = 1, maxP-1
                r(p,q) = b(p,q) - (                                   &
                          rdx2 * (-U(p+1,q) + 2*U(p,q) - U(p-1,q) )   &
                        + rdy2 * (-U(p,q+1) + 2*U(p,q) - U(p,q-1) ) ) 
                U(p,q) = U(p,q) + omega * r(p,q)
            end do
        end do
    end do

end subroutine ser_Richardson

subroutine par_Richardson(U, r, nu, omega, Dx, Dy, maxP, maxQ, b) &
    bind(c, name="par_Richardson")
    real(C_DOUBLE), intent(inout)     :: U(0:maxP,0:maxQ), r(maxP-1,maxQ-1)
    real(C_DOUBLE), intent(in), value :: omega, Dx, Dy
    integer(C_INT), intent(in), value :: nu, maxP, maxQ
    real(C_DOUBLE), intent(in)        :: b(maxP-1,maxQ-1)

    real(C_DOUBLE) :: rdx2, rdy2, resid
    integer :: k, p, q

    rdx2 = 1 / Dx**2
    rdy2 = 1 / Dy**2
    
    !$omp parallel num_threads(NPROCS) default(none) private(p, q) &
    !$omp shared(U, r, b, nu, maxP, maxQ, omega, rdx2, rdy2)
    do k = 1, nu
        !$omp do
        do q = 1, maxQ-1
            do p = 1, maxP-1
                r(p,q) = b(p,q) - (                                   &
                          rdx2 * (-U(p+1,q) + 2*U(p,q) - U(p-1,q) )   &
                        + rdy2 * (-U(p,q+1) + 2*U(p,q) - U(p,q-1) ) ) 
                U(p,q) = U(p,q) + omega * r(p,q)
            end do
        end do
        !$omp end do
    end do
    !$omp end parallel

end subroutine par_Richardson

end module elliptic
