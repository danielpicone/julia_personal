include("ieuler.jl")

const κ = 0.25
const L = 1.0
const T = 5.0

λ(j) = κ*(j*π/L)^2
ϕ(j,x) = sin(j*π*x/L)
v(x) = c1 * ϕ(1,x) + c3 * ϕ(3,x)
const c1 = 2.0
const c3 = -1.0
u(x,t) = c1 * exp(-λ(1)*t) * ϕ(1,x) + c3 * exp(-λ(3)*t) * ϕ(3,x)
f(x,t) = 0.0

function maxerror(U, x, t)
    err = 0.0
    N = length(t) - 1
    P = length(x) - 1
    for n = 2:N, p=1:P-1
        next = U[p+1,n+1] - u(x[p+1],t[n+1])
        err = max(err, abs(next))
    end
    return err
end

@printf("Errors for the implicit Euler method.\n\n")

nrows = 5
P = 1000
N = 40
err = zeros(nrows)
@printf("Mesh refinement in time.\n")
@printf("\n%4s  %8s  %6s  %8s  %6s\n\n", 
        "", "N", "P", "error", "rate")
for row = 1:nrows
    N *= 2
    U, x, t = implicit_Euler(L, T, κ, v, f, P, N)
    err[row] = maxerror(U, x, t)
    if row == 1
        @printf("%4d  %8d  %6d  %8.3e\n", 
                row, N, P, err[row])
    else
        rate = log2(err[row-1]/err[row])
        @printf("%4d  %8d  %6d  %8.3e  %6.3f\n", 
                row, N, P, err[row], rate)
    end
end

nrows = 5
N = 500_000
P = 5
err = zeros(nrows)
@printf("\n\nMesh refinement in space.\n")
@printf("\n%4s  %8s  %6s  %8s  %6s\n\n", 
        "", "N", "P", "error", "rate")
for row = 1:nrows
    P *= 2
    U, x, t = implicit_Euler(L, T, κ, v, f, P, N)
    err[row] = maxerror(U, x, t)
    if row == 1
        @printf("%4d  %8d  %6d  %8.3e\n", 
                row, N, P, err[row])
    else
        rate = log2(err[row-1]/err[row])
        @printf("%4d  %8d  %6d  %8.3e  %6.3f\n", 
                row, N, P, err[row], rate)
    end
end
