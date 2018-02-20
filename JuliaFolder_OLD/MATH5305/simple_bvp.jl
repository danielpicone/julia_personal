using PyPlot

"""
U, x = bvp(L, a, b, c, f, γ0, γL, P)
Uses central difference approximations on a grid with P subintervals
to solve the 2-point BVP

-a(x)u'' + b(x)u' + c(x)u = f(x)   for 0 ≤ x ≤ L,

u(0) = γ0,   u(L) = γL.
"""
function bvp(L, a, b, c, f, γ0, γL, P)
    Δx = L / P
    x = linspace(0, L, P+1)
    adx2 = Float64[ a(x[p+1]) / Δx^2 for p = 1:P-1 ]
    A = Tridiagonal(-adx2[2:P-1], 2adx2[1:P-1], -adx2[1:P-2])
    bdx = Float64[ b(x[p+1]) / (2Δx) for p = 1:P-1 ]
    B = Tridiagonal(-bdx[2:P-1], zeros(P-1), bdx[1:P-2])
    C = Diagonal(Float64[ c(x[p+1]) for p = 1:P-1 ])
    rhs = Float64[ f(x[p+1]) for p = 1:P-1 ]
    rhs[1]   += ( adx2[1]   + bdx[1]   ) * γ0
    rhs[P-1] += ( adx2[P-1] - bdx[P-1] ) * γL
    U = zeros(P+1)
    U[1]   = γ0
    U[2:P] = ( A + B + C ) \ rhs
    U[P+1] = γL
    return U, x
end

a(x) = exp(x)
b(x) = exp(-x)
c(x) = 1 + x^2
u(x) = 1 + sin(x)
f(x) =  a(x) * sin(x) + b(x) * cos(x) + c(x) * (1+sin(x))
γ0 = u(0.0)
L = π/2
γL = u(L)

P = 40

U, x = bvp(L, a, b, c, f, γ0, γL, P)

figure(1)
#plot(x,U)
plot(x, u(x), x, U, "o")
xlabel(L"$x$", fontsize=14)
grid(true)
legend(("exact soln", "finite difference approx");
       loc="lower right")
