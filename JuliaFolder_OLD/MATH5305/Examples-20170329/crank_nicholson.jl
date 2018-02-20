"""
    U, x, t = crank_nicholson(L, T, κ, v, f, P, N)

Use the implicit Euler method to solve the initial-boundary
value problem for the heat equation

    uₜ - κ uₓₓ = f(x,t)    for 0 < t < T and 0 < x < L,

    u(0,t) = 0 = u(L,t)    for 0 < t < T,

    u(x,0) = v(x)          for 0 < x < L.

The finite difference grid has `P` subintervals in x, and `N`
subintervals in t.  Note:

    x[p+1] = xₚ,  t[n+1] = tₙ,  U[p+1,n+1] ≈ u(xₚ, tₙ)

for 0 ≤ p ≤ P and 0 ≤ n ≤ N.
"""
function crank_nicholson(L, T, κ, v, f, P, N)
    Δx = L / P
    Δt = T / N
    x = linspace(0, L, P+1)
    t = linspace(0, T, N+1)
    U = zeros(P+1, N+1)
    ρ = κ * Δt / Δx^2
    d = fill(ρ, P-1)
    du = fill(-ρ/2, P-2)
    IplusA = SymTridiagonal((1+ρ)*ones(d), du)
    IminusA = SymTridiagonal((1-ρ)*ones(d), -du)
    F = factorize(IplusA)
    rhs = Array(Float64, P-1)
    for p = 0:P
        U[p+1,1] = v(x[p+1])
    end
    for n = 1:N         # Time step
        rhs = IminusA*U[2:P,n]
        for p = 1:P-1           # Distance step
            rhs[p] += 0.5*(f(x[p+1],t[n+1])+f(x[p+1],t[n])) * Δt
        end
        U[2:P,n+1] = F \ rhs
    end
    return U, x, t
end
