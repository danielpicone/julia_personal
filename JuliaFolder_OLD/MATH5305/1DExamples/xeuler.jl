"""
    U, x, t = explicit_Euler(L, T, κ, v, f, P, N)

Solve the initial-boundary value problem for the heat equation 

    uₜ - κ uₓₓ = f(x,t)    for 0 < t < T and 0 < x < L,

    u(0,t) = 0 = u(L,t)    for 0 < t < T,

    u(x,0) = v(x)          for 0 < x < L.

The finite difference grid has `P` subintervals in x, and `N`
subintervals in t.  Note:

    x[p+1] = xₚ,  t[n+1] = tₙ,  U[p+1,n+1] ≈ u(xₚ, tₙ)

for 0 ≤ p ≤ P and 0 ≤ n ≤ N.
"""
function explicit_Euler(L, T, κ, v, f, P, N)
    Δx = L / P
    Δt = T / N
    x = linspace(0, L, P+1)
    t = linspace(0, T, N+1)
    U = zeros(P+1, N+1)
    ρ = κ * Δt / Δx^2
    if ρ > 1/2
        warn("Stability condition violated")
    end
    for p = 1:P-1
        U[p+1,1] = v(x[p+1])
    end
    for n = 0:N-1, p = 1:P-1
        U[p+1,n+2] = f(x[p+1],t[n+1])*Δt + ρ * U[p,n+1] +
                     (1-2ρ) * U[p+1,n+1] + ρ * U[p+2,n+1]
    end
    return U, x, t
end

