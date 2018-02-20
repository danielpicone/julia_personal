using PyPlot

const κ = 0.5
const L = 1.0
const T = 4.0
const P = 20
const N = ceil(Integer, (κ*T/0.45)*(P/L)^2)

include("xeuler.jl")

v(x) = x .* ( 1 - x )
f(x,t) = exp(-t)

U, x, t = explicit_Euler(L, T, κ, v, f, P, N)

mesh(x, t, U', rstride=div(N,P))
xlabel(L"$x$", fontsize=14)
ylabel(L"$t$", fontsize=14)
ax = gca()
ax[:view_init](elev=30, azim=50)
show()
