using PyPlot

const κ = 0.2
const L = 1.0
const T = 0.8
const P = 30
const N = ceil(Integer, (κ*T/0.45)*(P/L)^2)

include("crank_nicholson.jl")

v(x) = (0 <= x <= L/2) ? 1.0 : 0.0
f(x,t) = 0


U, x, t = crank_nicholson(L, T, κ, v, f, P, N)

mesh(x, t, U', rstride=div(N,P))
xlabel(L"$x$", fontsize=14)
ylabel(L"$t$", fontsize=14)
ax = gca()
ax[:view_init](elev=30, azim=50)
show()
