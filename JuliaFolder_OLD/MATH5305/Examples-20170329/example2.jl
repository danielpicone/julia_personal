using PyPlot

const κ = 0.2
const L = 1.0
const T = 5.0
const P = 30
const N = ceil(Integer, (κ*T/0.45)*(P/L)^2)

include("ieuler.jl")

v(x) = 0.0
#f(x,t) = (1.0<=t<=3.0) ? (t-1.0)^2*(3.0-t)^2 : 0.0
f(x,t) = (1.0<=t<=3.0) ? (sin(t)^2-1.0) : 0.0


U, x, t = implicit_Euler(L, T, κ, v, f, P, N)

mesh(x, t, U', rstride=div(N,P))
xlabel(L"$x$", fontsize=14)
ylabel(L"$t$", fontsize=14)
ax = gca()
ax[:view_init](elev=30, azim=50)
show()
