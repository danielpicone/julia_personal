# This is a script to test the MCIntegration function

include("areaFind.jl")
using AreaFind
using PyPlot

a = 0
b = π/2
N = 10000
#f(x) = 3 + x + x^2
f(x) = sin(x)

MC, randX, randY,x = MCIntegration(a,b,N,f)
figure(1)
plot(randX,randY,".",x,f.(x),"k",linewidth = 2)
