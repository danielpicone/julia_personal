# This is the solution to Quesiton 5
# We will compute the Maximum error and the rate of convergence by using a grid
# which gets finer

include("AllenCahn.jl")
using AllenCahn
using PyPlot
using OffsetArrays

const Ω = (-0.5, 1.5)
const ϵ = 0.015
const c = 3 / (sqrt(2)*ϵ)
const T = 1/c

function u(x, t)
    c = sqrt(2)*ϵ
    y = ( x - 3t/c ) / (2c)
    return ( 1 - tanh(y) ) / 2
end

numiter = 5
maxerr = zeros(numiter)

@printf("   P       N     Max error     rate      seconds\n\n")


for j=1:numiter
  P = 2^(j+5)
  N = 8P
  gr = grid1D(Ω,T,P,N)
  x, t = gr.x, gr.t
  u0 = u.(x,0.0)
  start = time()
  U = semi_implicit(ϵ,gr,u0)
  stop = time()
  err_ = maxabs(U[1:P,0:N] - u.(x[1:P], gr.t[0:N]'),1)
  maxerr[j] = maximum(err_)
  seconds = stop - start
  if j==1
    @printf("%4d    %4d     %2.3e               %1.3f\n",P,N,maxerr[j],seconds)
  else
    rn = log2(maxerr[j-1]/maxerr[j])
    @printf("%4d    %4d     %2.3e     %1.3f     %1.3f\n",P,N,maxerr[j],rn,seconds)
  end

end
