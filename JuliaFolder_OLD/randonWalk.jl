# This script is to graph a discrete rand walk

using PyPlot
N = 100
t = linspace(0,1,N)
S0 = 4.0
r = 0.001
u = 1.5
K = 2
b = 1

function symRandWalk(N::Int64,S0,u,K::Int64,b::Int64,r)
  M = zeros(N)
  maxVal = -Inf
  for i=2:N
    rand1 = rand()
    if rand1 < 0.5
      M[i] = M[i-1] - 1
    elseif rand1 >= 0.5
      M[i] = M[i-1] + 1
    end
    maxVal = max(maxVal,M[i])
  end


  if maxVal >= b
    ST = max(S0*u^(M[N]) - S0*u^K,0)/(1+r)^N
  else
    ST = 0
  end
  #println(maxVal)
  return M, maxVal, ST
end

function randWalk(N::Int64,S0,u,K::Int64,b::Int64,r)
  p = (1+r-u^-1)/(u-u^-1)
  #println(p)
  M = zeros(N)
  maxVal = -Inf
  for i=2:N
    rand1 = rand()
    if rand1 < 1-p
      M[i] = M[i-1] - 1
    elseif rand1 >= 1-p
      M[i] = M[i-1] + 1
    end
    maxVal = max(maxVal,M[i])
  end


  if maxVal >= b
    ST = max(S0*u^(M[N]) - S0*u^K,0)/(1+r)^N
  else
    ST = 0
  end
  return M, maxVal, ST
end

avgPrice = 0
avgPricesym = 0
numIter = 100000
for iter=1:numIter
  Msym, maxVal, STsym = symRandWalk(N,S0,u,K,b,r)
  #println(Msym)
  M, maxVal, ST = randWalk(N,S0,u,K,b,r)
  avgPrice += ST
  avgPricesym += STsym
end
avgPrice /= numIter
avgPricesym /= numIter
@printf("The average price for the Up and In barrier option is:           %f\n",avgPrice)
@printf("The average price for the Up and In barrier option is (sym):     %f\n",avgPricesym)



Msym, maxVal, STsym = symRandWalk(N,S0,u,K,b,r)
println(STsym)
plot(t,Msym[1:end])
