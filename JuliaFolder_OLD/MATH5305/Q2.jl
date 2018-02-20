# v2.0 of Q2

# First we shall define a function which solves a bvp given general functions and domain

using PyPlot

function solve(a,b,c,f,γ0,γL,L,P)

  Δx = L/P
  x = linspace(0,L,P+1)
  axd2 = Float64[ a(x[p+1])/(Δx^2) for p=1:P-1]
  A = Tridiagonal(-axd2[2:P-1],2*axd2[1:P-1],-axd2[1:P-2])
  bxd = Float64[ b(x[p+1])/(2*Δx) for p=1:P-1]
  B = Tridiagonal(-bxd[2:P-1],zeros(P-1),bxd[1:P-2])
  C = Diagonal(Float64[ c(x[p+1]) for p=1:P-1 ])

  rhs = Float64[ f(x[p+1]) for p=1:P-1]
  rhs[1] += (axd2[1] + bxd[1])*γ0
  rhs[P-1] += (axd2[P-1] - bxd[P-1])*γL
  U = zeros(P+1)
  U[1] = γ0
  U[P+1] = γL
  U[2:P] = (A+B+C)\rhs
  return U, x
end

a(x) = 1
b(x) = 2
c(x) = -1
f(x) = sin(x)
γ0 = 1
γL = -1
L = 2
P = 10000
start = time()
U , x = solve(a,b,c,f,γ0,γL,L,P)
stop = time()
#println(U)

println(-start+stop)

#figure(1)
plot(x[1:Int64(P/4)],U[1:Int64(P/4)],"bD",x[Int64(P/4):Int64(3P/4)],U[Int64(P/4):Int64(3P/4)],"h",x[Int64(3P/4):P],U[Int64(3P/4):P],"H")
#@printf("%3.3lf",x[3])
