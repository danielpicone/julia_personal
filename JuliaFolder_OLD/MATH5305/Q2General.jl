# Q2 General practice

using PyPlot

function q2(a,b,c,f,γ0,γL,L,P)

  Δx = L/P
  x = linspace(0,L,P+1)
  adx2 = Float64[ a(x[p+1])/(Δx^2) for p=1:P-1 ]
  A = Tridiagonal(-adx2[2:P-1],2*adx2[1:P-1],-adx2[1:P-2])
  bdx = Float64[ b(x[p+1])/(2*Δx) for p=1:P-1 ]
  B = Tridiagonal(-bdx[2:P-1],zeros(P-1),bdx[1:P-2])
  C = Diagonal(Float64[ c(x[p+1]) for p=1:P-1 ])
  rhs = Float64[ f(x[p+1]) for p=1:P-1 ]
  rhs[1] += (adx2[1] + bdx[1])*γ0
  rhs[P-1] += (adx2[P-1] - bdx[P-1])*γL
  U = zeros(P+1)
  U[1] = γ0
  U[2:P] = (A+B+C)\rhs
  U[P+1] = γL

  return U,x
end

a(x) = 1
b(x) = 2
c(x) = -1
f(x) = 1
u(x) = -1 + 2.*exp(x) -x.*exp(x)
γ0 = 1
γL = -1
P = 6
L = 2

U,x = q2(a,b,c,f,γ0,γL,L,P)
xaxis = linspace(0,L,1000)
plot(x,U,"r",xaxis,u(xaxis))
@printf("x value    f.d. value     function value  error\n\n")
for i=1:P+1
  @printf("%2.2f       %2.5f        %2.5f         %2.10f\n",x[i],U[i],u(x[i]),U[i]-u(x[i]))
end
