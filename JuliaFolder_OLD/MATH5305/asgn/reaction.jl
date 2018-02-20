# This is the solution to Assignment 2 Question 2
#     TITLE: reaction.jl
#     By: Daniel Picone

using PyPlot

function q2(ϵ,u0,T,N,η,uexact)
  t = linspace(0,T,N+1)
  Δt = T/N
  U = zeros(N+1)
  U[1] = u0
  U[2] = u0*(1+Δt * ϵ^(-2)*(1-u0^2))
  #U[3:N+1] = Float64[ (U[n-1] + ϵ^(-2)*Δt*( 3/2 U[n-1](1-U[n-1]^2) -1/2 * U[n-2]*(1-U[n-2]^2) )) for n=2:N ]
  for n=3:N+1
    U[n] =  U[n-1] + ϵ^(-2)*Δt*( 3/2 * U[n-1]*(1-U[n-1]^2) -1/2 * U[n-2]*(1-U[n-2]^2) )
  end


  return t,U
end

ϵ = 0.01
η = 10.0^(-4)
u0 = 0.1
T = ϵ^2 * log(1/(u0*sqrt(2*η)))
N = 50

A = u0/(sqrt(1-u0^2))
uexact(t) = A.*exp(ϵ^(-2).*t)./(sqrt(1+A^2 .* exp(2.*ϵ^(-2).*t)))
taxis = linspace(0,T,1000)

t,U = q2(ϵ,u0,T,N,η,uexact)
# Plot the FD soltuion and the exact solution
figure(1)
plot(t,U,"ro",linewidth=5)
plot(taxis,uexact(taxis))
title(L"$N = 50$")
xlabel(L"$t$")
grid(linestyle = ":")
legend((L"$U^n$",L"$u(t)$"))

# Give the maximum error and convergence rate of the FD solution
numiter = 6
maxerr = zeros(numiter)
rn = zeros(numiter)
@printf("  N     error         rate\n\n")
for j=1:numiter
  N = 20*2^(j-1)
  t,U = q2(ϵ,u0,T,N,η,uexact)
  err = zeros(N+1)
  for i=1:N+1
    err[i] = abs(U[i] - uexact(t[i]))
    maxerr[j] = max(maxerr[j],err[i])
  end
  if (j >= 2)
    rn[j] = log2(maxerr[j-1]/maxerr[j])
  end
  if j==1
    @printf("%3d     %1.3e\n",N,maxerr[j])
  else
    @printf("%3d     %1.3e     %1.4f\n",N,maxerr[j],rn[j])
  end
end
