# SDE Solver

# We want to numerically solve an SDE with form
# dX_t = a(X_t)dt + b(X_t)dW_t
# with initial condition X0 = x0

using PyPlot

N = 1000
T = 15
t = linspace(0,T,N)
x0 = 1

X = zeros(N)
X[1] = x0

Y = zeros(N)
Y[1] = x0

# This is an example using the O-U Process
μ = 1.5     # This is the long term average of the process
Θ = 0.1     # Seems to determine how strongly the process is "pulled" back to μ
σ = 0.5
a(t,x) = Θ*(μ-x)*x
b(t,x) = σ*x

Γ = zeros(N)
Δt = T/N

for i=1:N-1
  ΔW = sqrt(Δt)*randn()
  X[i+1] = X[i] + a(t[i],X[i])*Δt + b(t[i],X[i])*ΔW
  Γ[i]   = Y[i] + a(t[i],Y[i])*Δt + b(t[i],Y[i])*Δt^(0.5)
  Y[i+1] = Y[i] + a(t[i],Y[i])*Δt + b(t[i],X[i])*ΔW + 0.5*(b(t[i],Γ[i]) - b(t[i],Y[i]))*(ΔW^2 - Δt)*Δt^(-0.5)
end

plot(t,X,t,Y,t,μ*ones(N),"k")
#plot(t,Y,t,μ*ones(N),"k")
