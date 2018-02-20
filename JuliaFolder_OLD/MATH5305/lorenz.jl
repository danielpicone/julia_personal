# This is a script to solve the Lorenz Equation using the Euler Method

# Let the Lorenz Equations be in the form
# dx/dt = f(t,x)
# dy/dt = g(t,x)
# dz/dt = h(t,x)

using PyPlot
i = 1
N = 1_000_00
σ = 10
ρ = 56
β = 8/3

x0 = 1
y0 = 1
z0 = 1

T = 20
t = linspace(0,T,N)
Δt = T/N

X = zeros(N,3)
X[1,1] = x0
X[1,2] = y0
X[1,3] = z0

f(t,x) = σ*x[i,2] - σ*x[i,1]
g(t,x) = ρ*X[i,1] - X[i,1]*X[i,3] - X[i,2]
h(t,x) = X[i,1]*X[i,2] - β*X[i,3]

for i=1:N-1
  X[i+1,1] = X[i,1] + f(t,X)*Δt
  X[i+1,2] = X[i,2] + g(t,X)*Δt
  X[i+1,3] = X[i,3] + h(t,X)*Δt
end


#=
f(t,x,y,z) = σ*y - σ*x
g(t,x,y,z) = ρ*x - x*z - y
h(t,x,y,z) = x*y - β*z



for i=1:N-1
  k1x = f(t[i],X[i,1],X[i,2],X[i,3])
  k2x = f(t[i] + Δt/2,X[i,1]+Δt*k1x/2,X[i,2],X[i,3])
  k3x = f(t[i] + Δt/2,X[i,1]+Δt*k2x/2,X[i,2],X[i,3])
  k4x = f(t[i] + Δt,X[i,1]+Δt*k3x,X[i,2],X[i,3])

  k1y = f(t[i],X[i,1],X[i,2],X[i,3])
  k2y = f(t[i] + Δt/2,X[i,1],X[i,2]+Δt*k1y/2,X[i,3])
  k3y = f(t[i] + Δt/2,X[i,1],X[i,2]+Δt*k2y/2,X[i,3])
  k4y = f(t[i] + Δt,X[i,1],X[i,2]+Δt*k3y,X[i,3])

  k1z = f(t[i],X[i,1],X[i,2],X[i,3])
  k2z = f(t[i] + Δt/2,X[i,1],X[i,2],X[i,3]+Δt*k1z/2)
  k3z = f(t[i] + Δt/2,X[i,1],X[i,2],X[i,3]+Δt*k2z/2)
  k4z = f(t[i] + Δt,X[i,1],X[i,2],X[i,3]+Δt*k3z)

  X[i+1,1] = X[i,1] + Δt/6 * (k1x + 2*k2x + 2*k3x + k4x)
  X[i+1,2] = X[i,2] + Δt/6 * (k1y + 2*k2y + 2*k3y + k4y)
  X[i+1,3] = X[i,3] + Δt/6 * (k1z + 2*k2z + 2*k3z + k4z)
end
=#
figure(1)
mesh(X[:,1],X[:,2],X[:,3])
figure(2)
plot(t,X[:,1],t,X[:,3])
