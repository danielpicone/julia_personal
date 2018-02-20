# This is a script to test the Levy flight

# First we draw two random variables from a normal distribution

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

β = 1

sigma_u = ( (gamma(1+β) * sin(π*β/2)) / (gamma((1+β)/2))*β*2^((β-1)/2)   )^(1/β)
sigma_v = 1

u = randn()*sigma_u
v = randn()*sigma_v

s = u/(abs(v)^(1/β))

numIter = 20

start = Array{Complex}(numIter,1)

start[1] = 0.0+0.0im
for i=1:numIter-1
    # choose a direction
    direction = (rand()-0.5)+(rand()-0.5)*im
    # Normalise the direction and multiply by the Levy flight length
    direction = direction./norm(direction)
    # Create the Levy flight
    u = randn()*sigma_u
    v = randn()*sigma_v

    # s = u/(abs(v)^(1/β))
    s = rand()^(-2.5)
    move = abs(s)*direction
    start[i+1] = start[i]+move
end

plot(real.(start),imag.(start))
