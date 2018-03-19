# Solving the traveling salesman problem (TSP)

# include("TSP_module.jl")

srand(100)

using JuMP
using CPLEX
# using TSP_module

num_of_cities = 400

# First create random points

cities = rand(1:100,num_of_cities,2)

# Now create the distance matrix

D = Array{Float64}(num_of_cities,num_of_cities)

for i=1:num_of_cities
    for j=i:num_of_cities
        D[i,j] = sqrt((cities[i,1] - cities[j,1])^2 + (cities[i,2] - cities[j,2])^2)
        # D[j,i] = D[i,j]
    end
end

D = Symmetric(D)

# Now create the IP to solve
TSP = Model(solver = CplexSolver())

# The variables are the edges
#@variable(TSP, x[i=1:num_of_cities,j=i+1:num_of_cities],Bin)
@variable(TSP, x[i=1:num_of_cities,j=1:num_of_cities],Bin)
@variable(TSP, u[i=1:num_of_cities],Int)
# The objective is to minimise the distance of the path
@objective(TSP, Min, sum( D[i,j]*x[i,j] for i=1:num_of_cities,j=i:num_of_cities ))
# @objective(TSP, Min, 1)
# Add the constraints
# Path entering one node must be 1

for i=1:num_of_cities
    @constraint(TSP, sum(x[i,j] for j=1:num_of_cities) == 2)
    @constraint(TSP, x[i,i] == 0)
    for j=i+1:num_of_cities
        @constraint(TSP, x[i,j] == x[j,i])
    end
end

solve(TSP)
println(getvalue(x))
println(getvalue(x[1,1]))
valx = getvalue(x)
objval = 0
for i = 1:num_of_cities
    for j=i:num_of_cities
        objval += D[i,j]*valx[i,j]
    end
end
println(objval)
println(getobjectivevalue(TSP))

# println("Next Iteration of TSP")

# @constraint(TSP, x[1,2] + x[2,3] >= 1)
# solve(TSP)

# println(getvalue(x))
# println(getobjectivevalue(TSP))
