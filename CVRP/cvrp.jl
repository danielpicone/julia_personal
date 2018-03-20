# This is a script which solves the cvrp

include("setup.jl")
using cvrp_setup
using JuMP
using CPLEX
using LightGraphs
using PlotRecipes

n = 27
K = 20
srand(10)
c,demands,capacity = create_inputs(n)

cvrp = Model(solver = CplexSolver())

@variable(cvrp, x[i=0:n+1,j=0:n+1], Bin)
@variable(cvrp, demands[i] <= y[i=0:n+1] <= capacity)

# First create the objective value
@objective(cvrp, Min, sum( c[i,j]*x[i,j] for i=0:n+1,j=i:n+1 ))

# Now create the flow constraints
for i=1:n
    @constraint(cvrp, sum( x[i,j] for j=1:n+1 ) == 1)
end
for h=1:n
    @constraint(cvrp, sum( x[i,h] for i=0:n ) - sum( x[h,j] for j=1:n+1 ) == 0)
end
@constraint(cvrp, sum( x[0,j] for j=1:n ) <= K)

for i=0:n+1
    for j=0:n+1
        @constraint(cvrp, y[j] >= y[i] + demands[j]*x[i,j] - capacity*(1-x[i,j]))
    end
    @constraint(cvrp, x[i,i]==0)
end
solve(cvrp)

# xmat = getvalue(x)
# adj_graph = Array{Int64}(n+2,n+2)
# for i=0:n+1
#     for j=0:n+1
#         adj_graph[i+1,j+1] = xmat[i,j]
#     end
# end
#
# graphplot(sparse(adj_graph),names=0:n)
