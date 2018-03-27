# This is a module to setup the problem
module cvrp_setup
using JuMP
using CPLEX
export create_inputs, cvrp, _unpack, solve_OG

using OffsetArrays

# srand(100)
struct cvrp
    c
    demands
    capacity
    customers
end


function create_inputs(n)

    # First create random points
    customers = rand(-100:100,n+2,2)
    customers[1,:] = [0,0]
    customers[end,:] = customers[1,:]

    # Now create the cost matrix

    c = OffsetArray{Float64}(0:n+1,0:n+1)

    for i=0:n+1
        for j=i:n+1
            c[i,j] = sqrt((customers[i+1,1] - customers[j+1,1])^2 + (customers[i+1,2] - customers[j+1,2])^2)
            c[j,i] = c[i,j]
        end
    end

    demands = rand(1:10,n+2)
    demands[1] = 0
    demands[end] = 0
    demands = OffsetArray(demands,0:n+1)
    capacity = 10
    return cvrp(c,demands,capacity,customers)
end

function solve_OG(cvrp_problem::cvrp)
    c,demands,capacity,n = _unpack(cvrp_problem)
    K = n
    cvrp = Model(solver = CplexSolver())

    @variable(cvrp, x[i=0:n+1,j=0:n+1], Bin)
    @variable(cvrp, demands[i] <= y[i=0:n+1] <= capacity)

    # First create the objective value
    @objective(cvrp, Min, sum( c[i,j]*x[i,j] for i=0:n+1,j=0:n+1 ))

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
    return getobjectivevalue(cvrp),getvalue(x)
end


function _unpack(cvrp_problem)
    c = cvrp_problem.c
    demands = cvrp_problem.demands
    capacity = cvrp_problem.capacity
    n = length(indices(cvrp_problem.demands)[1])-2
    return c,demands,capacity,n
end

end
