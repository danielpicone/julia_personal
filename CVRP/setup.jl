# This is a module to setup the problem
module cvrp_setup

export create_inputs

using OffsetArrays

srand(100)
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
    capacity = 30
    return c,demands,capacity
end



end
