# Dynamic programming algorithm for the ESPRP

using OffsetArrays

function espprc_solve(c,u,demands,capacity)
    start_node = 0
    end_node = n+1
    node_list = collect(0:n+1)
    u = OffsetArray(u,0:n+1)
    adjusted_cost = copy(c)
    for i=0:n+1
        adjusted_cost[:,i] = c[:,i] - u
    end
    # c = c'
    paths = []

    # Create first n paths
    for i=1:n
        new_path = ([0,i],adjusted_cost[0,i],demands[i])
        push!(paths,new_path)
    end
    deleteat!(node_list,1)

    complete_paths = []
    total_paths = paths
    paths_added = paths
    iter = 1
    while length(paths_added)!=0
        paths = paths_added
        paths_added=[]

        for new_node in node_list
            for (index,path) in enumerate(paths)
                if !(new_node in path[1])
                    new_path = [path[1];new_node]
                    new_cost = path[2]+adjusted_cost[path[1][end],new_node]
                    new_demand = path[3]+demands[new_node]
                    if new_demand <= capacity
                        # Add new_path to list of paths
                        if new_node==n+1
                            push!(complete_paths,(new_path,new_cost,new_demand))
                            push!(total_paths,(new_path,new_cost,new_demand))
                        else
                            # Check for dominance
                            dominated=0
                            push!(paths_added,(new_path,new_cost,new_demand))
                            push!(total_paths,(new_path,new_cost,new_demand))
                        end
                    end
                end

            end
            # filter!(e->e != path,paths)
            if length(paths)==0
                break
            end
            # println("Length of paths_added is: ",length(paths_added))

        end
        deleteat!(node_list,1)
        iter+=1

    end

    best_path = []
    best_value = Inf

    for path in unique(complete_paths)
        # if [0,1,10,5,3,n+1] == path[1] || [0,3,5,10,1,n+1] == path[1]
        if path[2] < best_value
            best_path = path
            best_value = path[2]
        end
    end
    return best_path
end
