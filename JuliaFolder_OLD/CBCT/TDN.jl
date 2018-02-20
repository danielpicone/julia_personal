# This is a file which will solve the 4D-CBCT Problem

# This solves a quadratic shortest path problem on a TDN grid

B = 1 # Number of bins
P = 3 # Number of projections per bin
num_t = 4 # Number of discrete time steps
num_theta = 5 # Number of discrete angles
Θ_0 = 0 # Initial angle
Θ_end = 180 # Final angle
ω = 2
# Import the chest height at certain times
# For now just use a function
chest_func(t) = sin(10*t)

final_time = 240

# We will now create the TDN grid
function createTDNGrid(num_t::Int64,num_theta::Int64,Θ_0,Θ_end,B=1,ω=1,final_time=240)
    # Create the discrete angles
    discAngles = linspace(Θ_0,Θ_end,num_theta)
    discTime = linspace(0,final_time,num_t)
    maxHeight = maximum(chest_func.(discTime))
    minHeight = minimum(chest_func.(discTime))
    # Create the bin boundaries
    binBoundaries = linspace(minHeight,maxHeight,B+1)
    function getBin(chestHeight,binBoundaries)
        bin = 1
        while ((chestHeight > binBoundaries[bin]))
            bin+=1
            if bin>B
                break
            end
        end
        bin-=1
        return bin
    end
    # Total number of nodes
    num_nodes = 2*num_t * num_theta + 1
    nodes = Array{Float64}(num_nodes,4)
    # Create first row of TDN grid
    # This will have the form:
    # [NodeNumber DiscTimeInd DiscAngleInd Bin]
    nodes[1,1:3] = [1 0 1]
    for i=2:2*num_t+1
        # nodes[i,:] = [i discTime[Int64(floor(i/2))] discAngles[1] getBin(chest_func(discTime[Int64(floor(i/2))]),B) ]
        nodes[i,:] = [i Int64(floor(i/2)) 1 getBin(chest_func(discTime[Int64(floor(i/2))]),B) ]
    end
    for j=2:num_theta
        for i=j*num_t+2:2*j*num_t+1
            println(i)
            belowNode = ((i-2) % (2*num_t)) +2
            # nodes[i,:] = [i nodes[belowNode,2] discAngles[j] getBin(chest_func(nodes[belowNode,2]),B) ]
            nodes[i,:] = [i nodes[belowNode,2] j getBin(chest_func(nodes[belowNode,2]),B) ]
        end
    end
    for i=1:num_nodes
        # nodes[i,4] = chest_func(nodes[i,2])
    end
    # Now create the edges
    l = 1
    maxSpeedNodes = zeros(1)
    while (2*(l + (l-1)*ω*num_t)) <= num_nodes
        if l==1
            maxSpeedNodes[l] = 2*(l + (l-1)*ω*num_t)
        else
            maxSpeedNodes = [maxSpeedNodes; 2*(l + (l-1)*ω*num_t)]
        end
        l+=1
    end
    l = 1
    minSpeedNodes = zeros(1)
    while (num_nodes-1 - 2 * (l-1)*ω*num_t - 2) >= 0
        if l==1
            minSpeedNodes[l] = num_nodes-1- 2*(l-1)*ω*num_t-2*(l-1)
        else
            minSpeedNodes = [minSpeedNodes; num_nodes-1- 2*(l-1)*ω*num_t-2*(l-1)]
        end
        l+=1
    end
    sort!(minSpeedNodes)

    num_edges = 10
    fromToNode = Array{Float64}(1,8)
    fromToNode[1,:] = [1 0 discAngles[1] chest_func(0) 2 0 discAngles[1] chest_func(0)]
    # Create all edges which go from node i to node i+1
    for i=1:num_nodes
        if (i-1) % (2*num_t) == 0
        else
            if (i >= maxSpeedNodes[min(Int64(nodes[i,3]),length(maxSpeedNodes))]) && (i <= minSpeedNodes[min(Int64(nodes[i,3]),length(maxSpeedNodes))])
                # fromToNode = [fromToNode; nodes[i,:]' nodes[i+1,:]']
            else
            end
        end
    end
    # Create all edges going from node i to node i+2
    for i=1:num_nodes-2
        if i % 2 == 0
            # println(i)
            # println(maxSpeedNodes[min(Int64(nodes[i,3]),length(maxSpeedNodes))])
            # println(minSpeedNodes[min(Int64(nodes[i,3])+1,length(maxSpeedNodes))])
            if (i >= maxSpeedNodes[min(Int64(nodes[i,3]),length(maxSpeedNodes))]) && (i <= minSpeedNodes[min(Int64(nodes[i,3])+1,length(maxSpeedNodes))])
                fromToNode = [fromToNode; nodes[i,:]' nodes[i+2,:]']
            end
        else
        end
    end

    return nodes,discTime,discAngles,fromToNode,minSpeedNodes
end
