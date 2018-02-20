# This is a script to create a functor for a neural network
# It takes in as an argument the points and then the weights as the argument of a function

function constraints(x; bounds = [-Inf*ones(dim) Inf*ones(dim)])
    # First the equality constraints
    # These are in the form g(x) = 0
    numeqconstraints = 0
    eq = Array{Float64}(numeqconstraints)
    if numeqconstraints!=0
        eq[1] = x[1]^2-1
        # eq[2] = sin(x[1])
    end
    # Now the inequality constraints
    # These are in the form g(x) ≤ 0
    numineqconstraints = 0
    ineq = Array{Float64}(numineqconstraints)
    if numineqconstraints!=0
        ineq[1] = 0.25 - x[2]
        # ineq[1] = 1 - (x[2]^2 * x[3])/(71785*x[1]^4)
        # ineq[2] = (4*x[2]^2 -x[1]*x[2])/(12566*(x[1]^3 * x[2] - x[1]^4)) + 1/(5108*x[1]^2) - 1
        # ineq[3] = 1 - (140.45*x[1])/(x[2]^2 * x[3])
        # ineq[4] = (x[1] + x[2])/1.5 - 1
    end
    lbbound = Array{Float64}(dim)
    ubbound = Array{Float64}(dim)
    for k=1:dim
        lbbound[k] = bounds[k,1] - x[k]
        ubbound[k] = -bounds[k,2] + x[k]
    end

    return 10000000*(sum(eq.^2)+sum(max.(ineq,0)) + sum(max.(lbbound,0) + max.(ubbound,0)))
end


struct neural_network{R,R_output}
    # data::Array{R,2}
    data::Array{R,2}
    output::Array{R_output,2}
end

function (nn::neural_network)(W)
    println("nn was called")
    # First we decide on how many neurons the network has in each layer
    final_layer = 2
    layers = [size(nn.data,2),2,final_layer]
    biggest_layer = maximum(layers)
    num_layers = length(layers)
    # Now create a for loop which loops over all the points
    num_points = size(nn.data,1)
    # Reshape W to make it a 3D matrix
    if length(W) != biggest_layer*(biggest_layer+1)*(num_layers-1)
        println("The length of W needs to be ", biggest_layer*(biggest_layer+1)*(num_layers-1))
        return
    end
    # println("Does this show??")
    # println(W)
    W = reshape(W, (biggest_layer,biggest_layer+1,num_layers-1))
    # Make sure all entries which go to an empty node are zero
    for layer=1:num_layers-1
        # Need to add biggest_layer-layers[layer]
        W[(layers[layer]+1):end,:,layer] .= 0
    end
    # println(W)
    # Also choose the eval function
    σ(x) = (1-exp(-x))/(1+exp(-x))
    # Define the loss function
    L(x,y) = norm(x-y)
    # loss = Vector{size(nn.data,1)}
    loss = Inf*ones(size(nn.data,1))
    # println(loss)
    for i=1:num_points
        # Initialise the first layer
        a0 = nn.data[i,:]
        z = W[:,end,1]+W[:,1:(end-1),1]*a0
        # Now create a for loop which loops over all layers
        for l=1:num_layers-1
            a = σ.(z)
            z = W[:,end,l] + W[:,1:(end-1),l]*a
            # println(i,l)
        end
        output = z

        loss[i] = L(output[size(nn.output,2)], nn.output[i])

    end
    return sum(loss)
end


# Create a simple classification problem to solve
# Generate uniform points and classify them as 1 if they are inside the unit circle
# Classify them as 0 otherwise
points = 20
fake_data = 2*rand(points,2)-1
fake_data_output = zeros(20,2)
for i=1:points
    if norm(fake_data[i,:]) < 1
        fake_data_output[i,1] = 1
    else
        fake_data_output[i,2] = 1
    end
end

nn = neural_network(fake_data, fake_data_output)
total_loss = nn(rand(12))
println("The total loss is: ",total_loss)

include("FA.jl")
println("Trying nn now")
fa(nn, bounds = bounds)
