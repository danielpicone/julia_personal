# This is a script to create a functor for a neural network
# It takes in as an argument the points and then the weights as the argument of a function

srand(100)
# global dim = 36

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


# struct neural_network{R,R_output,layers}
#     # data::Array{R,2}
#     data::typeof(R)
#     output::Array{R_output,2}
#     layers::Array{layers,1}
# end

struct neural_network
    data
    output
    layers
end

# struct neural_network_eval{R,layers}
#     data::Array{R,2}
#     layers::Array{layers,1}
# end

struct neural_network_eval
    data
    layers
end

function (nn::neural_network)(W)

    function cross_entropy(x,y)
        # return -y*log(x)-(1-y)*log(1-x)
        return -sum(x.*log.(y) + (1.-x).*log.(1.-y))
    end


    # println("nn was called")
    # First we decide on how many neurons the network has in each layer
    layers = nn.layers
    biggest_layer = maximum(layers)
    num_layers = length(layers)
    # Now create a for loop which loops over all the points
    num_points = size(nn.data,1)
    # Reshape W to make it a 3D matrix
    if length(W) != biggest_layer*(biggest_layer+1)*(num_layers-1)
        println("The length of W needs to be ", biggest_layer*(biggest_layer+1)*(num_layers-1))
        return
    end
    W = reshape(W, (biggest_layer,biggest_layer+1,num_layers-1))
    # Make sure all entries which go to an empty node are zero
    for layer=1:num_layers-1
        # Need to add biggest_layer-layers[layer]
        one = [ones(layers[layer+1],layers[layer]) zeros(layers[layer+1],biggest_layer-layers[layer])]
        zero = zeros(biggest_layer-layers[layer+1],biggest_layer)
        last = [ones(layers[layer+1]); zeros(biggest_layer-layers[layer+1])]
        c = [one; zero]
        c = [c last]
        W[:,:,layer] = c.*W[:,:,layer]
    end
    # Also choose the eval function
    σ(x) = (1-exp(-x))/(1+exp(-x))
    # Define the loss function
    L(x,y) = norm(x-y)
    # loss = Vector{size(nn.data,1)}
    loss = Inf*ones(size(nn.data,1))
    output = zeros(biggest_layer)
    for point=1:num_points
        # Initialise the first layer
        a0 = nn.data[point,:]
        z = W[:,end,1]+W[:,1:(end-1),1]*a0
        # Now create a for loop which loops over all layers
        for l=1:num_layers-1
            a = σ.(z)
            z = W[:,end,l] + W[:,1:(end-1),l]*a
        end
        # output .= z
        output .= exp.(z)/(sum(exp.(z)))
        loss[point] = L(output, nn.output[point,:])
        # loss[point] = cross_entropy(output, nn.output[point,:])
    end

    return 1000*sum(loss)
end

function (nn_eval::neural_network_eval)(W)
    # println("nn_eval was called")
    # First we decide on how many neurons the network has in each layer
    layers = nn_eval.layers
    biggest_layer = maximum(layers)
    num_layers = length(layers)
    # Now create a for loop which loops over all the points
    num_points = size(nn_eval.data,1)
    # Reshape W to make it a 3D matrix
    if length(W) != biggest_layer*(biggest_layer+1)*(num_layers-1)
        println("The length of W needs to be ", biggest_layer*(biggest_layer+1)*(num_layers-1))
        return
    end
    W = reshape(W, (biggest_layer,biggest_layer+1,num_layers-1))
    # Also choose the eval function
    σ(x) = (1-exp(-x))/(1+exp(-x))
    # Define the loss function
    L(x,y) = norm(x-y)
    # loss = Vector{size(nn.data,1)}
    loss = Inf*ones(size(nn.data,1))

    output = zeros(num_points, biggest_layer)
    for i=1:num_points
        # Initialise the first layer
        a0 = nn_eval.data[i,:]
        z = W[:,end,1]+W[:,1:(end-1),1]*a0
        # Now create a for loop which loops over all layers
        for l=1:num_layers-1
            a = σ.(z)
            z = W[:,end,l] + W[:,1:(end-1),l]*a
            # println(i,l)
        end
        # output[i,:] = z
        output[i,:] .= exp.(z)/(sum(exp.(z)))

        # loss[i] = L(output[size(nn.output,2)], nn.output[i])

    end
    return output
end

# Create a simple classification problem to solve
# Generate uniform points and classify them as 1 if they are inside the unit circle
# Classify them as 0 otherwise


input_layer = 1

layers = [input_layer,30,30,1]


points = 100
biggest_layer = maximum(layers)
num_layers = length(layers)
global dim = biggest_layer*(biggest_layer+1)*(num_layers-1)
# training_data = π*rand(points,input_layer)
# training_data = [training_data zeros(points,biggest_layer-input_layer)]
# training_data_output = zeros(points,biggest_layer)
# for i=1:points
#     training_data_output[i,1] = sin(training_data[i,1])+0.2*randn()
# end

training_data = 3*(rand(points,input_layer)-0.5)
training_data = [training_data zeros(points,biggest_layer-input_layer)]
training_data_output = zeros(points,biggest_layer)
for i=1:points
    # if norm(training_data[i,1]) < 0.8
    #     training_data_output[i,1] = 1.0
    # else
    #     training_data_output[i,2] = 1.0
    # end
    training_data_output[i,1] = sin(training_data[i,1])
end

nn = neural_network(training_data, training_data_output, layers)
println("Training now")
total_loss = nn(rand(dim))
println("The total loss is: ",total_loss)

include("FA.jl")
println("Trying nn now")
x, values, bestpos, bestval = fa(nn, α0 = 3.0, maxiter = 1000, popnum = 70, bounds = bounds, verbose = "yes", maxcounter = 10000000);
weights = bestpos
println("Training is complete")
# test_data = π*rand(points,input_layer)
test_data = 3*rand(points,input_layer)-1
test_data = [test_data zeros(points,biggest_layer-input_layer)]
test_data_output = zeros(points,biggest_layer)
for i=1:points
    # if test_data[i,1] < 0.7
    #     test_data_output[i,1] = 1.0
    # else
    #     test_data_output[i,2] = 1.0
    # end
    test_data_output[i,1] = sin(test_data[i,1])
end
best_network = neural_network_eval(test_data, layers)
final_output = best_network(weights)

println(size(final_output))
println(size(test_data_output))

test_loss = norm(final_output[:,1]-test_data_output[:,1])

println("The test loss is: ", test_loss)

weights = reshape(weights, (biggest_layer,biggest_layer+1,num_layers-1))

correctly_classified = 0
for i=1:points
    if norm(final_output[i,1] - test_data_output[i,1]) < 0.1
        correctly_classified += 1
    end
end
println("The number of correctly classified points is ", correctly_classified)

# x = linspace(0,π,100)
x = linspace(-1.5,1.5,100)
x = collect(x)
x = [x zeros(100,biggest_layer-input_layer)]
network = neural_network_eval(x,layers)
network(weights)
# plot(x,network(weights)[:,1],x,sin.(x))
plot(x,network(weights)[:,1])
plot(training_data[:,1],training_data_output[:,1],".")
