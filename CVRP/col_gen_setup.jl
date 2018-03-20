# This is a script which solves the cvrp via column generation
module col_gen_setup
# include("setup.jl")
# using cvrp_setup
export create_subtour,get_subtour_value,create_columns,get_subtour

# First we create some columns through enumeration
function create_subtour(c,demands,capacity)
    # Choose a customer to go to first
    num_customers = length(indices(c)[1])-2
    nodes = 1:num_customers
    length_of_subtour = rand(2:num_customers)
    rand_indices = randperm(num_customers)
    i = 1
    total_demand = 0
    while total_demand <= capacity && (i < length(rand_indices) && i < length_of_subtour)
        total_demand += demands[rand_indices[i]]
        i+=1
    end
    return nodes[rand_indices[1:i]]
end

function get_subtour_value(sub_tour,c)
    value = 0
    sub_tour = [0; sub_tour; 0]
    for i=1:length(sub_tour)-1
        value += c[sub_tour[i],sub_tour[i+1]]
    end
    a = Array{Int64}(length(indices(c)[1])-2)
    for i=1:length(a)
        if i in sub_tour
            a[i] = 1
        else
            a[i] = 0
        end
    end
    return value,a
end


function create_columns(sub_tours,c,n)
    # Create dictionary of lambda variables
    index_dict = Dict{Int64,Array{Int64}}()
    # Create cost array and a matrix
    cost_array = Array{Float64}(length(sub_tours))
    customer_matrix = Array{Int64}(n,length(sub_tours))
    for (index,sub_tour) in enumerate(sub_tours)
        value,a = get_subtour_value(sub_tour,c)
        index_dict[index] = sub_tour
        cost_array[index] = value
        customer_matrix[:,index] = a
    end
    return cost_array,customer_matrix,index_dict
end

function get_subtour(x,n)
    subtour = Array{Int64}(1)
    next_customer = 0
    while true
        next_customer = find(x[next_customer,:].==1)[1]-1
        subtour = [subtour; next_customer]
        if next_customer==n+1
            break
        end
    end
    return subtour[2:end-1]
end



end
