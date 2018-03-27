# This is a script which solves the cvrp via column generation
module col_gen_setup
using JuMP
using CPLEX
# include("setup.jl")
using cvrp_setup
using OffsetArrays
export create_subtour,get_subtour_value,create_initial_columns,get_subtour,solve_subprob,
    create_subprob,test_integrality,OG_solution_values

# First we create some columns through enumeration
function create_subtour(c,demands,capacity)
    # Choose a customer to go to first
    num_customers = length(indices(c)[1])-2
    nodes = 1:num_customers
    length_of_subtour = rand(2:num_customers)
    rand_indices = randperm(num_customers)
    i = 0
    total_demand = 0
    while total_demand <= capacity && (i < length(rand_indices) && i < length_of_subtour)
        i+=1
        total_demand += demands[rand_indices[i]]
    end
    if total_demand > capacity
        return nodes[rand_indices[1:i-1]]
    else
        return nodes[rand_indices[1:i]]
    end
end

function get_subtour_value(sub_tour,c)
    value = 0
    if sub_tour[1] == 0 || sub_tour[end]==0 || sub_tour[end]==length(indices(c))
    else
        sub_tour = [0; sub_tour; 0]
    end
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


function create_initial_columns(sub_tours,c,n)
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
    next_customer = 1
    while true
        next_customer = find(x[next_customer,:].==1)[1]
        subtour = [subtour; next_customer-1]
        if next_customer==n+2
            break
        end
    end
    return subtour[2:end-1]
end

function solve_subprob(u,σ,cvrp_problem::cvrp)
    n = length(indices(cvrp_problem.demands)[1])-2
    c = cvrp_problem.c
    demands = cvrp_problem.demands
    capacity = cvrp_problem.capacity
    subprob = Model(solver = CplexSolver(CPXPARAM_ScreenOutput=0))
    # subprob = Model(solver = CplexSolver(CPXPARAM_TimeLimit=200,CPXPARAM_Simplex_Limits_LowerObj=-10.0))

    @variable(subprob, x[i=0:n+1,j=0:n+1],Bin)
    @variable(subprob, demands[i] <= y[i=0:n+1] <= capacity)
    @objective(subprob, Min, sum( (c[i,j] - u[i+1])*x[i,j] for i=0:n+1,j=0:n+1 ) - σ)
    for i=1:n
        # @constraint(subprob, sum( x[i,j] for j=1:n if j!=i ) == 1)
        # @constraint(subprob, sum( x[i,j] for j=0:n+1 if j!=i ) == 1)
        # @constraint(subprob, sum( x[h,i] for h=0:n if h!=i ) - sum( x[i,j] for j=1:n+1 if j!=i ) == 0)
        @constraint(subprob, sum( x[h,i] for h=0:n+1 if h!=i) - sum( x[i,j] for j=0:n+1 if j!=i ) == 0)
    end

    @constraint(subprob, sum( x[0,j] for j=1:n ) == 1)
    # @constraint(subprob, sum( x[0,j] for j=1:n ) <= 3)
    @constraint(subprob, sum( x[i,n+1] for i=1:n ) == 1)
    # @constraint(subprob, sum( x[i,n+1] for i=1:n ) <= 3)

    for i=0:n+1
        for j=0:n+1
            @constraint(subprob, y[j] >= y[i] + demands[j]*x[i,j] - capacity*(1-x[i,j]))
        end
    end
    solve(subprob)
    x = getvalue(x)[:,:]

    x = round.(x)
    subprob_obj_val = getobjectivevalue(subprob)
    new_subtour = get_subtour(x,n)
    return subprob_obj_val,new_subtour

end


function create_subprob(u,σ,cvrp_problem::cvrp)
    n = length(indices(cvrp_problem.demands)[1])-2
    c = cvrp_problem.c
    demands = cvrp_problem.demands
    capacity = cvrp_problem.capacity
    u = OffsetArray(u,0:n+1)
    for i=0:n+1
        c[:,i] = c[:,i] - u
    end

    return c,demands,capacity
end

function perform_col_gen(cvrp_problem::cvrp,rmp)
    c,demands,capacity,n = _unpack(cvrp_problem)
    solve(rmp)
    u = getdual(all_customers_serviced)
    u = [0; u; 0]
    println("DOES")
    firstu = u
    σ = getdual(num_subtour_constraint)[1]
    firstσ = σ

    # subprob_obj_val,new_subtour = solve_subprob(u,σ,OGcvrp)
    # println(subprob_obj_val,new_subtour)
    new_subtour,subprob_obj_val = espprc_solve(c,u,demands,capacity)[1:2]
    λNewArray = Variable[]
    new_cost_array = []
    println("START COLUMN GENERATION")
    while subprob_obj_val < -10.0^-7
    # for l=1:1
        new_c,new_a_col = get_subtour_value(new_subtour,c)
        push!(new_a_col,1)
        @variable(rmp, 0 <= λ_new <= Inf,
            objective = new_c,
            inconstraints = [all_customers_serviced;num_subtour_constraint],
            coefficients = convert(Array{Float64,1},new_a_col))
        push!(λNewArray,λ_new)
        push!(new_cost_array,new_c)
        new_index = MathProgBase.numvar(rmp)
        index_dict[new_index] = new_subtour
        setname(λ_new,string("λ[",MathProgBase.numvar(rmp),"]"))
        redirect_stdout()
        solve(rmp)
        redirect_stdout(TT)
        u = getdual(all_customers_serviced)
        u = [0; u; 0]
        σ = getdual(num_subtour_constraint)[1]
        # tic()
        subprob_obj_val,new_subtour = solve_subprob(u,σ,OGcvrp)
        # toc()
        # println(new_subtour,subprob_obj_val)
        # tic()
        # println(espprc_solve(c,u,demands,capacity))
        # new_subtour,subprob_obj_val = espprc_solve(c,u,demands,capacity)[1:2]
        # toc()
        # println(new_subtour)
        println("The subprob objective value is: ",subprob_obj_val)
    end
    total_variables = append!(λ,λNewArray)
    solution = Array{Float64}(MathProgBase.numvar(rmp))
    println(λ)
    for i=1:MathProgBase.numvar(rmp)
        solution[i]=getvalue(λ[i])
    end

    return solution

end

function OG_solution_values(solution,index_dict,n)
    OGvariables = OffsetArray(zeros(n+2,n+2),0:n+1,0:n+1)
    for (index,value) in enumerate(solution)
        for i=1:length(index_dict[index])-1
            OGvariables[index_dict[index][i],index_dict[index][i+1]] += value
        end
    end
    return OGvariables

end

function test_integrality(solution,index_dict,n)
    most_fractional_value = 0
    most_fractional_index = NaN
    OGsolution = OG_solution_values(solution,index_dict,n)
    # for (index,value) in enumerate(OGsolution)
    for i=0:n+1
        for j=0:n+1
            value = OGsolution[i,j]
            if value - round(value) > most_fractional_value
                most_fractional_value = value
                most_fractional_index = (i,j)
                if most_fractional_value == 0.5
                    println(OGsolution[most_fractional_index...])
                    return most_fractional_index,most_fractional_value
                end
            end
        end
    end
    return most_fractional_index,most_fractional_value
end

# function set_variable_to_value_rmp(rmp,index_dict,index,value)
#     for (key,value) in enumerate(index_dict)
#         println("Value is: ",value[2])
#         println("index is: ",index)
#         if index[1] in value[2] && index[2] in value[2]
#             index1 = find(x->x==index[1],value[2])[1]
#             index2 = find(x->x==index[2],value[2])[1]
#             if abs(index1-index2)==1
#                 println("this occurs")
#                 @constraint(rmp,λ[key]==value[2])
#             end
#         end
#     end
#
#     return rmp
#
# end


end
