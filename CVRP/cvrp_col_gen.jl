# This is a script which solves the cvrp via column generation

using JuMP
using CPLEX
include("setup.jl")
include("col_gen_setup.jl")
include("espprc.jl")
using cvrp_setup
using col_gen_setup
using GLPKMathProgInterface
using PyPlot

TT = STDOUT
tic()
n = 10
K = 10
srand(13)
OGcvrp = create_inputs(n)
c,demands,capacity,n = _unpack(OGcvrp)
# First we create some columns through enumeration
sub_tours = Array{Array{Int64}}(n)
# Create the subtours which only have one stop first
for i=1:n
    sub_tours[i] = [i]
end

# sub_tours = [sub_tours; [[1,2]]]

# Now create random subtours
# for i=1:10000
#     sub_tours = [sub_tours; [create_subtour(c,demands,capacity)]]
# end
# Delete duplicate tours
sub_tours = unique(sub_tours)

cost_array,customer_matrix,index_dict = create_initial_columns(sub_tours,c,n)

rmp = Model(solver = CplexSolver(CPXPARAM_ScreenOutput=0,CPXPARAM_Preprocessing_Dual=-1))
# rmp = Model(solver = GLPKSolverLP(method=:Exact,msg_lev=0))

num_var = length(cost_array)

@variable(rmp, λ[1:num_var] >= 0)

@objective(rmp, Min,sum( cost_array[r]*λ[r] for r=1:num_var ))

constrrefs = Vector{ConstraintRef}(0)
@constraintref(all_customers_serviced[1:size(customer_matrix,1)])

@constraint(rmp, all_customers_serviced[i=1:size(customer_matrix,1)], sum( customer_matrix[i,r]*λ[r] for r=1:num_var ) == 1)
@constraintref(num_subtour_constraint[1:1])
@constraint(rmp, num_subtour_constraint[1:1],sum( λ[r] for r=1:num_var ) <= K)

function perform_col_gen(cvrp_problem::cvrp,rmp)
    c,demands,capacity,n = _unpack(cvrp_problem)
    solve(rmp)
    u = getdual(all_customers_serviced)
    u = [0; u; 0]
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



lambda_solution = perform_col_gen(OGcvrp,rmp)
println(test_integrality(lambda_solution,index_dict,n))

function set_variable_to_value_rmp(rmp,index_dict,index,value)
    for (key,dict_value) in enumerate(index_dict)
        println("Value is: ",dict_value[2])
        println("index is: ",index)
        if index[1] in dict_value[2] && index[2] in dict_value[2]
            index1 = find(x->x==index[1],dict_value[2])[1]
            index2 = find(x->x==index[2],dict_value[2])[1]
            if abs(index1-index2)==1
                println("this occurs")
                println(λ[key])
                println(value)
                @constraint(rmp,λ[key]==value)
            end
        end
    end

    return rmp

end
rmp1 = rmp
rmp1 = set_variable_to_value_rmp(rmp1,index_dict,(3,7),0)
rmp2 = set_variable_to_value_rmp(rmp,index_dict,(3,7),1)
solve(rmp1)
solve(rmp2)

# solve(rmp)
# println(getvalue(λ))
# u = getdual(all_customers_serviced)
# u = [0; u; 0]
# firstu = u
# σ = getdual(num_subtour_constraint)[1]
# firstσ = σ
#
# # subprob_obj_val,new_subtour = solve_subprob(u,σ,OGcvrp)
# # println(subprob_obj_val,new_subtour)
# new_subtour,subprob_obj_val = espprc_solve(c,u,demands,capacity)[1:2]
# λNewArray = Variable[]
# new_cost_array = []
# println("START COLUMN GENERATION")
# while subprob_obj_val < -10.0^-7
# # for l=1:1
#     new_c,new_a_col = get_subtour_value(new_subtour,c)
#     push!(new_a_col,1)
#     @variable(rmp, 0 <= λ_new <= Inf,
#         objective = new_c,
#         inconstraints = [all_customers_serviced;num_subtour_constraint],
#         coefficients = convert(Array{Float64,1},new_a_col))
#     push!(λNewArray,λ_new)
#     push!(new_cost_array,new_c)
#     new_index = MathProgBase.numvar(rmp)
#     index_dict[new_index] = new_subtour
#     setname(λ_new,string("λ[",MathProgBase.numvar(rmp),"]"))
#     redirect_stdout()
#     solve(rmp)
#     redirect_stdout(TT)
#     u = getdual(all_customers_serviced)
#     u = [0; u; 0]
#     σ = getdual(num_subtour_constraint)[1]
#     # tic()
#     subprob_obj_val,new_subtour = solve_subprob(u,σ,OGcvrp)
#     # toc()
#     # println(new_subtour,subprob_obj_val)
#     tic()
#     # println(espprc_solve(c,u,demands,capacity))
#     # new_subtour,subprob_obj_val = espprc_solve(c,u,demands,capacity)[1:2]
#     toc()
#     println(λ_new)
#     println(getvalue(λ_new))
#     # println(new_subtour)
#     println("The subprob objective value is: ",subprob_obj_val)
# end


function graph_cvrp(subtours,cvrp_problem::cvrp)
    println(cvrp_problem.customers)
    plot(cvrp_problem.customers[:,1],cvrp_problem.customers[:,2],"ro")
    # return cvrp_problem.customers
end
