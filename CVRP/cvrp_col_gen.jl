# This is a script which solves the cvrp via column generation

include("setup.jl")
include("col_gen_setup.jl")
using cvrp_setup
using col_gen_setup
using JuMP
using CPLEX
using GLPKMathProgInterface

n = 3
K = 3
srand(10)
c,demands,capacity = create_inputs(n)
# First we create some columns through enumeration
sub_tours = Array{Array{Int64}}(n)
# Create the subtours which only have one stop first
for i=1:n
    sub_tours[i] = [i]
end
# Now create random subtours
for i=1:5
    sub_tours = [sub_tours; [create_subtour(c,demands,capacity)]]
end
# Delete duplicate tours
sub_tours = unique(sub_tours)

cost_array,customer_matrix,index_dict = create_columns(sub_tours,c,n)

# rmp = Model(solver = CplexSolver())
rmp = Model(solver = GLPKSolverLP(method=:Exact))

num_var = length(cost_array)

@variable(rmp, λ[1:num_var] >= 0)

@objective(rmp, Min,sum( cost_array[r]*λ[r] for r=1:num_var ))

@constraint(rmp, all_customers_serviced[i=1:size(customer_matrix,1)],sum( customer_matrix[i,r]*λ[r] for r=1:num_var ) == 1)

@constraint(rmp, num_subtour_constraint,sum( λ[r] for r=1:num_var ) <= K)
solve(rmp)
println(getvalue(λ))
u = getdual(all_customers_serviced)
u = [0; u; 0]
σ = getdual(num_subtour_constraint)

subprob = Model(solver = CplexSolver())
@variable(subprob, x[i=0:n+1,j=0:n+1],Bin)
@variable(subprob, demands[i] <= y[i=0:n+1] <= capacity)
@objective(subprob, Min,sum( (c[i,j] - u[i+1])*x[i,j] for i=0:n+1,j=0:n+1 ) - σ)
for i=1:n
    # @constraint(subprob, sum( x[i,j] for j=1:n if j!=i ) == 1)
    @constraint(subprob, sum( x[h,i] for h=0:n if h!=i ) - sum( x[i,j] for j=1:n+1 if j!=i ) == 0)
end

@constraint(subprob, sum( x[0,j] for j=1:n ) == 1)
@constraint(subprob, sum( x[i,0] for i=1:n+1 ) == 0)

for i=0:n+1
    for j=0:n+1
        @constraint(subprob, y[j] >= y[i] + demands[j]*x[i,j] - capacity*(1-x[i,j]))
    end
end
solve(subprob)
x = getvalue(x)
println(getobjectivevalue(subprob))
new_subtour = get_subtour(x,n)

new_c,new_a_col = get_subtour_value(new_subtour,c)
xnew = Variable[]
@variable(rmp, 0 <= xnew,
    objective=new_c,
    inconstraints=all_customers_serviced,
    coefficients=new_a_col)
