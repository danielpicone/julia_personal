# Col gen test

# using JuMP
# using CPLEX
#
# model = Model(solver = CplexSolver())
#
# @variable(model, x >= 0)
#
# @objective(model, Min,3*x)
#
# @constraint(model, con1, 0.5<=x)
# solve(model)
#
# println(getvalue(x))
#
# @variable(model,0<=y,-2,ConstraintRef{LinearConstraint}[con1[1]],[1])
# solve(model)


using JuMP
using CPLEX
m = Model(solver = CplexSolver())
@variable(m, 0 <= x <= 1)
@variable(m, 0 <= y <= 1)
@objective(m, Max, 5x + 1y)
constrrefs = Vector{ConstraintRef}(0)
@constraintref(con[1:2])
# for i=1:2
#     con[i] = @constraint(m, i*x + y <= i+5)
# end

for i=1:2
    con[i] = @constraint(m, i*x + y <= i+5)
    push!(constrrefs,con[i])
end

# print(m)
# solve(m)  # x = 1, y = 1
# println("x= ", getvalue(x), " y= ", getvalue(y))

# @variable(m, 0 <= z <= 1, 10.0, ConstraintRef{LinearConstraint}[con[i] for i in [1:2]], [1.0;-2.00])
# @variable(m, 0 <= z <= 1, 10.0, [con], [1.0;-2.00])
@variable(m, 0 <= z <= 1, objective = 10.0, inconstraints = constrrefs, coefficients = [1.0;-2.00])
solve(m)


# using JuMP
# using CPLEX
# m = Model(solver = CplexSolver())
# @constraint(m, con1, 0 <= 1)
# @constraint(m, con2, 0 <= 1)
# constrrefs = Vector{ConstraintRef}(0)
# push!(constrrefs,con1)
# push!(constrrefs,con2)
# values = [1.0, 1.0]
# @variable(m, lowerbound = 0.0, upperbound = 1.0, objective = 10.0, inconstraints = constrrefs, coefficients = values)
# solve(m)
