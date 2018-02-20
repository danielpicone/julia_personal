# Column Generation of eta
# This is a simple example to illustrate how the column generation algorithm works
# The example is taken from page 190 of Wolsey

using JuMP
using CPLEX
using GLPKMathProgInterface

num_of_cities = 5

D = [0 7 2 1 5;
     7 0 3 6 8;
     2 3 0 4 2;
     1 6 4 0 9;
     5 8 2 9 0]

#LPM = Model(solver = CplexSolver())
LPM = Model(solver = GLPKSolverLP(method=:Exact))
#LPM = Model(solver = GLPKSolverLP())

c = [28,25,21,19,22,18,28]

@variable(LPM,lambda[1:7]>=0)
@objective(LPM, Min, sum( c[i]*lambda[i] for i=1:7))

A = [2 2 2 2 2 2 2;
     2 2 2 1 1 2 3;
     2 3 2 3 2 3 1;
     2 2 3 3 3 1 1;
     2 1 1 1 2 2 3]
b = 2*ones(1,5)
@constraintref consName[1:5]
for i=1:5
    consName[i] = @constraint(LPM, (A*lambda)[i] == b[i])
    #@constraint(LPM, (A*lambda)[i] == b[i])
end

solve(LPM)

uvec = getdual(consName)

LPM_dual = Model(solver = CplexSolver())

@variable(LPM_dual, u[1:5])
@objective(LPM_dual, Max, -sum(2*u[i] for i=1:5))

@constraintref consNameDual[1:7]
for i=1:7
    consNameDual[i] = @constraint(LPM_dual, c[i] + (A'*u)[i] >= 0)
end

solve(LPM_dual)

#A = rand(5,5)
for i=1:5
    #@constraint(LPM_dual, (A*u)[i] >= 0)
end

u = getvalue(u)
redCost = zeros(5,5)
for i=1:5
    for j=i:5
        redCost[i,j] = D[i,j] - uvec[i] - uvec[j]
        redCost[j,i] = redCost[i,j]
    end
end

# Solve for eta_1
num_of_cities = 5
eta = Model(solver = CplexSolver())
@variable(eta, x[i=1:num_of_cities,j=1:num_of_cities],Bin)
# The objective is to minimise the distance of the path
#@objective(eta, Min, sum( D[i,j]*x[i,j] for i=1:num_of_cities,j=i+1:num_of_cities ))
@objective(eta, Min, sum( redCost[i,j]*x[i,j] for i=1:num_of_cities,j=i+1:num_of_cities ))
# Add the constraints
# Path entering one node must be 1

for i=1:num_of_cities
    @constraint(eta, sum(x[i,j] for j=1:num_of_cities) == 2)
    @constraint(eta, x[i,i] == 0)
    for j=i+1:num_of_cities
        @constraint(eta, x[i,j] == x[j,i])
    end
end

solve(eta)

x = getvalue(x)
total = 0
for i=1:num_of_cities
    for j=i+1:num_of_cities
        total += D[i,j]*x[i,j]
    end
end
