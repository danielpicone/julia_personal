# This is a julia script to solve the binary brainsnack problem

# Mostly working, does not have the constraint which ensures that all cols and rows are unique
n = 5

using JuMP
using CPLEX

bbproblem = Model(solver = CplexSolver())

@variable(bbproblem, x[1:2*n,1:2*n,0:1],Bin)

for i=1:2*n
    for j=0:1
        @constraint(bbproblem,sum(x[i,:,j]) == n)
    end
end

for i=1:2*n
    for j=0:1
        @constraint(bbproblem,sum(x[:,i,j]) == n)
    end
end

for i=1:2*n
    for j=0:1
        for k=1:(2*n-2)
            @constraint(bbproblem,sum(x[i,k:k+2,j]) <= 2)
        end
    end
end

for i=1:2*n
    for j=0:1
        for k=1:(2*n-2)
            @constraint(bbproblem,sum(x[k:k+2,i,j]) <= 2)
        end
    end
end

for i=1:2*n
    for j=1:2*n
        @constraint(bbproblem,sum(x[i,j,:]) == 1)
    end
end

initial_sol = [2 1 2 1 2 0 2 2 1 2;
               0 2 0 2 2 2 2 2 2 2;
               2 2 2 2 2 2 1 1 2 2;
               2 2 1 2 2 0 2 2 2 2;
               0 2 2 2 2 2 2 2 2 2;
               2 2 2 2 2 2 2 0 0 2;
               1 2 2 2 2 2 2 1 2 1;
               2 2 0 2 2 2 0 2 2 2;
               2 2 2 2 2 2 0 2 1 2;
               2 2 2 2 0 2 2 0 2 2;
]
# initial_sol = [1 2 2 0 2 2;
#                2 2 0 0 2 1;
#                2 0 0 2 2 1;
#                2 2 2 2 2 2;
#                0 0 2 1 2 2;
#                2 1 2 2 0 0;
# ]

for i=1:2*n
    for j=1:2*n
        if initial_sol[i,j] != 2
            @constraint(bbproblem,x[i,j,initial_sol[i,j]] == 1)
        end
    end
end

@objective(bbproblem,Min,1)
solve(bbproblem)

x_val = getvalue(x)
# Create a matrix to store the solution
sol = zeros(Int,2*n,2*n)  # 9x9 matrix of integers
for i in 1:2*n, j in 1:2*n, k in 0:1
    # Integer programs are solved as a series of linear programs
    # so the values might not be precisely 0 and 1. We can just
    # round them to the nearest integer to make it easier
    #if iround(x_val[i,j,k]) == 1
    if x_val[i,j,k] >= 0.9
        sol[i,j] = k
    end
end
