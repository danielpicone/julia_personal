# Solving a Sudoku puzzle with Integer Programming

using JuMP
using GLPKMathProgInterface
using Cbc
start = time()
sudoku = Model(solver = GLPKSolverMIP())
#sudoku = Model(solver = CbcSolver())
#@defVar(sudoku, x[i=1:9, j=1:9, k=1:9], Bin)
@variable(sudoku, x[i=1:9, j=1:9, k=1:9], Bin)

for i=1:9, j=1:9 # Each row and column
  # Sum across all possible digits
  # This means that only one digit can belong in each position
  #@addConstraint(sudoku, sum(x[i,j,k] for k=1:9) == 1)
  @constraint(sudoku, sum(x[i,j,k] for k=1:9) == 1)
end

for ind = 1:9
  for k = 1:9
    # Sum across columns (j)
    #@addConstraint(sudoku, sum(x[ind,j,k] for j=1:9) == 1)
    @constraint(sudoku, sum(x[ind,j,k] for j=1:9) == 1)
    # Sum across rows (i)
    #@addConstraint(sudoku, sum(x[i,ind,k] for i=1:9) == 1)
    @constraint(sudoku, sum(x[i,ind,k] for i=1:9) == 1)
  end
end


for i = 1:3:7, j = 1:3:7, k = 1:9
    # i is the top left row, j is the top left column
    # We'll sum from i to i+2, e.g. i=4, r=4, 5, 6
    #@addConstraint(sudoku, sum(x[r,c,k] for r=i:i+2 for c=j:j+2) == 1)
    @constraint(sudoku, sum(x[r,c,k] for r=i:i+2 for c=j:j+2) == 1)
end

# The given digits
init_sol = [ 5 3 0 0 7 0 0 0 0;
             6 0 0 1 9 5 0 0 0;
             0 9 8 0 0 0 0 6 0;
             8 0 0 0 6 0 0 0 3;
             4 0 0 8 0 3 0 0 1;
             7 0 0 0 2 0 0 0 6;
             0 6 0 0 0 0 2 8 0;
             0 0 0 4 1 9 0 0 5;
             0 0 0 0 8 0 0 7 9]
for i = 1:9, j = 1:9
    # If the space isn't empty
    if init_sol[i,j] != 0
        # Then the corresponding variable for that digit
        # and location must be 1
        #@addConstraint(sudoku, x[i,j,init_sol[i,j]] == 1)
        @constraint(sudoku, x[i,j,init_sol[i,j]] == 1)
    end
end

solve(sudoku)
#solver = GLPKSolverMIP(sudoku)
#CbcSolver(sudoku)

# Extract the values of x
x_val = getvalue(x)
# Create a matrix to store the solution
sol = zeros(Int,9,9)  # 9x9 matrix of integers
for i in 1:9, j in 1:9, k in 1:9
    # Integer programs are solved as a series of linear programs
    # so the values might not be precisely 0 and 1. We can just
    # round them to the nearest integer to make it easier
    #if iround(x_val[i,j,k]) == 1
    if x_val[i,j,k] >= 0.7
        sol[i,j] = k
    end
end
stop = time()
println(stop - start)
# Display the solution
println(sol)
