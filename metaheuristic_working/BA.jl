# This is a script for the bat algorithm

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# Define the function you want to minimise

# f(x,y) = x^2 + y^2
# f(x) = x^2
# f(x,y,z) = (1-x)^2 + 100*(y-x^2)^2 + 100*(z-y^2)^2
# f(x...) = x[1]^2*x[2]*(2+x[3])
# f(x,y,z) = (1-x+y^2)^2 + 4*(y-z^2+x-4)^4
# f(x,y) = -( sin(x)*(sin(x^2/π))^20 + sin(y)*(sin(2*y^2/π))^20 )
# f(x) = 1/14 * (x + 4)*(x + 1)*(x - 1)*(x - 3)
# f(x,y) = (1-x*y)^2 + x^2
# f(x,y) = -( exp(-(x-4)^2-(y-4)^2) + exp(-(x+4)^2 - (y-4)^2) + 2*(exp(-x^2 - y^2) + exp(-x^2 - (y+4)^2)) )
# f(x,y) = x^2 + y^2 + 25*(sin(x)^2 + sin(y)^2)
# global dim = 2
# bounds = Array{Float64}(dim,2)
# bounds[:,1] = -2*π*ones(dim)
# bounds[:,2] = 2*π*ones(dim)

# Constraints by the penalty method
# function constraints(x...; bounds = [-Inf*ones(dim) Inf*ones(dim)])
#     # First the equality constraints
#     # These are in the form g(x) = 0
#     numeqconstraints = 0
#     eq = Array{Float64}(numeqconstraints)
#     if numeqconstraints!=0
#         eq[1] = x[1]^2-1
#         # eq[2] = sin(x[1])
#     end
#     # Now the inequality constraints
#     # These are in the form g(x) ≤ 0
#     numineqconstraints = 0
#     ineq = Array{Float64}(numineqconstraints)
#     if numineqconstraints!=0
#         ineq[1] = 0.25 - x[2]
#         # ineq[1] = 1 - (x[2]^2 * x[3])/(71785*x[1]^4)
#         # ineq[2] = (4*x[2]^2 -x[1]*x[2])/(12566*(x[1]^3 * x[2] - x[1]^4)) + 1/(5108*x[1]^2) - 1
#         # ineq[3] = 1 - (140.45*x[1])/(x[2]^2 * x[3])
#         # ineq[4] = (x[1] + x[2])/1.5 - 1
#     end
#     lbbound = Array{Float64}(dim)
#     ubbound = Array{Float64}(dim)
#     for k=1:dim
#         lbbound[k] = bounds[k,1] - x[k]
#         ubbound[k] = -bounds[k,2] + x[k]
#     end
#
#     return 10000000*(sum(eq.^2)+sum(max.(ineq,0)) + sum(max.(lbbound,0) + max.(ubbound,0)))
# end
"""
    Takes in a function and attempts to find the global minimum
        kwargs:
        * bounds = zeros(dim,2)         lower and upper bounds of the variables
        * popnum = 10                   number of bats used by the algorithm
        * maxiter = 4000                maximum number of iterations performed
        * A = 0.25                      increase to increase change of better sol being chosen
        * pulse = 0.5                   increase to decrease chance of random sol being chosen
        * graph="no/yes"                graphs the contour plot for a 3d function
        * verbose="no/yes"              prints output of the solution

        return x, values, bestsol, bestval
"""
function ba(f; bounds=zeros(dim,2), popnum = 10, maxiter = 4000, A = 0.25, pulse = 0.5, graph="no", verbose="no")
    # Initialise bats
    x = Array{Array{Float64,1}}(popnum,maxiter)
    x[:,1] = [(bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]
    values = Array{Float64}(popnum,maxiter)

    values[:,1] = [f(x[i,1]...)+constraints.(x[i,1]..., bounds = bounds) for i=1:popnum]
    bestsol = x[find(minimum(values[:,1] .== values[:,1]))[1]]
    bestval = values[find(minimum(values[:,1] .== values[:,1]))[1]]
    # Maximum and minimum frequencies
    qmin = 0
    qmax = 2

    q = zeros(popnum)
    velocity = Array{Array{Float64,1}}(popnum,maxiter)
    [velocity[i,1] = zeros(dim) for i=1:popnum]

    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end


    newsol = 0
    iter = 1
    while iter < maxiter
        for i=1:popnum
            q[i] = qmin + (qmin - qmax)*rand()
            # println("This is run i= ",i,", iter = ",iter)
            velocity[i,iter+1] = velocity[i,iter] + (x[i,iter] - bestsol)*q[i]
            newsol = x[i,iter]+velocity[i,iter+1]

            # newsol[i,iter+1] = x[i,iter]+velocity[i,iter+1]
            # Pulse rate
            if rand() > pulse
                newsol=bestsol+0.01*randn(dim)
            end

            # Evaluate the new solutions
            newval = f(newsol...)+constraints.(newsol..., bounds = bounds)

            if (newval <= values[i,iter] && rand() < A)
                x[i,iter+1] = newsol
                if graph == "yes" # && (iter < 50 || iter % 50 == 0) && iter < 1000
                    plot(newsol[1],newsol[2],"k.",markersize=3)
                    sleep(0.001)
                end
            else
                x[i,iter+1] = x[i,iter]
            end

            # Update the current best
            if newval <= bestval
                bestsol = newsol
                bestval = newval
            end
        end

        if verbose=="yes"
            if iter % Int64(0.1*maxiter) == 0
                println("We are up to iteration $iter out of $maxiter")
            end
        end
        iter += 1
    end

    for j=1:maxiter
        values[:,j] = [f(x[i,j]...)+constraints.(x[i,j]..., bounds = bounds) for i=1:popnum]
    end
    if verbose=="yes"
        println("The best solution found is $bestval located at $bestsol")
    end
    return x, values, bestsol, bestval
end

ba(f,maxiter = 10);
