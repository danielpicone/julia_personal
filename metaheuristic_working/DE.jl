# This is a script for the differential evolution metaheuristic

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# Define the function which you want to minimise

# f(x,y) = x^2 + y^2
# f(x,y,z) = (1-x)^2 + 100*(y-x^2)^2 + 100*(z-y^2)^2
# f(x,y,z) = (1-x+y^2)^2 + 4*(y-z^2+x-4)^4
# f(x,y) = -( sin(x)*(sin(x^2/π))^20 + sin(y)*(sin(2*y^2/π))^20 )
# f(x) = 1/14 * (x + 4)*(x + 1)*(x - 1)*(x - 3)
# f(x,y) = (1-x*y)^2 + x^2
# f(x) = -(x^4 + 4*x^3 - 6*x^2 - 4*x + 12)
# f(x...) = x[1]^2*x[2]*(2+x[3])
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

function de(f; bounds=zeros(dim,2), popnum=10, maxgen=10000, scheme="rand", maxnomovement = 100, graph="no", verbose="no")
    # First we create the initial population
    # popnum = 10             # Population number at each generation
    # maxgen = 1000000      # Maximum number of generations
    x = Array{Array{Float64,1}}(popnum,maxgen)
    x[:,1] = [(bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]
    values = [f(x[i,1]...)+constraints.(x[i,1]...) for i=1:popnum]

    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end

    # Set the parameters
    F = 0.8
    Cr = 0.9                # Probability of crossover
    gen = 1
    global criteria = 1
    while (gen < maxgen && criteria <= maxnomovement)
        for i=1:popnum
            # Choose three vectors at random
            indexes = randperm(popnum)[1:3]
            if scheme=="rand"
                donorvec = x[indexes[1],gen] + F*(x[indexes[2],gen] - x[indexes[3],gen])
            elseif scheme=="best"
                currentvalues = [f(x[i,gen]...)+constraints.(x[i,gen]...) for i=1:popnum]
                if gen==1
                else
                    currentvalues = [f(x[i,gen]...)+constraints.(x[i,gen]...) for i=1:popnum]
                end
                bestval = minimum(currentvalues)
                solution = x[find(currentvalues .== bestval)[1],gen]
                donorvec = solution + F*(x[indexes[2],gen] - x[indexes[3],gen])
            end
            randindex = rand(1:dim)
            randnum = rand()
            newsol = Array{Float64}(dim)
            for j=1:dim
                if (randnum <= Cr || j==randindex)
                    newsol[j] = donorvec[j]
                else
                    newsol[j] = x[i,gen][j]
                end
            end
            # Update the solution if it is better
            if f(newsol...)+constraints.(newsol...) <= f(x[i,gen]...)+constraints.(x[i,gen]...)
                x[i,gen+1] = newsol
            else
                x[i,gen+1] = x[i,gen]
            end
        end
        # Stop the algorithm if the best solution is not moving much
        if gen==1
        else
            if norm(minimum([f(x[i,gen]...)+constraints.(x[i,gen]...) for i=1:popnum])-minimum([f(x[i,gen-1]...)+constraints.(x[i,gen-1]...) for i=1:popnum]))<=10.0^-10
                criteria += 1
                if criteria>=maxnomovement

                    values = [f(x[i,j]...)+constraints.(x[i,j]...) for i=1:popnum,j=1:gen]
                    bestval = minimum(values)
                    solution = x[find(values .== bestval)][1]
                    if graph == "yes"
                        plot(solution[1],solution[2],"ro",markersize=5)
                    end
                    finalx = x[1:popnum,1:gen]
                    if verbose=="yes"
                        println("Stopped due to convergence")
                        println("Generations: $gen")
                        println("Local minimum of $bestval found at $solution")
                    end
                    return finalx, solution, values, bestval
                end

            else
                criteria = 1
            end
        end
        # For plotting purposes only
        if graph == "yes"
            xpos = []
            for i=1:popnum
                xpos = push!(xpos,x[i,gen][1])
            end
            ypos = []
            for i=1:popnum
                ypos = push!(ypos,x[i,gen][2])
            end
            plot(xpos,ypos,"k.",markersize=3)
            sleep(0.01)
        end

        gen += 1
    end
    # Use the best vector
    values = [f(x[i,j]...)+constraints.(x[i,j]...) for i=1:popnum,j=1:maxgen]
    bestval = minimum(values)
    solution = x[find(values .== bestval)][1]
    if graph == "yes"
        plot(solution[1],solution[2],"ro",markersize=5)
    end
    finalx = x[1:popnum,1:gen-1]
    if verbose=="yes"
        println("Stopped due to maximum generations reached")
        println("Generations: $gen")
        println("Local minimum of $bestval found at $solution")
    end
    return finalx, solution, values, bestval;
end

de(f, maxgen = 10);
