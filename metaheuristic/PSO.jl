# This is a script for the particle swarm method of optimisation

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# Define the function you want to minimise

# f(x,y) = x^2 + y^2
# f(x) = x^2
# f(x,y,z) = (1-x)^2 + 100*(y-x^2)^2 + 100*(z-y^2)^2
# f(x,y,z) = (1-x+y^2)^2 + 4*(y-z^2+x-4)^4
# f(x,y) = -( sin(x)*(sin(x^2/π))^20 + sin(y)*(sin(2*y^2/π))^20 )
# f(x) = 1/14 * (x + 4)*(x + 1)*(x - 1)*(x - 3)
# f(x,y) = (1-x*y)^2 + x^2
# f(x) = -(x^4 + 4*x^3 - 6*x^2 - 4*x + 12)
# f(x...) = x[1]^2*x[2]*(2+x[3])
# global dim = 3

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
#     numineqconstraints = 4
#     ineq = Array{Float64}(numineqconstraints)
#     if numineqconstraints!=0
#         ineq[1] = 1 - (x[2]^2 * x[3])/(71785*x[1]^4)
#         ineq[2] = (4*x[2]^2 -x[1]*x[2])/(12566*(x[1]^3 * x[2] - x[1]^4)) + 1/(5108*x[1]^2) - 1
#         ineq[3] = 1 - (140.45*x[1])/(x[2]^2 * x[3])
#         ineq[4] = (x[1] + x[2])/1.5 - 1
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

function pso(f; bounds=zeros(dim,2), popnum = 10, maxiter = 1000, α=0.7, β=0.5, graph="no", maxcounter = 200, verbose="no")
    x = Array{Array{Float64,1}}(popnum,maxiter)
    x[:,1] = [(bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]

    values = Array{Float64}(popnum,maxiter)
    values[:,1] = [f(x[i,1]...) + constraints.(x[i,1]..., bounds = bounds) for i=1:popnum]
    globalbestvalue = minimum(values[:,1])
    # globalbest is a vector with the particle and iteration the best value was achieved
    globalbest = [find(globalbestvalue .== values[:,1])[1], 1]

    # currentbest is a vector which holds the iteration of the best value for each particle
    currentbest = Array{Int64}(popnum)
    currentbest .= 1
    currentbestvalues = [f(x[i,1]...) + constraints.(x[i,1]..., bounds = bounds) for i=1:popnum]
    # velocity = Array{Array{Float64,1}}(popnum,maxiter)
    # velocity[:,1] = [zeros(dim) for i=1:popnum]
    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end
    iter = 1; counter = 1;

    while (iter < maxiter && counter < maxcounter )
        for i=1:popnum
            # Update the velocity
            # velocity[i,iter+1] = 0.7*velocity[i,iter] + α*rand(dim).*(x[globalbest...]-x[i,iter])+β*rand(dim).*(x[i,currentbest[i]]-x[i,iter])
            # x[i,iter+1] = x[i,iter] + velocity[i,iter+1]
            x[i,iter+1] = (1 - β)*x[i,iter] + β*x[globalbest...] + α*randn(dim)
            # Update current best for particle i
            if f(x[i,iter+1]...) + constraints.(x[i,iter+1]..., bounds = bounds) < f(x[i,currentbest[i]]...)+constraints.(x[i,currentbest[i]]..., bounds = bounds)
                currentbest[i] = iter+1
                currentbestvalues[i] = f(x[i,currentbest[i]]...) + constraints.(x[i,currentbest[i]]..., bounds = bounds)
            end

        end
        values[:,iter+1] = [f(x[i,iter+1]...) + constraints.(x[i,iter+1]..., bounds = bounds) for i=1:popnum]

        if any(currentbestvalues .< globalbestvalue)
            counter = 1
            globalbestvalue = minimum(currentbestvalues)
            globalbest[1] = find(globalbestvalue .== currentbestvalues)[1]
            globalbest[2] = iter+1
        else
            counter += 1
            # println(counter)
        end

        # Find best current particle
        currentbestparticle = find(minimum(currentbestvalues) .== currentbestvalues)[1]
        values[:,iter+1] = [f(x[i,iter+1]...) + constraints.(x[i,iter+1]..., bounds = bounds) for i=1:popnum]

        if graph == "yes"
            xpos = []
            for i=1:popnum
                xpos = push!(xpos,x[i,iter][1])
            end
            ypos = []
            for i=1:popnum
                ypos = push!(ypos,x[i,iter][2])
            end
            plot(xpos,ypos,"k.",markersize=3)
            sleep(0.01)
        end
        iter += 1
        α *= 0.9
    end
    if graph == "yes"
        plot(x[globalbest...]...,"ro",markersize=5)
    end
    bestval = f(x[globalbest...]...)+constraints.(x[globalbest...]..., bounds = bounds)
    bestpos = x[globalbest...]
    if verbose=="yes"
        println("The best value is $bestval")
        println("Found at $bestpos")
    end
    return x[:,1:iter], values[:,1:iter], globalbest, globalbestvalue;
end

pso(f, maxiter = 10);
