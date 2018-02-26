# This is a script for the firefly algorithm

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# Define the function you want to minimise

global dim = 12
# f(x,y) = x^2 + y^2
# f(x,y) = ( 1.5 - x + x * y )^2 + (2.25 - x + x * y^2)^2 + (2.625 - x + x * y^3)^2
# f(x,y...) = x[1] + y[1]^2 + y[2]
# g(x) = x^2 + 1
bounds = Array{Float64}(dim,2)
bounds[:,1] = -4.5*ones(dim)
bounds[:,2] = 4.5*ones(dim)

function constraints(x...; bounds = [-Inf*ones(dim) Inf*ones(dim)])
    # First the equality constraints
    # These are in the form g(x) = 0
    numeqconstraints = 0
    eq = Array{Float64}(numeqconstraints)
    if numeqconstraints!=0
        eq[1] = x[1]^2-1
        # eq[2] = sin(x[1])
    end
    # Now the inequality constraints
    # These are in the form g(x) ≤ 0
    numineqconstraints = 0
    ineq = Array{Float64}(numineqconstraints)
    if numineqconstraints!=0
        ineq[1] = 0.25 - x[2]
        # ineq[1] = 1 - (x[2]^2 * x[3])/(71785*x[1]^4)
        # ineq[2] = (4*x[2]^2 -x[1]*x[2])/(12566*(x[1]^3 * x[2] - x[1]^4)) + 1/(5108*x[1]^2) - 1
        # ineq[3] = 1 - (140.45*x[1])/(x[2]^2 * x[3])
        # ineq[4] = (x[1] + x[2])/1.5 - 1
    end
    lbbound = Array{Float64}(dim)
    ubbound = Array{Float64}(dim)
    for k=1:dim
        lbbound[k] = bounds[k,1] - x[k]
        ubbound[k] = -bounds[k,2] + x[k]
    end

    return 10000000*(sum(eq.^2)+sum(max.(ineq,0)) + sum(max.(lbbound,0) + max.(ubbound,0)))
end

function intensity(a,b,γ; β0 = 1, βmin = 0.3)
    return max(β0*exp(-γ*norm(a-b)^2),βmin)

end

"""
    Takes in a function and attempts to find the global minimum
        * bounds = zeros(dim,2)         lower and upper bound of the variables
        * popnum = 20                   number of fireflies used
        * maxiter = 5000                maximum number of iterations
        * γ = 1                         rate of decay of the intensity of fireflies
        * α0 = 0.5                      initial α
        * αInf = 0.3                    final α which decreases randomness from α0
        * maxcounter = 100              algorithm will stop if the best sol hasn't changed in this many iterations
        * graph="no/yes"                plots the contour plot for a 3d function
        * verbose="no/yes"              prints the output of the solution

        return x, values, bestpos, bestval
"""
function fa(f; bounds = zeros(dim,2), popnum = 20, maxiter = 5000, γ = 1, α0 = 0.5, αInf = 0.3, maxcounter = 100, graph="no", verbose="no")
    # Initialise fireflies
    x = Array{Array{Float64,1}}(popnum,maxiter)
    # bounds = Array{Float64}(dim,2)
    # bounds[:,1] = -2*π*ones(dim)
    # bounds[:,2] = 2*π*ones(dim)
    # bounds[:,1] = [1.0 2.0]
    # bounds[:,2] = [2.0 5.0]
    # bounds[:,1] = [0.05 0.25 0.2]
    # bounds[:,2] = [2.0 1.3 15.0]


    x[:,1] = [(bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]

    # x[:,1] = [([5,5]-[-1,-1]).*rand(dim)+[-1,-1] for i=1:popnum]

    values = Array{Float64}(popnum,maxiter)
    println("Does this work in fa??")
    println(bounds)
    println(constraints.(x[1,1]..., bounds = bounds))
    println("hello")
    # values[:,1] = [f(x[i,1]...)+constraints.(x[i,1]..., bounds = bounds) for i=1:popnum]
    values[:,1] = [f(x[i,1])+constraints.(x[i,1], bounds = bounds) for i=1:popnum]
    println("Does this work in fa??2")
    α = α0
    iter = 1; counter = 1

    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end


    while (iter < maxiter && counter < maxcounter)
        for i=1:popnum
            for j=1:popnum

                if values[i,iter] > values[j,iter]
                    counter = 1
                    x[i,iter+1] = x[i,iter] + intensity(x[i,iter],x[j,iter],γ)*(x[j,iter]-x[i,iter]) + α*(rand(dim)-0.5)
                else
                    # println(isassigned(x,i+popnum*(iter)))

                    if isassigned(x,i+popnum*(iter))
                    else
                        x[i,iter+1] = x[i,iter]
                        counter += 1
                    end
                end
            end
        end
        if graph == "yes" && (iter < 50 || iter % 100 == 0) && iter < 1000
            # xpos = []
            # for i=1:popnum
            #     xpos = push!(xpos,x[i,iter][1])
            # end
            # ypos = []
            # for i=1:popnum
            #     ypos = push!(ypos,x[i,iter][2])
            # end
            # plot(xpos,ypos,"k.",markersize=3)
            for i=1:popnum
                plot(x[i,iter][1],x[i,iter][2],"k.",markersize=3)
            end
            sleep(0.001)
        end

        values[:,iter+1] = [f(x[i,iter+1]...)+constraints.(x[i,iter+1]..., bounds = bounds) for i=1:popnum]
        α = αInf + (α0 - αInf)*exp(-iter)
        if verbose=="yes"
            if iter % Int64(0.1*maxiter) == 0
                println("We are up to iteration $iter out of $maxiter")
            end
        end
        iter+=1


    end
    bestval = minimum(values[:,iter])
    bestpos = x[find(bestval .== values)[1]]

    if verbose=="yes"
        println("The best objective value was $bestval located at $bestpos")
    end
    return x, values, bestpos, bestval;
end

# fa(f, maxiter = 100);
