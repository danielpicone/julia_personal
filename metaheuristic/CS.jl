# This is a script for the cuckoo search algorithm

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
#         # ineq[1] = 0.25 - x[2]
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

function levy(β)
    σ_u = ( (gamma(1+β) * sin(π*β/2)) / (gamma((1+β)/2))*β*2^((β-1)/2)   )^(1/β)
    σ_v = 1

    u = randn()*σ_u
    v = randn()*σ_v

    s = u/(abs(v)^(1/β))

    return s
end

function boundscheck(x,bounds)
    x = min.(max.(x,bounds[:,1]),bounds[:,2])
    return x
end


function cs(f; bounds=zeros(dim,2), popnum = 10, maxiter = 1000, maxcounter = 250, α = 1, pa = 0.25, graph="no", β = 0.5, verbose = "no")
    # Initialise nests
    # x = Array{Array{Float64,1}}(popnum,maxiter)
    x = Array{Array{Float64,1}}(popnum)
    # bounds = Array{Float64}(dim,2)
    # bounds[:,1] = 0*ones(dim)
    # bounds[:,2] = 5*ones(dim)

    [x[i] = (bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]
    # values = Array{Float64}(popnum,maxiter)
    # values[:,1] = [f(x[i,1]...)+constraints.(x[i,1]..., bounds = bounds) for i=1:popnum]

    values = [f(x[i]...)+constraints.(x[i]..., bounds = bounds) for i=1:popnum]
    bestsol = x[find(minimum(values) .== values)[1]]
    bestval = values[find(minimum(values) .== values)[1]]

    abandon = (Int64(round(popnum*(1-pa))))

    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end

    iter = 1; counter = 0;
    while (iter < maxiter && counter < maxcounter)
        # Choose a cuckoo randomly
        index = Int64(floor( popnum*rand() ) + 1)
        direction = (rand(dim).-0.5)
        direction /= norm(direction)
        newsol = x[index] + α*direction*levy(β)
        newsol = boundscheck(newsol,bounds)
        newval = f(newsol...)+constraints.(newsol..., bounds = bounds)
        index = Int64(floor( popnum*rand() ) + 1)
        if graph == "yes" && (iter < 100 || iter % 5 == 0) # && iter < 1000
            plot(newsol[1],newsol[2],"k.",markersize=3)
            sleep(0.001)
        end
        if newval < values[index]
            counter = 1
            x[index] = newsol
            values[index] = newval
            # if graph == "yes" # && (iter < 50 || iter % 50 == 0) && iter < 1000
            #     plot(newsol[1],newsol[2],"k.",markersize=3)
            #     sleep(0.001)
            # end
        end
        sortindex = sortperm(values)
        values = values[sortindex]
        x = x[sortindex]
        for i=popnum:abandon
            direction = (rand(dim).-0.5)
            direction /= norm(direction)
            x[i] = α*direction*levy(β)
            x[i] = boundscheck(x[i],bounds)
            values[i] = f(x[i]...)+constraints.(x[i]..., bounds = bounds)
        end
        if minimum(values) < bestval
            bestval = minimum(values)
            # println("The best value was changed")
        else
            # println(counter)
            counter += 1;
        end
        iter += 1

    end
    bestval = values[find(minimum(values) .== values)[1]]
    bestsol = x[find(minimum(values) .== values)[1]]
    if verbose=="yes"
        if counter >= maxcounter
            println("The algorithm stopped as convergence was reached")
        end
        println("The solution of $bestval was found at $bestsol in $iter iterations")
    end
    return x, values, bestsol, bestval;
end

cs(f, maxiter = 10);
