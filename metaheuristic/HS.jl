# This is a script for harmony search

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# Define the function which you want to minimise

# f(x) = x^2
# f(x,y) = x^2 + y^2
# f(x,y) = (1-x)^2+100*(y-x^2)^2
# f(x,y,z) = (1-x)^2 + 100*(y-x^2)^2 + 100*(z-y^2)^2
# f(x,y,z) = (1-x+y^2)^2 + 4*(y-z^2+x-4)^4
# f(x,y) = -( sin(x)*(sin(x^2/π))^20 + sin(y)*(sin(2*y^2/π))^20 )
# f(x) = 1/14 * (x + 4)*(x + 1)*(x - 1)*(x - 3)
# f(x,y) = (1-x*y)^2 + x^2
# f(x) = -(x^4 + 4*x^3 - 6*x^2 - 4*x + 12)
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
#     numineqconstraints = 0
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

function hs(f; bounds=zeros(dim,2), popnum = 20, maxiter = 10000, HMacceptrate = 0.95, PArate = 0.7, maxcounter = 200, graph = "no", verbose = "no")
    # Initialise the points
    # HM = Array{Array{Float64,1}}(popnum,maxiter)
    HM = Array{Float64}(popnum,dim)
    pabounds = 200*ones(dim)
    [HM[i,:] = (bounds[:,2]-bounds[:,1]).*rand(dim)+bounds[:,1] for i=1:popnum]
    HMacceptrate = 0.95
    PArate = 0.7
    values = [f(HM[i,:]...)+constraints.(HM[i,:]..., bounds = bounds) for i=1:popnum]
    iter = 1; counter = 1; timer = 0
    newsol = Array{Float64}(dim)

    # Plot the function
    if graph=="yes"
        size = 1000
        xgrid = linspace(bounds[1,1],bounds[1,2],size)
        ygrid = linspace(bounds[2,1],bounds[2,2],size)
        z = [f(xgrid[i],ygrid[j]) for i=1:size,j=1:size]
        contour(xgrid,ygrid,z)
    end



    while (iter < maxiter && counter < maxcounter)
        for j=1:dim
            # Randomly search
            if (rand() >= HMacceptrate)
                newsol[j] = bounds[j,1]+(bounds[j,2]-bounds[j,1])*rand()
            else
                # Harmony Memory Accepting rate
                newsol[j] = HM[Int64(trunc(popnum*rand()+1)),j]
                if (rand() <= PArate)
                    # Pitch Adjusting in a given bounds
                    pa = (bounds[j,2]-bounds[j,1])/pabounds[j]
                    newsol[j] = newsol[j] + pa*(rand()-0.5)
                end
            end
        end
        if graph == "yes"
            if iter % Int64(floor(maxiter/500)) == 0
                xpos = []
                for i=1:popnum
                    xpos = push!(xpos,newsol[1])
                end
                ypos = []
                for i=1:popnum
                    ypos = push!(ypos,newsol[2])
                end
                # plot(newsol[1],newsol[2],"k.",markersize=3)
                
                plot(xpos,ypos,"k.",markersize=3)
                sleep(0.005)
            end
        end
        # Evaluate the new solution
        newval = f(newsol...)+constraints.(newsol...,  bounds = bounds)
        if any(newval .<= values)
            counter = 1
            index = find(values .== maximum(values))[1]
            values[index] = newval
            HM[index,:] = newsol
            # if graph == "yes"
            #     xpos = []
            #     for i=1:popnum
            #         xpos = push!(xpos,HM[index,1])
            #     end
            #     ypos = []
            #     for i=1:popnum
            #         ypos = push!(ypos,HM[index,2])
            #     end
            #     plot(xpos,ypos,"k.",markersize=3)
            #     sleep(0.01)
            # end
        end

        # if graph == "yes"
        #     xpos = []
        #     for i=1:popnum
        #         xpos = push!(xpos,HM[i,1])
        #     end
        #     ypos = []
        #     for i=1:popnum
        #         ypos = push!(ypos,HM[i,2])
        #     end
        #     plot(xpos,ypos,"k.",markersize=3)
        #     sleep(0.01)
        # end


        iter += 1
        counter += 1

        if iter % Int64(floor(maxiter/5)) == 0 && verbose=="yes"
            timer += 1
            completed = timer*20
            println("We are $completed% of the way through")
        end
    end
    bestpos = HM[find(values .== minimum(values))[1],:]
    if graph=="yes"
        plot(bestpos[1],bestpos[2],"ro",markersize=5)
    end
    bestval = f(bestpos...)
    if verbose=="yes"
        println("The best solution was $bestval found at $bestpos")
    end
    return HM, bestpos, bestval
end


hs(f,maxiter = 10);
