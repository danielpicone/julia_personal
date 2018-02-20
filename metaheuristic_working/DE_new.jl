# This is a script for the differential evolution metaheuristic

ENV["MPLBACKEND"]="qt4agg"

using PyPlot
include("DE.jl")

# Define the function which you want to minimise

# f(x,y) = x^2 + y^2
f(x,y,z) = (1-x)^2 + 100*(y-x^2)^2 + 100*(z-y^2)^2
# f(x,y,z) = (1-x+y^2)^2 + 4*(y-z^2+x-4)^4
# f(x,y) = -( sin(x)*(sin(x^2/π))^20 + sin(y)*(sin(2*y^2/π))^20 )
# f(x) = 1/14 * (x + 4)*(x + 1)*(x - 1)*(x - 3)
# f(x,y) = (1-x*y)^2 + x^2
# f(x) = -(x^4 + 4*x^3 - 6*x^2 - 4*x + 12)
global dim = 3

# Constraints by the penalty method
function constraints(x...)
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
    numineqconstraints = 1
    ineq = Array{Float64}(numineqconstraints)
    if numineqconstraints!=0
        ineq[1] = -x[1] - x[2] - x[3] + 4
    end
    return 10000000*(sum(eq.^2)+sum(max.(ineq,0)))
end

function de_new(f; popnum=10, maxgen=10000, scheme="rand", maxnomovement = 100, graph="no")
    # First we create the initial population
    # popnum = 10             # Population number at each generation
    # maxgen = 1000000      # Maximum number of generations
    # x = Array{Array{Float64,1}}(popnum,maxgen)
    x = Array{Float64}(popnum,dim,maxgen)
    lb = -10
    ub = 10
    # x[:,1] = [(ub-lb).*rand(dim)+lb for i=1:popnum]
    x[:,:,1] = rand(popnum,dim)
    values = [f(x[i,:,1]...)+constraints.(x[i,:,1]...) for i=1:popnum]

    if graph=="yes"
        size = 1000
        xgrid = linspace(lb,ub,size)
        ygrid = linspace(lb,ub,size)
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
                donorvec = x[indexes[1],:,gen] + F*(x[indexes[2],:,gen] - x[indexes[3],:,gen])
            elseif scheme=="best"
                currentvalues = [f(x[i,:,gen]...)+constraints.(x[i,:,gen]...) for i=1:popnum]
                if gen==1
                else
                    currentvalues = [f(x[i,:,gen]...)+constraints.(x[i,:,gen]...) for i=1:popnum]
                end
                bestval = minimum(currentvalues)
                solution = x[find(currentvalues .== bestval)[1],:,gen]
                donorvec = solution + F*(x[indexes[2],:,gen] - x[indexes[3],:,gen])
            end
            randindex = rand(1:dim)
            randnum = rand()
            newsol = Array{Float64}(dim)
            for j=1:dim
                if (randnum <= Cr || j==randindex)
                    newsol[j] = donorvec[j]
                else
                    newsol[j] = x[i,:,gen][j]
                end
            end
            # Update the solution if it is better
            if f(newsol...)+constraints.(newsol...) <= f(x[i,:,gen]...)+constraints.(x[i,:,gen]...)
                x[i,:,gen+1] = newsol
            else
                x[i,:,gen+1] = x[i,:,gen]
            end
        end
        # Stop the algorithm if the best solution is not moving much
        if gen==1
        else
            if norm(minimum([f(x[i,:,gen]...)+constraints.(x[i,:,gen]...) for i=1:popnum])-minimum([f(x[i,:,gen-1]...)+constraints.(x[i,:,gen-1]...) for i=1:popnum]))<=10.0^-10
                criteria += 1
                # println("Does this run??")
                if criteria>=maxnomovement
                    println("Stopped due to convergence")
                    println("Generations: $gen")
                    values = [f(x[i,:,j]...)+constraints.(x[i,:,j]...) for i=1:popnum,j=1:gen]
                    bestval = minimum(values)
                    println(find(values .== bestval)[1])
                    solutionindex = find(values .== bestval)[1]
                    solutionpop = solutionindex % popnum
                    if solutionpop==0
                        solutionpop = popnum
                    end
                    # solution = x[solutionindex:solutionindex+dim-1]
                    solution = x[solutionpop,:,max(Int64(floor(solutionindex/(popnum))),1)]
                    finalx = x[1:popnum,1:dim,1:gen]
                    println("Local minimum of $bestval found at $solution")
                    return finalx, solutionindex, values, bestval
                end

            else
                criteria = 1
            end
        end
        # For plotting purposes only
        if graph == "yes"
            xpos = []
            for i=1:popnum
                xpos = push!(xpos,x[i,:,gen][1])
            end
            ypos = []
            for i=1:popnum
                ypos = push!(ypos,x[i,:,gen][2])
            end
            plot(xpos,ypos,"k.",markersize=3)
            sleep(0.01)
        end

        gen += 1
    end
    # Use the best vector
    values = [f(x[i,:,j]...)+constraints.(x[i,:,j]...) for i=1:popnum,j=1:maxgen]
    bestval = minimum(values)
    solutionindex = find(values .== bestval)[1]
    solution = x[solutionindex % popnum,:,max(Int64(floor(solutionindex/(popnum*dim))),1)]
    finalx = x[1:popnum,1:dim,1:gen-1]
    println("Stopped due to maximum generations reached")
    println("Generations: $gen")
    println("Local minimum of $bestval found at $solution")
    return finalx, solutionindex, values, bestval
end
tic()
x, solution, values, bestval = de_new(f)
toc()

tic()
x, solution, values, bestval = de(f)
toc()

counter = 0
for k = 1:100
    newdetime1 = time()
    de_new(f);
    newdetime2 = time()
    newdetime = newdetime2 - newdetime1

    detime1 = time()
    de(f);
    detime2 = time()
    detime = detime2 - detime1

    if detime < newdetime
        counter += 1
    end
end

println("The old de version was faster in $counter many runs")
