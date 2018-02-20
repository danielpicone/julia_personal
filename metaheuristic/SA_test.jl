# This is a script to test the simulated annealing algorithm

ENV["MPLBACKEND"]="qt4agg"

using PyPlot

# f(x,y) = (1-x)^2 + 100(y-x^2)^2
# f(x,y) = 8*x^2+4*x*y+3*y^2
# f(x) = -cos(x)*exp(-(x-π)^2)
f(x) = x^2
penalty = 0*10.0^10
# Input the equality constraints as an array of functions
# Should be of the form q(x) = 0 or q(x) <= 0
function constraints(x)
    eq = Array{Float64}(1)
    eq[1] = 0
    ineq = Array{Float64}(2)
    ineq[1] = x-10
    ineq[2] = -x-10
    return eq,ineq
end

function getH(x)
    eq,ineq = constraints(x...)
    eqH = (eq.!=0)
    ineqH = (ineq.>=0)
    return eq,eqH,ineq,ineqH
end

function fun(x)
    eq,eqH,ineq,ineqH = getH(x...)
    return f(x...) + penalty*sum(eqH.*(eq.^2)) + penalty*sum(ineqH.*max.(0,ineq))
end

## Here is the code for Simulated Annealing (SA)

# Global scale of randoms
global scale = 1
global dim = 1

function sa(f; α=0.95, initial=scale*randn(dim))
    # Cooling schedule default is 0.95
    # Initial temperature
    initialtemperature = 10.0^5
    # Maximum iterations
    N = 1000
    # Final stopping temperature
    finaltemperature = 10.0^(-10)
    # Maximum number of rejections
    maxrej = 500
    # Maximum number of runs
    maxrun = 150
    # Maximum number of accepts
    maxaccept = 50
    # Inital search period
    initialsearch = 500
    # Initalise the counters
    i = 0; j = 0; accept = 0; totaleval = 0
    temperature = initialtemperature
    # Initial guess
    # x = initial
    x = [0.0]
    xbest = x
    fbest = fun(xbest...)
    Einital = fun(x...)
    Eold = Einital; Enew = Eold;
    # while ( temperature > finaltemperature && n < N)
    # while i < N
    while (temperature > finaltemperature && totaleval < 1000000)
    # while temperature > finaltemperature
        i += 1
        if mod(totaleval,10000)==0
            # println(temperature)
            # println(totaleval)
        end
        # Check if max number of run/accept are met
        if (i >= maxrun || accept >= maxaccept)
            # reset the counters
            i = 1; accept = 1
            # Cooling according to schedule
            temperature = cooling(temperature,α)
            @printf(".")
        end
        # Function evaluations at new locations
        if totaleval < initialsearch
            initflag = 1
            ns = newsolution(x,initflag)
        else
            initflag = 0
            ns = newsolution(xbest,initflag)
        end
        totaleval += 1
        Enew = fun(ns...)
        # Decide to accept the new solution
        ΔE = Enew-Eold
        # Accept if improved
        if ΔE < 0
            xbest = ns; Eold = Enew
            accept += 1; j = 0
        end
        # Accept with a probability if not improved
        if ( ΔE >= 0 && exp(-ΔE/temperature) > rand() )
            xbest = ns; Eold = Enew;
            accept += 1
        else
            j+=1
        end
        # Update the estimated optimal solution
        fbest = Eold
    end
    # println(maxj)
    @printf("\n")
    println(temperature)
    println(totaleval)
    return xbest,fbest
end

# Cooling schedule
function cooling(t,α)
    return t*α
end

# Create a new solution
function newsolution(x,initflag)
    # Either search around at the start
    if initflag == 1
        s = scale*(rand(dim,1)-0.5)
    else
        stepsize = 0.01
        s = x + stepsize*scale*rand(dim,1)
    end
    return s
end



println("Starting sa now")
output = sa(f)
