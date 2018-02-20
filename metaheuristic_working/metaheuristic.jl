# This is a script which keeps all the metaheuristic algorithms

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

include("DE.jl");    #Differential Evolution
include("PSO.jl");   #Particle Swarm Optimisation
include("HS.jl");    #Harmony Search
include("FA.jl");    # Firefly Algorithm
include("BA.jl");    # Bat Algorithm
include("CS.jl");    # Cuckoo Search

# de(f, verbose="yes")

# struct Polynomial{R}
#     coeffs::Vector{R}
# end
#
# function (p::Polynomial)(x)
#     v = p.coeffs[end]
#     for i = (length(p.coeffs)-1):-1:1
#         v = v*x + p.coeffs[i]
#     end
#     return v
# end
#
# struct Point{R}
#     points::Vector{R}
# end
#
# function (r::Point)(w)
#     # println("w is: ", w)
#     return sum(r.points.*w)^2
# end
#
# r = Point([2,4])
# println(r([2,3]))
# println(bounds)

# p = Polynomial([0,-2,1])
# q = Polynomial([0,-2,1])
# fa(f,bounds = bounds, verbose="yes")

# println("DE")
# de(f, bounds = bounds, verbose="yes", maxnomovement = 1000, maxgen = 10000)
# println("PSO")
# pso(f, bounds = bounds, verbose="yes")
# println("HS")
# hs(f, bounds = bounds, verbose="yes")
# println("FA")
# fa(f, bounds = bounds, verbose="yes")
# println("BA")
# ba(f, bounds = bounds, verbose="yes", maxiter = 10000)
# println("CS")
# cs(f, bounds = bounds, verbose="yes")
