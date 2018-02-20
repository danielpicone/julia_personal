# Kernel Density Estimation
# This is a script which estimates a function f(x) using kernel density estimation

using PyPlot
println("Done loading PyPlot!")
f(x) = sin(4*x)*x-x^2+exp(-x^3 + cos(x))

num_points = 100
data = randn(num_points,2)
a = 2
b = -2
xi = abs(a-b)*(rand(num_points)-0.5)+(a+b)/2
# yi = sin.(4*xi)+randn(num_points)
yi = f.(xi) + randn(num_points)/2
# yi = data[:,2]

λ = abs(a-b)/10

function D(u)
    if abs(u) <= 1
        return 0.75*(1-u^2)
        # return 1-abs(u)
        # return 15/16 * (1-u^2)^2
        # return exp(-1/λ * norm(u)^2)
    else
        return 0
    end
end

K(x0,x) = D(abs(x-x0)/λ)
# function K(x0,x)
#     if abs(x0-x) <= 1
#         return 1-abs(x0-x)
#     else
#         return 0
#     end
# end

# f(x) = sum(K.(x,xi).*yi)/sum(K.(x,xi))
f̂(x) = sum(K.(x,xi).*yi)/sum(K.(x,xi))

# xaxis = linspace(minimum(xi),maximum(xi),1000)
xaxis = linspace(a,b,1000)
plot(xaxis,f̂.(xaxis),xaxis,f.(xaxis),xi,yi,"ro")
