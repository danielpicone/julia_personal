# Mean shift algorithm



# Create the xi's
num_points = 100
xi = randn(num_points)

# Choose the kernel function
c = 4
# K(x) = exp(-c*norm(x)^2)
function K(x)
    lambda = 0.5
    if x <= lambda
        return 1
    else
        return 0
    end
end



x = 0
prev_mx = 10
mx = prev_mx

# for i=1:10
while (abs(x-prev_mx) > 10^(-10.0))
    prev_mx = mx
    mx = sum(K.(xi-x).*xi)/sum(K.(xi-x))
    println(mx)
    x = mx
end
true_mean = mean(xi)
println("Final mx: $mx")
println("True mean: $true_mean")
