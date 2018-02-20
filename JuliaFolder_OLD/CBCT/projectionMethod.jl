# Project onto a vector containing only {0,1} and satisfying Ax=b

function projectBinary(x)
    round(x)
    for i=1:size(x)[1]
        if x[i] <= 0.5
            x[i] = 0
        else
            x[i] = 1
        end
    end
    return x
end

function projectHyper(x,A,b)
    n = size(A)[1]
    # For stochastic method:
    # u = randperm(n)[1:Int64(n/10)]
    # for i in u
    for i=1:n
        # println(i)
        # println(((A[i,:]'*x - b[i])/(sum(A[i,:].^2))))
        x = x - ((A[i,:]'*x - b[i])/( sum(A[i,:].^2)) ) .* A[i,:]
    end
    return x
end

dimension = 11500
numConstraints = 3418
x = zeros(dimension)
A = randn(numConstraints,dimension)
b = rand(numConstraints)
# A = [1 2]
# b = [4]
# projectBinary(x)
# A = [1 0 1;
#      1 -1 -1]
# b = [1,0]
# for i=1:10000
#     # println(x)
#     x = projectHyper(x,A,b)
#     err = norm(A*x-b)
#     println("The error is: $err")
#     # x = projectBinary(x)
# end
err = 1
numIter = 1
tic()
while err > 10.0^(-3)
    # println(x)
    x = projectHyper(x,A,b)
    err = norm(A*x-b)
    numIter += 1
    if numIter % 1000 == 0
        println("Continuing... $numIter")
        println("The error is: $err")
    end
    # x = projectBinary(x)
end
timeTaken = toc()
println("The number of iterations is:   $numIter")
println("The error is:                  $err")
println("The total time taken is:       $timeTaken")
