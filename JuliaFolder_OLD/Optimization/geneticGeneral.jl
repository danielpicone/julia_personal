# This is a script which implements a genetic algorithm to solve an unconstrained optimisation problem
# This is for the single dimension case
dim = 2
Q = [1 2; 2 5]
# Q = rand(dim,dim)
f(x) = x' * Q * x

# Number of points which can breed
numIndividuals = 5
numBreed = 3

function breed!(dominantPoint,otherPoint)
    # Set a 1% chance for mutation of dominantPoint
    # 25% chance for average of parents
    # 25% chance for half of each parent
    # 25% chance for dom - other
    # 25% chance for dom + other
    if rand() < 0.05
        newPoint = dominantPoint
        newPoint[1:end-1] = dominantPoint[1:end-1]+10*rand()
        newPoint[end] = f(newPoint[1:end-1])[1]
        println("We mutated!!")
        return newPoint
    else
        chance = rand()
        if chance <= 0.25
            newPoint = (dominantPoint+otherPoint)./2
            newPoint[end] = f(newPoint[1:end-1])[1]
            if newPoint[end] < otherPoint[end]
                return newPoint
            else
                return otherPoint
            end
        elseif chance <= 0.50
            swapDim = randperm(dim)[Int32(floor(dim/2))]
            newPoint = otherPoint
            newPoint[swapDim] = dominantPoint[swapDim]
            newPoint[end] = f(newPoint[1:end-1])[1]
            if newPoint[end] < otherPoint[end]
                return newPoint
            else
                return otherPoint
            end
        elseif chance <= 0.75
            newPoint = dominantPoint-otherPoint
            newPoint[end] = f(newPoint[1:end-1])[1]
            if newPoint[end] < otherPoint[end]
                return newPoint
            else
                return otherPoint
            end
        else
            newPoint = dominantPoint+otherPoint
            newPoint[end] = f(newPoint[1:end-1])[1]
            if newPoint[end] < otherPoint[end]
                return newPoint
            else
                return otherPoint
            end
        end
    end
    return newPoint
end




# Choose the point with lowest obj value and let it breed with the numBreed
points = Array(Float64,numIndividuals,dim+1)

for j=1:numIndividuals
    points[j,1:dim] = 1*rand(1,dim)
    points[j,end] = f(points[j,1:end-1])[1]
end

iterations = 100
for i=1:iterations
    # Choose numBreed points at random
    # Indicies of breeding points
    ind = randperm(numIndividuals)[1:numBreed]
    # Find the index of the dominant point
    dominantInd = find(x->x==minimum(points[ind,2]),points[:,2])[1]

    for j in ind
        # Now pairwise breed between the dominant and the other points
        if j==dominantInd
        else
            points[j,:] = breed!(points[dominantInd,:],points[j,:])
        end

    end
    # Display the best current point
    bestInd = find(x->x==minimum(points[:,2]),points[:,2])[1]
    println(points[bestInd,:])

end
dominantInd = find(x->x==minimum(points[:,2]),points[:,2])[1]
bestPoint = points[dominantInd,:]
