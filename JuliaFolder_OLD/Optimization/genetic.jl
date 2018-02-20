# This is a script which implements a genetic algorithm to solve an unconstrained optimisation problem
# This is for the single dimension case
dim = 1
# f(x,y) = x^2+y^2
# f(x) = x^2+1
f(x) = x^4 -3*x^3 - 2*x^2 + 6*x

# Number of points which can breed
numIndividuals = 10
numBreed = 5

function breed!(dominantPoint,otherPoint)
    # Set a 1% chance for mutation of dominantPoint
    # 25% chance for average of parents
    # 25% chance for half of each parent
    # 25% chance for dom - other
    # 25% chance for dom + other
    if rand() < 0.01
        otherPoint = dominantPoint
        otherPoint[1:end-1] = dominantPoint[1:end-1]+randn()
        otherPoint[end] = f(otherPoint[1:end-1]...)
        println("We mutated!!")
        return otherPoint
    else
        chance = rand()
        if chance <= 0.333333
            otherPoint = (dominantPoint+otherPoint)./2
            otherPoint[end] = f(otherPoint[1:end-1]...)
            return otherPoint
        elseif chance <= 0.666666
            otherPoint = dominantPoint-otherPoint
            otherPoint[end] = f(otherPoint[1:end-1]...)
            return otherPoint
        else
            otherPoint = dominantPoint+otherPoint
            otherPoint[end] = f(otherPoint[1:end-1]...)
            return otherPoint
        end
    end
    return otherPoint
end




# Choose the point with lowest obj value and let it breed with the numBreed
points = Array(Float64,numIndividuals,2)

for j=1:numIndividuals
    points[j,1:dim] = 100*rand(1,dim)
    points[j,dim+1] = f(points[j,1]...)
end

iterations = 1000
for i=1:iterations
    # Choose numBreed points at random
    # Indicies of breeding points
    ind = randperm(numIndividuals)[1:numBreed]
    # Find the index of the dominant point
    dominantInd = find(x->x==minimum(points[ind,2]),points[:,2])[1]
    # breedPoints = sortrows(points[ind,:],by = x->x[end])

    # Create the dominant point which has the lowest obj value
    # dominantPoint = breedPoints[1,:]

    for j in ind
        # Now pairwise breed between the dominant and the other points
        if j==dominantInd
        else
            points[j,:] = breed!(points[dominantInd,:],points[j,:])
        end

    end

end
dominantInd = find(x->x==minimum(points[:,2]),points[:,2])[1]
bestPoint = points[dominantInd,:]
