# This is a program to estimate the area under a curve
# It will use both Riemann sums and the trapezoidal method
module AreaFind
#using PyPlot

export riemannSum,trapezoidalSum,simpsonSum,MCIntegration
# This will use the left Riemann sum to find the area
function riemannSum(a,b,N,f)
  L = b - a
  Δx = L/N
  x = linspace(a,b,N+1)
  riemannSum = 0
  for i=1:N+1
    riemannSum += f(x[i])*Δx
  end
  #@printf("The left Riemann area is:        %f\n",riemannSum)
  return riemannSum
end

# This will use the trapezoidal method to find the area
function trapezoidalSum(a,b,N,f)
  L = b - a
  Δx = L/N
  x = linspace(a,b,N+1)
  trapezoidalSum = 0
  for i=1:N
    trapezoidalSum += 0.5*Δx*(f(x[i])+f(x[i+1]))
  end
  #@printf("The trapezoidal area is:         %f\n",trapezoidalSum)
  return trapezoidalSum
end

# This will use simpson's rule to find the area
function simpsonSum(a,b,N,f)
  L = b - a
  Δx = L/N
  x = linspace(a,b,N+1)
  simpsonSum = 0
  simpsonSum += Δx/3 * (f(x[1]) + 4f(x[2]) + f(x[3]))
  if N >= 6
    for i=4:2:N-2
      simpsonSum += Δx/3 * (f(x[i-1]) + 4f(x[i]) + f(x[i+1]))
    end
  end
  simpsonSum += Δx/3 * (f(x[N-1]) + 4f(x[N]) + f(x[N+1]))
  @printf("The simpson area is:             %f\n",simpsonSum)
  return simpsonSum
end

function MCIntegration(a,b,N,f)
  x = linspace(a,b,10000)
  """maxFunc = -Inf
  minFunc = Inf
  for i=1:10000
    maxFunc = max(maxFunc,f(x[i]))
    minFunc = min(minFunc,f(x[i]),0)
  end"""
  maxFunc = maximum(f.(x))
  minFunc = minimum(f.(x))
  h = maxFunc - minFunc
  L = abs(b-a)
  Area = h*L
  randX = L*rand(N) + min(a,b)
  randY = h*rand(N) + min(f(a),f(b),0)
  below = 0
  for i=1:N
    if randY[i] <= f(randX[i])
      below += 1
    end
    #MCIntegration += f(L*rand() + min(a,b))
  end
  #MCIntegration = (MCIntegration/N)*rectangleArea
  #MCIntegration = L*MCIntegration/N
  #println(below/N)
  #@printf("h = %f     L = %f    \n",h,L)
  #println(Area)
  MCIntegration = (below/N)*Area
  @printf("The Monte Carlo Area is:       %f\n",MCIntegration)
  #plot(randX,randY,x,f.(x))
  #plot(randX,randY)
  return MCIntegration,randX,randY,x
end

end
