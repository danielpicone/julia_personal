# Calculate the Fourier series of a function

using PyPlot

# Define the function which is periodic with period L

L = π
N = 10
points = 500
simpsonPoints = 50

function f(x,L)
  if -L <= x <= L
    #return x.^2
  elseif x > L
    while x > L
      x -= 2*L
    end
    #return x/π
  elseif x < -L
    while x < -L
      x += 2*L
    end
    #return x/π
  end
  return x.^2/π^2
  #return x/π
  """if 0<= x <= L
    return 1
  else
    return -1
  end"""
end

function simpsonRule(L,coeff,n,x)
  simpsonSum = 0
  Δx = 2*L/simpsonPoints
  if coeff == "cos"
    simpsonSum += Δx/3 * (f(x[1],L).*cos(n*x[1]) + 4f(x[2],L).*cos(n*x[2]) + f(x[3],L).*cos(n*x[3]))
    if simpsonPoints >= 6
      for i=4:2:length(x)-2
        simpsonSum += Δx/3 * (f(x[i-1],L).*cos(n*x[i-1]) + 4f(x[i],L).*cos(n*x[i]) + f(x[i+1],L).*cos(n*x[i+1]))
      end
    end
    simpsonSum += Δx/3 * (f(x[simpsonPoints-1],L).*cos(n*x[simpsonPoints-1]) + 4f(x[simpsonPoints],L).*cos(n*x[simpsonPoints]) + f(x[simpsonPoints+1],L).*cos(n*x[simpsonPoints+1]))
  elseif coeff == "sin"
    simpsonSum += Δx/3 * (f(x[1],L).*sin(n*x[1]) + 4f(x[2],L).*sin(n*x[2]) + f(x[3],L).*sin(n*x[3]))
    if simpsonPoints >= 6
      for i=4:2:length(x)-2
        simpsonSum += Δx/3 * (f(x[i-1],L).*sin(n*x[i-1]) + 4f(x[i],L).*sin(n*x[i]) + f(x[i+1],L).*sin(n*x[i+1]))
      end
    end
    simpsonSum += Δx/3 * (f(x[simpsonPoints-1],L).*sin(n*x[simpsonPoints-1]) + 4f(x[simpsonPoints],L).*sin(n*x[simpsonPoints]) + f(x[simpsonPoints+1],L).*sin(n*x[simpsonPoints+1]))
  end
  return simpsonSum/π
end

function an(L,x,n)
  return 0
end

function bn(L,x,n)
  return 2*(-1)^(n+1)/(π*n)
end

function fourier(x)
  areax = linspace(-L,L,simpsonPoints+1)
  #s = an(L,x,0)/2
  s = simpsonRule(L,"cos",0,areax)*π
  #s += sum(an(L,x,n)*cos(n*x) for n=1:N)
  s += sum(simpsonRule(L,"cos",n,areax)*cos(n.*x) for n=1:N)
  s += sum(simpsonRule(L,"sin",n,areax)*sin(n.*x) for n=1:N)
  #s += sum(bn(L,x,n)*sin(n*x) for n=1:N)
  return s
end

xaxis = linspace(-4π,4π,points)
#fourier.(xaxis)
#plot(xaxis,f.(xaxis,L),xaxis,fourier(xaxis)[1:end-1])
plot(xaxis,f.(xaxis,L),xaxis,fourier.(xaxis))
