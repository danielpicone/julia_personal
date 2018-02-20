# Taylor series

using PyPlot

function taylor(x,N)
  a = 1
  func = a
  err = 0
  for i = 1:N
    a *= x/(i)
    func += a
    err = abs(func-exp(x))
  end
  return func,err
end

"""k = zeros(points)
xaxis = linspace(0,10,points)
for i=1:points
  k[i] = taylor(1,xaxis[i],Nk)
end"""

x = 0.1
N = 20
err = zeros(N)
tol = convert(Float64,10)^(-10)
funcval,err = taylor(u,N)

iter = 0

@printf("Iteration    Taylor Series       Error\n\n")
for i = 1:N
  funcval,err = taylor(x,i)
  if (N <= 50)&&(i!=N)
    @printf("%3d          %3.10f        %1.10e\n",i,funcval,err)
  end
  if (N > 50)&&(i!=N)
    if i%10 == 0
      @printf("%3d          %3.10f        %1.10e\n",i,funcval,err)
    end
  end
  if i==N
    @printf("%3d          %3.10f        %1.10e    %3.10f\n",i,funcval,err,e^u)
  end
  if err<tol
    break
  end
  iter += 1
end


#plot(xaxis,k,xaxis,l)
