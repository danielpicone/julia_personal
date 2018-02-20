# Hill climbing algorithm
using PyPlot

#f(x) = (x-3).^2
#f(x,y) = x.^2 - 2.*x.*y + y.^2 - sin.(x)
f(x,y) = x.^2 + 4y.^2+3x - sin.(x)
(x,y) = (3.0,1.0)

μ0 = 100
ρ = 0.8

Δx = μ0
Δy = μ0
iternum = 0
#while (f(x+Δx,y) <= f(x,y))||(f(x-Δx,y) <= f(x,y)||f(x,y+Δy) <= f(x,y))||(f(x,y-Δy) <= f(x,y))
for i=1:1000
  if (iternum % 2 == 1)
    if (f(x+Δx,y) < f(x,y))
      x += Δx
    else (f(x-Δx,y) < f(x,y))
      x -= Δx
    end
    Δx = μ0 * ρ^iternum
  elseif (iternum % 2 == 0)
    if (f(x,y+Δy) < f(x,y))
      y += Δy
    else (f(x,y-Δy) < f(x,y))
      y -= Δy
    end
    Δy = μ0 * ρ^iternum
  end
  println("x value is:",x)
  println("y value is:",y)
  iternum += 1
  if iternum > 100000
    break
  end
end

@printf("The local minimum is:  %f\n",f(x,y))
@printf("Achieved at:           (%f,%f)\n",x,y)
@printf("In %d iterations\n",iternum)

xaxis = linspace(Int64(round(x))-10,Int64(round(x))+10,10000)
yaxis = linspace(Int64(round(y))-10,Int64(round(y))+10,10000)
#figure(1)
#plot3D(xaxis,yaxis,f(xaxis,yaxis))
#plot(xaxis,f(xaxis))
