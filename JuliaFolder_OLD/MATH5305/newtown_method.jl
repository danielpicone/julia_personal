# Newtown's method of finding roots
N = 20
coeff = [1 6 -7 -4 3]

using PyPlot

x = 1.0*ones(N,1)
x_err = 1.0*ones(N,1)
x[1] = -1
tol = 0.00000001
iter = 1


for j = 1:N-1
  x[j+1] = x[j]-f(x[j],coeff)/g(x[j],coeff)
  x_err[j+1] = x[j+1]-x[j]
  iter += 1
  if abs(x_err[j+1]) <= tol
    break
  end
end

@printf("  j  x               f(x)               error\n")
@printf("%3d  %3.12f  %3.12f\n",1,x[1],f(x[1],coeff))
for j = 2:iter
  @printf("%3d  %3.12f  %3.12f    %3.12f\n",j,x[j],f(x[j],coeff),abs(x_err[j]))
end

#plot([1,2,3,4,5,6],[8,2,9,2,4,3],"rD")

#equation = ones(N,1)
#for i=1:10
#  equation = f(i,coeff)
#end
xaxis = linspace(-1.1,1,100)
plot(xaxis,f(xaxis,coeff),x,f(x,coeff),"D",markersize="3")
#plot(x,f(x,coeff),'rD')
grid(true)
#plot(1,2)
