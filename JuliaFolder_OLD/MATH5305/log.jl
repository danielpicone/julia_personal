# This is a function to calculate the digits of the log function

using PyPlot

N = 10
x = 0.2
xsqr = x^2
a = x
func = a
func_err = ones(N-1,1)
actual = log(1-x)
start = time()
"""for i = 1:N-1
  #@printf("%1.16f  %f\n",a,func)
  a *= x*i/(i+1)
  func += a
  func_err[i] = actual+func
end"""
stop = time()
#println("The value of the function is ",-func)
#println("The value of the function is ",-taylorlog(x,N))
println("The value of Julia function is ", actual)

@printf("  N   Error\n\n")
for i = 1:N-1
  #@printf("%3d   %1.16f\n", i, abs(func_err[i]))
end
println(stop-start)

xaxis = linspace(-2.25,0.2,100)'
#xaxis = [0.2 0.3 0.4]
#plot(xaxis,log(1-xaxis),xaxis',taylorlog(xaxis,N))
println(xaxis)
p = taylorlog(xaxis',N)
#println(length(p),length(xaxis))
plot(xaxis',log(1-xaxis)',xaxis',p,"g")
