# Estimating pi by throwing darts at a board
using PyPlot

n = 10000
points = 2*rand(2,n)-ones(2,n)

incircle = 0
for i=1:n
  if points[1,i]^2+points[2,i]^2<=1
    incircle += 1
  end
end
x = linspace(-1,1,10000)
func1 = sqrt(1-x.^2)
func2 = -sqrt(1-x.^2)

figure(1)
plot(x,func1,"black",x,func2,"black",points[1,:],points[2,:],"bo",markersize=0.1,x,-1*ones(10000,1),"black",x,1*ones(10000,1),"black",-1*ones(10000,1),x,"black",1*ones(10000,1),x,"black")
axis("equal")
ax = gca()
ax[:set_xlim]([-1,1])
ax[:set_ylim]([-1,1])


est = 4*incircle/n
@printf("The value of pi is: %f",est)
