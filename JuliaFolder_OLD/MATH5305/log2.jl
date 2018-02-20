# This is a function to calculate the digits of the log function

p = 0.5
q = 0.5
r = 0.25
iter = 1000

x = [16 4 4 1]
val = [11 0 0 0]

v(x) = max(x-5,0)

sum = 0
bigsum = 0
for j=1:iter
  for i=1:iter
    num = round(Int64,floor(4*rand(1)+1))
    sum += v(x[num])
    #println(sum)
  end
  sum /= iter
  bigsum += sum
end
v0 = bigsum/(iter*((1+r)^2))
#@printf('The average value is: %f\n',sum)
println("The average value is: $v0")
