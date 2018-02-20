num_points = 20

theta = linspace(-π/18,π+π/18,num_points)

sum = 0

for i = 1:num_points
  for j = i+1:num_points
    sum += 1/(theta[i]-theta[j])^2
  end
end
println(sum)
