using JLD

a = rand(1000)
b = rand(100,100)
c = sprand(1000,1000, 0.01)

save("mydata.jld", "a", a, "b", b, "c", c)
