const N = 10
x = linspace(0, 1, N+1)
writedlm("table.txt", [x besselj(0,x)], ' ')
