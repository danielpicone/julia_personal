# Writes a table to standard output

const N = 10
x = linspace(0, 1, N+1)
@printf("%5s  %18s\n", "#   x", "J₀(x)")
for j = 1:N+1
    #@printf("%5.2f  %18.14f\n", x[j], besselj(0,x[j]))
    println(x[j], ' ', besselj(0,x[j]))
end
