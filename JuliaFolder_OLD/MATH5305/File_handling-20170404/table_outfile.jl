const N = 10
open("table.txt", "w") do outfile
    x = linspace(0, 1, N+1)
    @printf(outfile, "%5s  %18s\n", "#   x", "J₀(x)")
    for j = 1:N+1
        @printf(outfile, "%5.2f  %18.14f\n", x[j], besselj(0,x[j]))
    end
end
