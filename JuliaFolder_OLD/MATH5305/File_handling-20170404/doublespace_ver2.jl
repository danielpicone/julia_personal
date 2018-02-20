open("load_vars.jl", "r") do infile
    for line in eachline(infile)
        println(line)
    end
end
