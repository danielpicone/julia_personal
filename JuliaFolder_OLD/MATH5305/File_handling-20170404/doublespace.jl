open(infile -> begin
                   for line in eachline(infile)
                       println(line)
                   end
               end,
     "load_vars.jl", "r")
