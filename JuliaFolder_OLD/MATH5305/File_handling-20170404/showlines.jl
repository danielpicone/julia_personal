function showlines(filename::String)
    infile = open(filename, "r")
    for line in eachline(infile)
        write(STDOUT, line)
    end
    close(infile)
    return nothing
end
