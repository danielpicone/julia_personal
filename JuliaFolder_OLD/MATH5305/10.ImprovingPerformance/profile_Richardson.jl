using Elliptic
using OffsetArrays

include("problem1.jl")

const P = 1024
const Q = 1024
gr = rectgrid(Lx, Ly, P, Q)
ω = Richardson_parameter(gr, :solver)
U = OffsetArray(Float64, 0:P, 0:Q)
insert_boundary_values!(U, g, gr)
r = zeros(P-1,Q-1)
b = f.(gr.x[1:P-1], gr.y[1:Q-1]')

function iterate(func::Function, maxits=50, rtol=1.0e-8)
    U[1:P-1,1:Q-1] = 0.0
    resid = Float64[]
    for k = 1:maxits
        func(U, r, 10, ω, gr, b)
        push!(resid, L2norm(r, gr))
        if resid[k] < rtol
            break
        end
    end
    return resid
end
