using Elliptic

include("problem1.jl")

const P0 = 32
const Q0 = 32
const L = 5


function smooth!(l, ν, mg)
    ω = Richardson_parameter(mg[l].gr, :smoother)
    Richardson!(mg[l].U, mg[l].r, ν, ω, mg[l].gr, mg[l].b)
end

function mgsolver(maxits=100, rtol=1.0e-8)
    mg = setup_vcycle(Lx, Ly, P0, Q0, L, f, g)
    resid = Float64[]
    for k = 1:maxits
        vcycle!(mg, smooth!; presmooth=2, postsmooth=2)
        push!(resid, L2norm(mg[L].r, mg[L].gr))
        if resid[k] < rtol
            break
        end
    end
    return resid
end
