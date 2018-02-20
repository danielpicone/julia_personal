using Elliptic
using OffsetArrays

include("problem1.jl")

const P0 = 32
const Q0 = 32
const L = 5
const maxits = 40
const ΔU_tol = 1.0e-9
const r_tol = 1.0e-14

name = "Richardson"
#name = "Red-black Gauss-Seidel"

mg = setup_vcycle(Lx, Ly, P0, Q0, L, f, g)
fine = mg.greq[L].gr
x, y, P, Q = fine.x, fine.y, fine.P, fine.Q
exactU = OffsetArray(Float64, 0:P, 0:Q)
for q = 0:Q, p = 0:P
    exactU[p,q] = exact_u(x[p], y[q])
end

function Richardson_smooth!(l, ν, mg)
    ω = Richardson_parameter(mg[l].gr, :smoother)
    Richardson!(mg[l].U, mg[l].r, ν, ω, mg[l].gr, mg[l].b)
end

rbGaussSeidel_smooth!(l, ν, mg) = rbGaussSeidel!(mg[l].U, ν, mg[l].gr, mg[l].b)

if name == "Richardson"
    smooth! = Richardson_smooth!
else
    smooth! = rbGaussSeidel_smooth!
end



println("\nMultigrid V-cycle solver with $(L+1) levels\n")
println("\tUsing $name smoother.")
P0, Q0 = mg[0].gr.P, mg[0].gr.Q
println("\tLevel 0: $P0-by-$Q0.")
println("\tLevel $L: $P-by-$Q\n")

@printf("%4s  %10s  %8s  %10s  %8s\n\n",
        "k", "||ΔUₖ||₂", "ratio", "||rₖ||₂", "ratio")

ΔU = Float64[]
resid = Float64[]
prev_U = OffsetArray(Float64, 0:P, 0:Q)
for k = 1:maxits
    copy!(prev_U, mg[L].U)
    vcycle!(mg, smooth!; presmooth=3, postsmooth=3)
    err = prev_U - mg[L].U
    push!(ΔU, L2norm(err, fine))
    residual!(mg[L].r, mg[L].b, mg[L].U, fine)
    push!(resid, L2norm(mg[L].r, mg[L].gr))
    if k == 1
        @printf("%4d  %10.2e  %8s  %10.2e\n", 1, ΔU[1], "", resid[1])
    else
        ΔU_ratio = ΔU[k] / ΔU[k-1]
        r_ratio  = resid[k] / resid[k-1]
        @printf("%4d  %10.2e  %8.4f  %10.2e  %8.4f\n", 
                k, ΔU[k], ΔU_ratio, resid[k], r_ratio)
    end
    if ΔU[k] < ΔU_tol || resid[k] < r_tol
        break
    end
end

err = exactU - mg[L].U
errnrm = L2norm(err[0:P,0:Q], mg[L].gr)
@printf("\n\nFinal error || U - u ||₂ = %10.2e\n", errnrm)
