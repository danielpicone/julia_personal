using Elliptic
using OffsetArrays

include("problem1.jl")

const P0 = 32
const Q0 = 32
ν = (1, 2)
iterations = 15

mg = multigrid(Lx, Ly, P0, Q0, 1)
coarse = mg[0].gr
fine   = mg[1].gr
P, Q = coarse.P, coarse.Q

name = "Richardson"
#name = "Red-black Gauss-Seidel"

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

println("\nTwo-grid iteration.")
println("\t$name smoother on $(2P)-by-$(2Q) fine grid.")
println("\tUsing ν pre-smoothing and ν post-smoothing steps.")
println("\tBand Cholesky solve on $P-by-$Q coarse grid.\n")

exactU = OffsetArray(Float64, 0:2P, 0:2Q)
for q = 0:2Q, p = 0:2P
    exactU[p,q] = exact_u(fine.x[p], fine.y[q])
end
for q = 1:2Q-1, p = 1:2P-1
    mg[1].b[p,q] = f(fine.x[p], fine.y[q])
end

resid = zeros(iterations, 2)
err = zero(resid)

for case = 1:2
    fill!(mg[1].U, 0.0)  # initial guess
    insert_boundary_values!(mg[1].U, g, fine)
    for k = 1:iterations
        smooth!(1, ν[case], mg)
        residual!(mg[1].r, mg[1].b, mg[1].U, fine)
        restrict!(mg[0].b, mg[1].r)
        band_Cholesky_solve!(mg.A0, vec(mg[0].b))
	for q=1:Q-1, p=1:P-1
            mg[0].U[p,q] = mg[0].b[p,q]
        end
        correct!(mg[1].U, mg[0].U) 
	smooth!(1, ν[case], mg)
        residual!(mg[1].r, mg[1].b, mg[1].U, fine)
	ϵ = exactU - mg[1].U
        err[k,case] = L2norm(ϵ, fine)
        resid[k,case] = L2norm(mg[1].r, fine)
    end
end

@printf("     |           %10s            |           %10s\n", 
        "ν = $(ν[1])", "ν = $(ν[2])")
@printf("%4s |%10s  %10s  %8s |%10s  %10s  %8s\n", 
        "", "", "", "", "", "", "")
@printf("%4s |%10s  %10s  %8s |%10s  %10s  %8s\n", 
        "k", "||ϵₖ||₂", "||rₖ||₂", "ratio", "||ϵₖ||₂", "||rₖ||₂", "ratio")
@printf("%4s |%10s  %10s  %8s |%10s  %10s  %8s\n", 
        "", "", "", "", "", "", "")
@printf("%4d |%10.2e  %10.2e  %8s |%10.2e  %10.2e  %8s\n", 
        1, err[1,1], resid[1,1], "", err[1,2], resid[1,2], "")
for k = 2:iterations
    ratio = Float64[resid[k,case] / resid[k-1,case] for case = 1:2]
    @printf("%4d |%10.2e  %10.2e  %8.4f |%10.2e  %10.2e  %8.4f\n", 
            k, err[k,1], resid[k,1], ratio[1], err[k,2], resid[k,2], ratio[2])
end
