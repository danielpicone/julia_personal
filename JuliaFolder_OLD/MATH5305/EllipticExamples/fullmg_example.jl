using Elliptic

include("problem1.jl")

const P0 = 8
const Q0 = 8
const L = 7
const maxits = 100
const ΔU_tol = 1.0e-9
const r_tol = 1.0e-9

mg = setup_fullmg(Lx, Ly, P0, Q0, L, f, g)
fine = mg.greq[L]
x, y = fine.gr.x, fine.gr.y
exactU = exact_u.(x, y')

println("\nFull multigrid solver with $(L+1) levels.\n")
P0, Q0 = mg.gr0.P, mg.gr0.Q
println("\tLevel 0: $P0-by-$Q0.")
P, Q = fine.gr.P, fine.gr.Q
println("\tLevel $L: $P-by-$Q\n.")

fullmg!(mg)

ϵ = exactU - fine.U
err = L2norm(ϵ, fine.gr)
@printf("\n\nFull mg error || U - u ||₂ = %10.2e\n", err)
residual!(fine.r, fine.b, fine.U, fine.gr)
fmg_resid = L2norm(fine.r, fine.gr)
@printf("\nFull mg residual || r ||₂ = %10.2e\n", fmg_resid)

@printf("\n%4s  %10s  %10s  %8s\n\n",
        "k", "||error||₂", "||rₖ||₂", "ratio")

ΔU = Float64[]
resid = Float64[]

for k = 1:maxits
    vcycle!(mg)
    residual!(fine.r, fine.b, fine.U, fine.gr)
    ϵ .= exactU - fine.U
    push!(ΔU, L2norm(ϵ, fine.gr))
    push!(resid, L2norm(fine.r, fine.gr))
    if k == 1
        @printf("%4d  %10.2e  %10.2e\n", 
	        1, ΔU[1], resid[1])
    else
        r_ratio  = resid[k] / resid[k-1]
        @printf("%4d  %10.2e  %10.2e  %8.4f\n", 
                k, ΔU[k], resid[k], r_ratio)
    end
    if ΔU[k] < ΔU_tol || resid[k] < r_tol
        break
    end
end

err = L2norm(fine.U-exactU, fine.gr)
@printf("\n\nFinal error || U - u ||₂ = %10.2e\n", err)
