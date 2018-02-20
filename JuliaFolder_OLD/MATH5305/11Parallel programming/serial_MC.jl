using MC

const trials = 50_000
const L = 2.0
const s = 1.2
const C = 2 / zeta(s)
const M = 30
const P = 4M
const γ = 0.5
const ω = Array(Float64, M)
const A = SymTridiagonal(zeros(P), zeros(P-1))
const b = Array(Float64, P)
const U = Array(Float64, P)
const u0 = Array(Float64, trials)

κ = diffusivity(L, s, C)
x, Δx = gridpoints(L, P)

@printf("Serial Monte-Carlo with %d trials.\n", trials)
@printf("\t%d random variables and %d finite difference grid points.\n", 
        M, P)

start = time()
Eu0 = 0.0
for n = 1:trials
    rand!(ω)
    for m = 1:M
        ω[m] -= 1/2
    end
    u0[n] = u_at_zero!(A, b, x -> κ(x, ω), γ, Δx)
    Eu0 += u0[n]
end
Eu0 /= trials
finish = time()

σκ = stdκ(L, s)
σu0 = std(u0; mean=Eu0)

@printf("\tStandard deviation in κ(x) = %0.4f\n", σκ)
@printf("\tu(0) = %0.5f ± %0.5f\n", Eu0, σu0)
@printf("\tElapsed time = %0.2f seconds\n", finish-start)
