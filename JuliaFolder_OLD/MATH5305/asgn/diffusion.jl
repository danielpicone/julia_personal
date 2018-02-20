# Answer to quesiton 8
# Edited by Daniel Picone z3459589

include("AllenCahn.jl")

using AllenCahn
using PyPlot
using OffsetArrays

const a = 0.0
const b = 1.0
const c = 0.0
const d = 1.0
const Ω = RectangularDomain(a, b, c, d)
const ϵ = Inf
const k = 2
const l = 1
const λ = (k*pi/(b-a))^2 + (l*pi/(d-c))^2
const T = 2.0 / λ

phi(x,y) = cos(k*π*(x-a)/(b-a)) * cos(l*π*(y-c)/(d-c))
function exactu(x, y, t)
    λ = (k*pi/(b-a))^2 + (l*pi/(d-c))^2
    return exp(-λ*t) * phi(x, y)
end

P = 64
Q = 64
N = 512

gr = grid2D(Ω, T, P, Q, N)
x, y, t = gr.x, gr.y, gr.t
#u0 = phi.(x, y')
u0 = Float64[phi(x[p],y[q]) for p=1:P, q=1:Q]
U = semi_implicit(ϵ, gr, u0)

figure(1)
#plot_surface(gr.x, gr.y, U[:,:,end]')
mesh(x, y, U[1:P,1:Q,N]')
xlabel(L"$x$")
ylabel(L"$y$")

figure(2)
mesh(gr.x, gr.y, exactu.(x,y',T)')
xlabel(L"$x$")
ylabel(L"$y$")

figure(3)
Uerr = U[:,:,N]-exactu.(x,y',T)
contourf(x, y, Uerr[1:P,1:Q]')
xlabel(L"$x$")
ylabel(L"$y$")
colorbar()

P = Int64(32/2)
Q = Int64(32/2)
N = Int64(128/2)
err = Float64[]
@printf("%5s  %5s  %5s  %12s  %8s  %8s\n\n",
        "P", "Q", "N", "Max error", "rate", "seconds")
for row = 1:5
    start = time()
    gr = grid2D(Ω, T, P, Q, N)
    x, y, t = gr.x, gr.y, gr.t
    u0 = phi.(x, y')
    U = semi_implicit(ϵ, gr, u0)
    push!(err, maxabs(U[1:P,1:Q,N]-exactu(x,y',T)))
    finish = time()
    elapsed = finish - start
    if row == 1
        @printf("%5d  %5d  %5d  %12.3e  %8s  %8.4f\n",
                P, Q, N, err[row], "", elapsed)
    else
        rate = log2(err[row-1]/err[row])
        @printf("%5d  %5d  %5d  %12.3e  %8.4f  %8.4f\n",
                P, Q, N, err[row], rate, elapsed)
    end
    P *= 2
    Q *= 2
    N *= 2
end
