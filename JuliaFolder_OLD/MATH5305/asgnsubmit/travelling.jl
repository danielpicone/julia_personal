include("AllenCahn.jl")
using AllenCahn
using PyPlot
using OffsetArrays

const Ω = (-0.5, 1.5)
const ϵ = 0.015
const c = 3 / (sqrt(2)*ϵ)
const T = 1/c
const P = 256
const N = 8P

gr = grid1D(Ω, T, P, N)

function u(x, t)
    c = sqrt(2)*ϵ
    y = ( x - 3t/c ) / (2c)
    return ( 1 - tanh(y) ) / 2
end

x, t = gr.x, gr.t
u0 = u.(x, 0.0)
U = semi_implicit(ϵ, gr, u0)

idx = 0:div(N,4):N
figure(1)
plot(x[1:P], u.(x[1:P], gr.t[idx]'))
xlabel(L"$x$")
grid(true)
legend((L"$t=0$", L"$t=\frac{1}{4}T$", L"$t=\frac{1}{2}T$",
                  L"$t=\frac{3}{4}T$", L"$t=T$"))
savefig("travelling-with-border.png")
run(`convert travelling-with-border.png -trim travelling.png`)
run(`rm travelling-with-border.png`)


figure(2)
plot(x[1:P], U[1:P,idx])
legend((L"$t=0$", L"$t=\frac{1}{4}T$", L"$t=\frac{1}{2}T$",
                  L"$t=\frac{3}{4}T$", L"$t=T$"))

err_ = maxabs(U[1:P,0:N] - u.(x[1:P], gr.t[0:N]'),1)
err = OffsetArray(err_[:], 0:N)
figure(3)
idx = 0:div(N,100):N
plot(t[idx], err[idx])
xlabel(L"$t$")
ylabel("Max error")
@printf("Max error = %10.3e\n", maximum(err))
