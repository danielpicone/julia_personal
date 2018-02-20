using Elliptic
using OffsetArrays
using PyPlot

const P = 20
const Q = 20
steps = 2

gr = rectgrid(1.0, 1.0, P, Q)
ω = Richardson_parameter(gr, :smoother)

f = zeros(P-1,Q-1)
r = zeros(P-1,Q-1)

ϵ = OffsetArray(zeros(P+1,Q+1), 0:P, 0:Q)
ϵ0 = rand(P-1, Q-1)

function error_plots(smoother!, figno, filename)
    ϵ[1:P-1,1:Q-1] = ϵ0
    figure(figno)
    subplot(2, 2, 1, projection="3d")
    surf(gr.x[0:P], gr.y[0:Q], ϵ[0:P,0:Q]', 
         rstride=1, cstride=1, cmap="autumn")
    xlabel(L"x")
    ylabel(L"y")
    title("Initial error")

    ν = 0
    for j = 1:3
        smoother!()
        ν += steps
        subplot(2, 2, j+1, projection="3d")
        surf(gr.x[0:P], gr.y[0:Q], ϵ[0:P,0:Q]', 
             rstride=1, cstride=1, cmap="autumn")
        zlim(0,1)
        xlabel(L"x")
        ylabel(L"y")
        title("After $ν iterations")
    end
    savefig("temp.png")
    run(`convert temp.png -trim $filename.png`)
end

error_plots(() -> Richardson!(ϵ, r, steps, ω, gr, f),
            "Richardson", "Richardson_smoothing")
error_plots(() -> rbGaussSeidel!(ϵ, steps, gr, f),
            "Red-black Gauss-Seidel", "rbGaussSeidel_smoothing")
