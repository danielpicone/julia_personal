using Elliptic
using OffsetArrays
using PyPlot

const P0 = 8
const Q0 = 8
const its = 3

include("problem1.jl")
mg = multigrid(Lx, Ly, P0, Q0, 1)
coarse = mg[0].gr
fine   = mg[1].gr
ω = Richardson_parameter(fine, :smoother)

smoother!() = Richardson!(mg[1].U, mg[1].r, 3, ω, fine, mg[1].b)
lo = [-0.1, -0.05, -0.005]
hi = [2.0, 0.2, 0.03]

#smoother!() = rbGaussSeidel!(mg[1].U, 3, fine, mg[1].b)
#lo = [-0.1, -0.2, -0.005]
#hi = [2.0, 0.05, 0.015]


P, Q = coarse.P, coarse.Q
insert_boundary_values!(mg[1].U, g, fine)
for q = 1:2Q-1, p = 1:2P-1
    mg[1].b[p,q] = f(fine.x[p], fine.y[q])
end
residual!(mg[1].r, mg[1].b, mg[1].U, fine)

Af = band_Poisson_matrix(fine)
band_Cholesky!(Af)
Mf = (2P-1)*(2Q-1)
Vf = Poisson_RHS(fine, f, g)
band_Cholesky_solve!(Af, Vf)
exactUf = copy(mg[1].U)
exactUf[1:2P-1,1:2Q-1] = reshape(Vf, 2P-1, 2Q-1)

figure(0, figsize=(13,8))

subplot(2, 2, 1, projection="3d")
surf(fine.x[0:2P], fine.y[0:2Q], exactUf[0:2P,0:2Q]', 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x", fontsize=14)
ylabel(L"y", fontsize=14)
title(L"U_{\mathrm{f}}")

subplot(2, 2, 2, projection="3d")
surf(fine.x[0:2P], fine.y[0:2Q], mg[1].U[0:2P,0:2Q]', 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x", fontsize=14)
ylabel(L"y", fontsize=14)
title(L"U^0_{\mathrm{f}}")

subplot(2, 2, 3, projection="3d")
ϵ = exactUf - mg[1].U
surf(fine.x[0:2P], fine.y[0:2Q], ϵ[0:2P,0:2Q]',
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x", fontsize=14)
ylabel(L"y", fontsize=14)
title(L"U-U^0_{\mathrm{f}}")

subplot(2, 2, 4, projection="3d")
surf(fine.x[1:2P-1], fine.y[1:2Q-1], (mg[1].r)', 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x", fontsize=14)
ylabel(L"y", fontsize=14)
title(L"r^0_{\mathrm{f}}")

Mc = (P-1)*(Q-1)

for k = 1:its
    figure(k, figsize=(19,8))

    smoother!()
    ϵ1 = exactUf - mg[1].U
    # Error after pre-smoothing
    subplot(2, 3, 1, projection="3d")
    surf(fine.x[0:2P], fine.y[0:2Q], ϵ1[0:2P,0:2Q]',
         rstride=1, cstride=1, cmap="autumn")
    s = latexstring("U_{\\mathrm{f}}-U_{\\mathrm{f}}^{$(k-1)+1/3}")
    title(s)
    zlim(lo[k], hi[k])

    # Fine-grid residual
    residual!(mg[1].r, mg[1].b, mg[1].U, fine)
    subplot(2, 3, 2, projection="3d")
    surf(fine.x[1:2P-1], fine.y[1:2Q-1], (mg[1].r)',
         rstride=1, cstride=1, cmap="autumn")
    s = latexstring("r_{\\mathrm{f}}^{$(k-1)+1/3}")
    title(s)

    # Coarse-grid residual
    restrict!(mg[0].b, mg[1].r)
    subplot(2, 3, 3, projection="3d")
    surf(coarse.x[1:P-1], coarse.y[1:Q-1], (mg[0].b)',
         rstride=1, cstride=1, cmap="autumn")
    s = latexstring("r_{\\mathrm{c}}^{$(k-1)+1/3}")
    title(s)

    # Coarse-grid approximation to error
    band_Cholesky_solve!(mg.A0, vec(mg[0].b))
    mg[0].U[1:P-1,1:Q-1] .= mg[0].b
    ϵ2 = copy(mg[0].U)
    subplot(2, 3, 4, projection="3d")
    surf(coarse.x[0:P], coarse.y[0:Q], ϵ2[0:P,0:Q]',
         rstride=1, cstride=1, cmap="autumn")
    zlim(lo[k], hi[k])
    title(L"$V_{\mathrm{c}}$")

    # Error after coarse-grid correction.
    correct!(mg[1].U, mg[0].U)
    ϵ3 = exactUf - mg[1].U
    subplot(2, 3, 5, projection="3d")
    surf(fine.x[0:2P], fine.y[0:2Q], ϵ3[0:2P,0:2Q]',
         rstride=1, cstride=1, cmap="autumn")
    zlim(lo[k], hi[k])
    s = latexstring("U_{\\mathrm{f}}-U_{\\mathrm{f}}^{$(k-1)+2/3}")
    title(s)

    # Error after post-smoothing
    smoother!()
    ϵ4 = exactUf - mg[1].U
    subplot(2, 3, 6, projection="3d")
    surf(fine.x[0:2P], fine.y[0:2Q], ϵ4[0:2P,0:2Q]',
         rstride=1, cstride=1, cmap="autumn")
    zlim(lo[k], hi[k])
    s = latexstring("U_{\\mathrm{f}}-U_{\\mathrm{f}}^{$k}")
    title(s)
end
