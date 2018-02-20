using Elliptic
using OffsetArrays
using PyPlot

const P0 = 8
const Q0 = 8
const Lx = 1.0
const Ly = 1.0

mg = multigrid(Lx, Ly, P0, Q0, 1)
coarse = mg[0].gr
fine   = mg[1].gr
v(x, y) = x * y + sin(2π*x) * cos(3π*y)

P, Q = coarse.P, coarse.Q
for q = 0:Q, p=0:P
    mg[0].U[p,q] = v(coarse.x[p], coarse.y[q])
end
PV = OffsetArray(Float64, 0:2P, 0:2Q)
fill!(mg[1].U, 0.0)
insert_boundary_values!(mg[1].U, v, fine)
correct!(mg[1].U, mg[0].U)

figure(1)
surf(coarse.x[0:P], coarse.y[0:Q], mg[0].U[0:P,0:Q]'; 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x")
ylabel(L"y")
title("Grid function on coarse grid")
savefig("V.png")
run(`convert V.png -trim coarse_grid_func.png`)

figure(2)
surf(fine.x[0:2P], fine.y[0:2Q], mg[1].U[0:2P,0:2Q]'; 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x")
ylabel(L"y")
title("Prolongation to fine grid")
savefig("PV.png")
run(`convert PV.png -trim prolongated_coarse_grid_func.png`)
