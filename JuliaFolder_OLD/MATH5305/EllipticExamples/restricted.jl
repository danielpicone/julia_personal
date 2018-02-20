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
for q = 0:2Q, p = 0:2P
    mg[1].U[p,q] = v(fine.x[p], fine.y[q])
end
restrict!(mg[0].U, mg[1].U)
insert_boundary_values!(mg[0].U, v, coarse)

figure(1)
surf(fine.x[0:2P], fine.y[0:2Q], mg[1].U[0:2P,0:2Q]'; 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x")
ylabel(L"y")
title("Grid function on fine grid")
savefig("V.png")
run(`convert V.png -trim fine_grid_func.png`)

figure(2)
surf(coarse.x[0:P], coarse.y[0:Q], mg[0].U[0:P,0:Q]'; 
     rstride=1, cstride=1, cmap="autumn")
xlabel(L"x")
ylabel(L"y")
title("Restriction to coarse grid")
savefig("RV.png")
run(`convert RV.png -trim restricted_fine_grid_func.png`)
