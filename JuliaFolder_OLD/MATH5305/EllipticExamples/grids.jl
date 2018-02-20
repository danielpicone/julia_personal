using PyPlot

const P0 = 4
const Q0 = 4

function gridlines(x, y)
    P = length(x) - 1
    Q = length(y) - 1
    xx = zeros(3(P+Q+2))
    yy = zeros(3(P+Q+2))
    zz = zeros(3(P+Q+2))
    for q = 0:Q
        xx[1+3q] = x[1]
	xx[2+3q] = x[P+1]
	xx[3+3q] = NaN
	yy[1+3q] = y[q+1]
	yy[2+3q] = y[q+1]
	yy[3+3q] = NaN
    end
    for p = 0:P
        xx[3(Q+1)+1+3p] = x[p+1]
        xx[3(Q+1)+2+3p] = x[p+1]
        xx[3(Q+1)+3+3p] = NaN
        yy[3(Q+1)+1+3p] = y[1]
        yy[3(Q+1)+2+3p] = y[Q+1]
        yy[3(Q+1)+3+3p] = NaN
    end
    return xx, yy
end

figure(1)
ax = gca(projection="3d")
P, Q, L = P0, Q0, 2
for l = 0:L
    x = linspace(0, 1, P+1)
    y = linspace(0, 1, Q+1)
    xx, yy = gridlines(x, y)
    zz = fill(l/L, length(xx))
    plot3D(xx, yy, zz)
    P *= 2
    Q *= 2
end
ax[:view_init](20,60)
axis("off")
savefig("multigrids_with_border.png")
run(`convert multigrids_with_border.png -trim multigrids.png`)
