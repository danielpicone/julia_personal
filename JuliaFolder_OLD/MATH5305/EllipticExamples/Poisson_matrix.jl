using Elliptic
using PyPlot

const P = 8
const Q = 6
const Dx = 1.0
const Dy = 1.0

gr = rectgrid(1.0, 1.0, 8, 6)
A = sparse_Poisson_matrix(gr)

spy(A)
title("Matrix for the $(P)x$(Q) discrete Poisson problem")
savefig("Poisson_matrix.png")
