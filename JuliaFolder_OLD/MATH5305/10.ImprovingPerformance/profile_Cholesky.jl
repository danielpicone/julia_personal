using Elliptic
using OffsetArrays

const P = 256
const Q = 256
include("problem1.jl")

gr = rectgrid(Lx, Ly, P, Q)
AB_save = band_Poisson_matrix(gr)
b_save = Poisson_RHS(gr, f, g)
AB = zero(AB_save)
b = zero(b_save)

function band_solve!(gr)
    copy!(AB, AB_save)
    copy!(b, b_save)
    band_Cholesky!(AB)
    band_Cholesky_solve!(AB, b)
    return nothing
end
