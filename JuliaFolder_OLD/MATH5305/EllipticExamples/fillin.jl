using Elliptic
using PyPlot

const P = 8
const Q = 6
const M = (P-1)*(Q-1)

gr = rectgrid(1.0, 1.0, P , Q)
A = sparse_Poisson_matrix(gr)

AB = band_Poisson_matrix(gr)
band_Cholesky!(AB)
R = band_to_upper(AB)
println("Error using band Cholesky routine:")
err = norm(A-R'*R,1)
println("\t|| A - RᵀR ||₁ = $err")
m = 0
for j=1:M, i=1:M
    if R[i,j] != 0.0
        m += 1
    end
end
println("\tnnz(R) = $m")

figure(1)
subplots_adjust(hspace=0.0)
subplot(1,2,1)
text(10.0, 40.0, L"Poisson matrix $A$")
spy(A)
subplot(1,2,2)
spy(R)
text(8.0, 40.0, L"Cholesky factor $R$")
savefig("fillin_with_border.png")
run(`convert fillin_with_border.png -trim fillin.png`)

F = cholfact(A)
L = sparse(F[:L])
p = F[:p]
RR = L[invperm(p),:]'
println("Error using Suitesparse")
err = norm(A-RR'*RR,1)
println("\t|| A - RᵀR ||₁ = $err")
println("\tnnz(R) = $(nnz(RR))")

figure(2)
subplot(1,2,1)
text(5.0, 40.0, L"Lower triangular factor $L$")
spy(L)
subplot(1,2,2)
text(5.0, 40.0, L"Cholesky factor $R=L^TP$")
spy(R)
savefig("reorder_fillin_with_border.png")
run(`convert reorder_fillin_with_border.png -trim reorder_fillin.png`)
