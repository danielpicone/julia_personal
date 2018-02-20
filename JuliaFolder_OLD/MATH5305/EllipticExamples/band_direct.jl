using Elliptic
using OffsetArrays

include("problem1.jl")
P0 = 2
Q0 = 2
nrows = 7

function band_solve(P, Q)
    gr = rectgrid(Lx, Ly, P, Q)
    start = time()
    AB = band_Poisson_matrix(gr)
    band_Cholesky!(AB)
    b = Poisson_RHS(gr, f, g)
    band_Cholesky_solve!(AB, b)
    U = OffsetArray(Float64, 0:P, 0:Q)
    U[1:P-1,1:Q-1] = reshape(b, P-1, Q-1)
    insert_boundary_values!(U, g, gr)
    finish = time()
    elapsed = finish - start
    maxerror = maxabs(U - exact_u.(gr.x',gr.y)')
    nonzeros = 0
    for elt in AB[:]
        nonzeros += (elt==0.0) ? 0 : 1
    end
    return elapsed, maxerror, nonzeros
end

err = zeros(nrows)
elapsed = zeros(nrows)
nonzeros = zeros(Int64, nrows)

# Don't want to include compile time in first row.
dummy_elapsed, dummy_err, dummy_nonzeros = band_solve(P0,Q0)

P, Q = P0, Q0
@printf("Poisson problem on PxQ grid; band direct solver.\n\n")
@printf("%5s %5s %8s %10s %6s %6s %6s %10s %6s %7s\n\n",
        "P", "Q", "M", "error", "rate", "secs", "rate",
        "nnz", "rate", "nnz/M^2")
for row = 1:nrows
    P *= 2
    Q *= 2
    M = (P-1)*(Q-1)
    elapsed[row], err[row], nonzeros[row] = band_solve(P, Q)
    fracnz = nonzeros[row]/M^2
    if row == 1
        @printf("%5d %5d %8d %10.2e %6s %6.3f %6s %10.2e %6s %6.2f%%\n",
                P, Q, M, err[row], " ", elapsed[row], " ",
                nonzeros[row], " ", 100*fracnz)
    else
        err_rate = log2(err[row-1]/err[row])
        elapsed_rate = log(elapsed[row]/elapsed[row-1]) / log(4.0)
        nnz_rate = log(nonzeros[row]/nonzeros[row-1]) / log(4.0)
        @printf("%5d %5d %8d %10.2e %6.3f %6.3f %6.3f %10.2e %6.3f %6.2f%%\n",
                P, Q, M, err[row], err_rate, elapsed[row], elapsed_rate,
                nonzeros[row], nnz_rate, 100*fracnz)
    end
end

