module Elliptic

using ArgCheck
using OffsetArrays

export rectgrid, d2matrix, sparse_Poisson_matrix, band_Poisson_matrix
export band_Cholesky!, band_to_upper, band_Cholesky_solve!
export Poisson_RHS, insert_boundary_values!, residual!, restrict!, correct!
export Richardson!, Richardson_parameter, rbGaussSeidel!, L2norm
export setup_vcycle, vcycle!, setup_fullmg, fullmg!, multigrid
export OMatrix, OVector, ser_Richardson!, par_Richardson!, C_Richardson!
export Richardson_ver2!

import Base.getindex

typealias OMatrix{T} OffsetArray{T, 2}
typealias OVector{T} OffsetArray{T, 1}

immutable RectGrid
    P  :: Int64
    Q  :: Int64
    x  :: OVector{Float64}
    y  :: OVector{Float64}
    Δx :: Float64
    Δy :: Float64
end

immutable GridEqn
    U  :: OMatrix{Float64}
    b  :: Matrix{Float64}
    r  :: Matrix{Float64}
    gr :: RectGrid
end

immutable MG
    L     :: Int64
    greq  :: OVector{GridEqn}
    A0    :: Matrix{Float64}
end

function rectgrid(Lx, Ly, P, Q)
    x_ = linspace(0, Lx, P+1)
    y_ = linspace(0, Ly, Q+1)
    x = OffsetArray(x_, 0:P)
    y = OffsetArray(y_, 0:Q)
    Δx = Lx / P
    Δy = Ly / Q
    return RectGrid(P, Q, x, y, Δx, Δy)
end

function grideqn(Lx, Ly, P, Q)
    U = OffsetArray(Float64, 0:P, 0:Q)
    fill!(U, 0.0)
    b = zeros(P-1, Q-1)
    r = zeros(P-1, Q-1)
    gr = rectgrid(Lx, Ly, P, Q)
    return GridEqn(U, b, r, gr)
end

getindex(mg::MG, l::Int64) = mg.greq[l]

function multigrid(Lx, Ly, P0, Q0, L)
    greq = OffsetArray(GridEqn, 0:L)
    greq[0] = grideqn(Lx, Ly, P0, Q0)
    P, Q = P0, Q0
    for l = 1:L
        P *= 2
        Q *= 2
        greq[l] = grideqn(Lx, Ly, P, Q)
    end
    A0 = band_Poisson_matrix(greq[0].gr)
    band_Cholesky!(A0)
    M0 = (P0-1)*(Q0-1)
    return MG(L, greq, A0)
end

function setup_vcycle(Lx, Ly, P0, Q0, L, f, g)
    mg = multigrid(Lx, Ly, P0, Q0, L)
    fine = mg[L].gr
    insert_boundary_values!(mg[L].U, g, fine)
    x, y, P, Q = fine.x, fine.y, fine.P, fine.Q
    b = mg[L].b
    for q = 1:Q-1, p = 1:P-1
       b[p,q] = f(x[p], y[q])
    end
    return mg
end

function setup_fullmg(Lx, Ly, P0, Q0, Levels, f, g)
    mg = multigrid(Lx, Ly, P0, Q0, Levels)
    Poisson_RHS!(vec(mg.[0].b), mg.gr0, f, g)
    U0, x, y = mg.U0, mg.gr0.x, mg.gr0.y
    insert_boundary_values!(U0, g, x, y)
    level = mg.greq
    for l = 1:Levels
        gr = level[l].gr
        x, y, P, Q = gr.x, gr.y, gr.P, gr.Q
        insert_boundary_values!(level[l].U, g, x, y)
	for q = 1:Q-1, p = 1:P-1
            level[l].b[p,q] = f(x[p], y[q])
	end
    end
    return mg
end

"""
    d2matrix(P, Δx)

Second derivative matrix (1D)
"""
function d2matrix(P, Δx)
    maindiag = fill(2/Δx^2, P-1)
    offdiag  = fill(-1/Δx^2, P-2)
    return spdiagm((offdiag, maindiag, offdiag), (-1,0,1))
end

function sparse_Poisson_matrix(gr::RectGrid)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    Ax = d2matrix(P, Δx)
    Ix = speye(P-1)
    Ay = d2matrix(Q, Δy)
    Iy = speye(Q-1)
    return kron(Iy, Ax) + kron(Ay, Ix)
end


function band_Poisson_matrix(gr::RectGrid)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    M = (P-1)*(Q-1)
    AB = zeros(P,M)
    a = -1 / Δy^2
    b = -1 / Δx^2
    c = 2/Δx^2 + 2/Δy^2
    for j = 1:M
        AB[1,j]   = a
        AB[P-1,j] = b
        AB[P,j]   = c
    end
    AB[1,1:P-1] = 0.0
    AB[P-1,1:P-1:(P-1)*(Q-2)+1] = 0.0
    return AB
end

function band_Cholesky!(AB::Matrix{Float64})
    s = size(AB, 1) - 1
    M = size(AB, 2)
    for j = 1:M
        for i = max(1,j-s):j-1
            Σ = 0.0
            @simd for k = max(1,j-s):i-1
                @inbounds Σ += AB[s+1+k-i,i] * AB[s+1+k-j,j]
            end
            AB[s+1+i-j,j] = ( AB[s+1+i-j,j] - Σ ) / AB[s+1,i]
        end
        Σ = 0.0
        for k = max(1,j-s):j-1
            Σ += AB[s+1+k-j,j]^2
        end
        AB[s+1,j] = sqrt( AB[s+1,j] - Σ )
    end
    return nothing
end

function band_to_upper(AB::Matrix{Float64})
    s = size(AB,1) - 1
    M = size(AB,2)
    A = zeros(M,M)
    for j = 1:M
        for i = max(1,j-s):j
            A[i,j] = AB[s+1+i-j,j]
        end
    end
    return A
end

"""
    band_Cholesky_solve!(Rb, rhs)

Solves the linear system with Cholesky factor `Rb` in Lapack
band storage mode, as computed by `band_Cholesky!`.  The right-hand 
side vector `rhs` is overwritten by the solution vector.
"""
function band_Cholesky_solve!(Rb::Matrix{Float64}, rhs::Vector{Float64})
    s = size(Rb, 1) - 1
    M = size(Rb, 2)
    for i = 1:M
        Σ = rhs[i]
        for j = max(1,i-s):i-1
            Σ -= Rb[s+1+j-i,i] * rhs[j]
        end
        rhs[i] = Σ / Rb[s+1,i]
    end
    for i = M:-1:1
        Σ = rhs[i]
        for j = i+1:min(M,i+s)
            Σ -= Rb[s+1+i-j,j] * rhs[j]
        end
        rhs[i] = Σ / Rb[s+1,i]
    end
    return nothing
end

function Poisson_RHS!(rhs::Matrix{Float64}, gr::RectGrid, 
                      f::Function, g::Function)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    x, y = gr.x, gr.y
    @argcheck size(rhs) == (P-1, Q-1)
    for q = 1:Q-1, p = 1:P-1
        rhs[p,q] = f(x[p], y[q])
    end
    rΔx2 = 1 / Δx^2
    rΔy2 = 1 / Δy^2
    for p = 1:P-1
        rhs[p,1]   += rΔy2 * g(x[p], y[0])
    end
    for q = 1:Q-1
        rhs[1,q]   += rΔx2 * g(x[0], y[q])
	rhs[P-1,q] += rΔx2 * g(x[P], y[q])
    end
    for p = 1:P-1
	rhs[p,Q-1] += rΔy2 * g(x[p], y[Q])
    end
    return nothing
end

function Poisson_RHS!(rhs::Vector{Float64}, gr::RectGrid, 
                      f::Function, g::Function)
    P, Q = gr.P, gr.Q
    @argcheck length(rhs) == (P-1)*(Q-1)
    rhs_mat = reshape(rhs, P-1, Q-1)
    Poisson_RHS!(rhs_mat, gr, f, g)
    return nothing
end

function Poisson_RHS(gr::RectGrid, f::Function, g::Function)
    M = (gr.P-1)*(gr.Q-1)
    rhs = zeros(M)
    Poisson_RHS!(rhs, gr, f, g)
    return rhs
end

function insert_boundary_values!(U::OMatrix{Float64},
                                 g::Function, gr::RectGrid)
    P, Q, x, y = gr.P, gr.Q, gr.x, gr.y
    iU = indices(U)
    @argcheck length(iU) == 2
    @argcheck indices(U) == (0:P, 0:Q)
    for p = 0:P
        U[p,0] = g(x[p],y[0])
    end
    for q = 1:Q-1
        U[0,q] = g(x[0], y[q])
        U[P,q] = g(x[P], y[q])
    end
    for p = 0:P
        U[p,Q] = g(x[p], y[Q])
    end
    return nothing
end

"""
    residual!(r, f, U, gr)
    
Computes the residual `r = f - AU`.
"""
function residual!(r::Matrix{Float64}, f::Matrix{Float64}, 
                            U::OMatrix{Float64}, gr::RectGrid)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(f) == (P-1, Q-1)
    @argcheck indices(U) == (0:P, 0:Q)
    for q = 1:Q-1 
        @simd for p = 1:P-1
            r[p,q] = f[p,q] - ( 
                          ( -U[p+1,q] + 2U[p,q] - U[p-1,q] ) / Δx^2 
	                + ( -U[p,q+1] + 2U[p,q] - U[p,q-1] ) / Δy^2 )
        end
    end
    return nothing
end

"""
    Richardson!(U, r, ν, ω, gr, f)

Perform `ν` Richardson iterations with relaxation parameter `ω`
on the grid `gr`.
Note: `U`  is (0:P)x(0:Q) but `r` and `f` are (P-1)x(Q-1).
"""
function Richardson!(U::OMatrix{Float64}, r::Matrix{Float64},
                     ν::Int64, ω::Float64, 
                     gr::RectGrid, b::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(b) == (P-1, Q-1)
    for k = 1:ν
	for q = 1:Q-1 
	    for p = 1:P-1
                r[p,q] = b[p,q] - ( 
                          ( -U[p+1,q] + 2U[p,q] - U[p-1,q] ) / Δx^2 
	                + ( -U[p,q+1] + 2U[p,q] - U[p,q-1] ) / Δy^2 )
                U[p,q] += ω * r[p,q]
            end
        end
    end
    return nothing
end

function Richardson_ver2!(U::OMatrix{Float64}, r::Matrix{Float64},
                     ν::Int64, ω::Float64, 
                     gr::RectGrid, b::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(b) == (P-1, Q-1)
    U_ = U.parent
    for k = 1:ν
	for q = 1:Q-1 
	    @simd for p = 1:P-1
            @inbounds begin
            r[p,q] = b[p,q] - ( 
                    ( -U_[p+2,q+1] + 2U_[p+1,q+1] - U_[p,q+1] ) / Δx^2 
	              + ( -U_[p+1,q+2] + 2U_[p+1,q+1] - U_[p+1,q] ) / Δy^2 )
            U_[p+1,q+1] += ω * r[p,q]
            end
        end
	end
    end
    return nothing
end

function ser_Richardson!(U::OMatrix{Float64}, r::Matrix{Float64},
                     ν::Int64, ω::Float64, gr::RectGrid, b::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(b) == (P-1, Q-1)
    ccall((:ser_Richardson, "./libelliptic.so"), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Cdouble, Cdouble, 
	   Cint, Cint, Ptr{Cdouble}),
	   U.parent, r, ν, ω, Δx, Δy, P, Q, b)
    return nothing
end

function par_Richardson!(U::OMatrix{Float64}, r::Matrix{Float64},
                     ν::Int64, ω::Float64, gr::RectGrid, b::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(b) == (P-1, Q-1)
    ccall((:par_Richardson, "./libelliptic.so"), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Cdouble, Cdouble, 
	   Cint, Cint, Ptr{Cdouble}),
	   U.parent, r, ν, ω, Δx, Δy, P, Q, b)
    return nothing
end

function C_Richardson!(U::OMatrix{Float64}, r::Matrix{Float64},
                       ν::Int64, ω::Float64, gr::RectGrid, b::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(r) == (P-1, Q-1)
    @argcheck size(b) == (P-1, Q-1)
    ccall((:C_Richardson, "./libelliptic.so"), Void,
         (Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, 
          Cdouble, Cdouble, Cdouble, Ptr{Cdouble}), 
         P, Q, U.parent, r, ν, ω, Δx, Δy, b)
    return nothing
end

function Richardson_parameter(gr::RectGrid, case::Symbol)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    τ = (Δy/Δx) + (Δx/Δy)
    if case == :smoother
        τPQ = (Δy/Δx) * cos(π/P) + (Δx/Δy) * cos(π/Q)
        ω = Δx * Δy / (2τ + τPQ)
    elseif case == :solver
        ω = Δx * Δy / (2τ)
    else
        ArgumentError("Unrecognised case")
    end
    return ω
end

function rbGaussSeidel!(U::OMatrix{Float64}, ν::Int64, 
                        gr::RectGrid, f::Matrix{Float64})
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    @argcheck indices(U) == (0:P, 0:Q)
    @argcheck size(f) == (P-1, Q-1)
    ΔxΔy = Δx * Δy
    τx = Δy / Δx
    τy = Δx / Δy
    τ = τx + τy
    for k = 1:ν
        for q = 1:Q-1
	    r = 2 - mod(q, 2)
	    for p = r:2:P-r
	        U[p,q] = ( f[p,q] * ΔxΔy 
                         + τx * ( U[p+1,q] + U[p-1,q] )
			 + τy * ( U[p,q+1] + U[p,q-1] )
			 ) / ( 2τ )
	    end
	    r = 2 - mod(q-1, 2)
	    for p = r:2:P-r
	        U[p,q] = ( f[p,q] * ΔxΔy 
		         + τx * ( U[p+1,q] + U[p-1,q] )
			 + τy * ( U[p,q+1] + U[p,q-1] )
			 ) / ( 2τ )
	    end
	end
    end
    return nothing
end

function L2norm(V::OMatrix{Float64}, gr::RectGrid)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    Σ = 0.0
    for q = 0:Q, p = 0:P
        Σ += V[p,q]^2
    end
    return sqrt(Σ*Δx*Δy)
end

function L2norm(V::Matrix{Float64}, gr::RectGrid)
    P, Q, Δx, Δy = gr.P, gr.Q, gr.Δx, gr.Δy
    Σ = 0.0
    for q = 1:Q-1, p = 1:P-1
        Σ += V[p,q]^2
    end
    return sqrt(Σ*Δx*Δy)
end

"""
    restrict!(RV, V)

Restriction operator from a `(2P)x(2Q)` grid to a `PxQ` grid.
Either `RV` is `(P-1)x(Q-1)` and V` is `(2P-1)x(2Q-1)`, or else
`RV` is `(0:P)x(0:Q)` and `V` is `(0:2P)x(0:2Q)`.  
Even in the latter case, however, only the values of `RV` at 
the interior nodes are computed, so the calling program 
should assign the desired boundary values.
"""
function restrict!(RV::Matrix{Float64}, V::Matrix{Float64})
    P = size(RV, 1) + 1
    Q = size(RV, 2) + 1
    @argcheck size(V) == (2P-1, 2Q-1)
    for q = 1:Q-1, p = 1:P-1
        RV[p,q] = (  V[2p-1,2q+1] + 2V[2p,2q+1] +  V[2p+1,2q-1]
                  + 2V[2p-1,2q]   + 4V[2p,2q]   + 2V[2p+1,2q-1]
                  +  V[2p-1,2q-1] + 2V[2p,2q-1] +  V[2p+1,2q-1] 
                  ) / 16
    end
    return nothing
end

function restrict!(RV::OMatrix{Float64}, V::OMatrix{Float64})
    i = indices(RV)
    P, Q = i[1].stop, i[2].stop
    @argcheck i[1].start == 0
    @argcheck i[2].start == 0
    @argcheck indices(V) == (0:2P, 0:2Q)
    for q = 1:Q-1, p = 1:P-1
        RV[p,q] = (  V[2p-1,2q+1] + 2V[2p,2q+1] +  V[2p+1,2q-1]
                  + 2V[2p-1,2q]   + 4V[2p,2q]   + 2V[2p+1,2q-1]
                  +  V[2p-1,2q-1] + 2V[2p,2q-1] +  V[2p+1,2q-1] 
                  ) / 16
    end
    return nothing
end

"""
    correct!(U, V)

Overwrites `U` with `U+PV`.  Here, `U` is `(0:2P)x(0:2Q)`
and `V` is `(0:P)x(0:Q)`, corresponding to
functions on a `(2P)x(2Q)` grid and a `PxQ` grid, respectively.
The boundary values of `U` are not affected by the prolongation `PV`.
"""
function correct!(U::OMatrix{Float64}, V::OMatrix{Float64})
    i = indices(V)
    @argcheck length(i) == 2
    P, Q = i[1].stop, i[2].stop
    @argcheck indices(U) == (0:2P, 0:2Q)
    U[1,1] += ( V[0,0] + V[1,0] + V[0,1] + V[1,1] ) / 4
    for p = 1:P-1
        U[2p,1]   += ( V[p,0] + V[p,1] ) / 2
        U[2p+1,1] += ( V[p,0] + V[p+1,0] + V[p,1] + V[p+1,1] ) / 4
    end
    for q = 1:Q-1
   	U[1,2q] += ( V[0,q] + V[1,q] ) / 2
        for p = 1:P-1
            U[2p,2q]   += V[p,q]
   	    U[2p+1,2q] += ( V[p,q] + V[p+1,q] ) / 2
        end
        U[1,2q+1] += ( V[0,q] + V[1,q] + V[0,q+1] + V[1,q+1] ) / 4
        for p = 1:P-1
	    U[2p,2q+1]   += ( V[p,q] + V[p,q+1] ) / 2
	    U[2p+1,2q+1] += ( V[p,q] + V[p+1,q] + V[p,q+1] + V[p+1,q+1] ) / 4
	end
    end
    return nothing
end

function vcycle!(mg::MG, smooth!::Function; 
                 presmooth::Int64=3, postsmooth::Int64=3)
    vcycle_!(mg, mg.L, smooth!; pre=presmooth, post=postsmooth)
    return nothing
end

function vcycle_!(mg::MG, L::Int64, smooth!::Function; 
                  pre::Int64=3, post::Int64=3)
    @argcheck L <= mg.L
    for l = L:-1:1
        smooth!(l, pre, mg)
        residual!(mg[l].r, mg[l].b, mg[l].U, mg[l].gr)
        restrict!(mg[l-1].b, mg[l].r)
	fill!(mg[l-1].U, 0.0)
    end
    band_Cholesky_solve!(mg.A0, vec(mg[0].b))
    P0, Q0 = mg[0].gr.P, mg[0].gr.Q
    U0, b0 = mg[0].U, mg[0].b
    for q = 1:Q0-1, p = 1:P0-1
        U0[p,q] = b0[p,q]
    end
    for l = 1:L 
        correct!(mg[l].U, mg[l-1].U)
        smooth!(l, post, mg)
    end
    return nothing
end

function fullmg!(mg::MG, smooth!::Function; 
                 presmooth::Int64=3, postsmooth::Int64=3)
    band_Cholesky_solve!(mg.A0, vec(mg[0].b))
    U0[2:P0,2:Q0] .= matb0
    P0, Q0 = mg[0].gr.P, mg[0].gr.Q
    U0, b0 = mg[0].U, mg[0].b
    for q = 1:Q0-1, p = 1:P0-1
        U0[p,q] = b0[p,q]
    end
    correct!(mg[1].U, mg[0].U)
    for l = 1:L-1
        vcycle_!(mg, l; pre=presmooth, post=postsmooth)
        correct!(mg[l+1].U, mg[l].U)
        smooth!(l+1, postsmooth, mg)
    end
    return nothing
end

end # module
