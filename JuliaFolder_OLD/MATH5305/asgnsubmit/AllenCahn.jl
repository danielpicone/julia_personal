module AllenCahn

using OffsetArrays
using ArgCheck

export grid1D, Crank_Nicolson_matrix, semi_implicit
export semi_implicit_test
export grid2D, RectangularDomain

typealias OMatrix{T} OffsetArray{T,2}
typealias OVector{T} OffsetArray{T,1}

immutable Grid1D
    Ω  :: Tuple{Float64,Float64}
    T  :: Float64
    P  :: Int64
    N  :: Int64
    Δx :: Float64
    Δt :: Float64
    x  :: Vector{Float64}  # holds x_1, x_2, ..., x_P
    t  :: OVector{Float64} # holds t_0, t_1, ..., t_N
end # Grid1D

function grid1D(Ω::Tuple{Float64,Float64}, T::Float64, P::Int64, N::Int64)
        a, b = Ω
        Δx = (b-a) / P
        Δt = T / N
        x = linspace(a+Δx/2, b-Δx/2, P)
        t_ = linspace(0, T, N+1)
	      t = OffsetArray(t_, 0:N)
        return Grid1D(Ω, T, P, N, Δx, Δt, x, t)
end

f(u) = u * (1-u^2)

function Crank_Nicolson_matrix(P, α)
    maindiag = fill(1+2α, P)
    maindiag[1] = 1+α
    maindiag[P] = 1+α
    offdiag  = fill(-α, P-1)
    return SymTridiagonal(maindiag, offdiag)
end

"""
U = semi_implicit(ϵ, gr, u0)
Semi-implicit finite difference scheme for the Allen-Cahn equation.
"""
function semi_implicit(ϵ::Float64, gr::Grid1D, u0::Vector{Float64})
    P, N, Δx, Δt = gr.P, gr.N, gr.Δx, gr.Δt
    @argcheck Δx < 10ϵ
    α = Δt/(2*Δx^2)
    σ = Δt/(ϵ^2)
    f(x) = x*(1-x^2)

    # A U^(n+1) = B U^n + σ(3/2 * f(U^n) - 1/2 * f(U^(n-1)))
    # Is what we are trying to solve
    # Create the matrix A

    A = Crank_Nicolson_matrix(P,α)

    # Create the matrix B
    B = SymTridiagonal((1-2*α)*ones(P),α*ones(P-1))
    B[1,1] = 1-α
    B[P,P] = 1-α

    # Initialise U
    U_ = zeros(P+1,N+1)
    U = OffsetArray(U_,0:P,0:N)
    U[1:P,0] = u0
    U[1:P,1] = A \ (B*U[1:P,0] + σ*f.(U[1:P,0]))
    for n=1:N-1
      U[1:P,n+1] = A \ (B*U[1:P,n] + σ*((3/2)*f.(U[1:P,n]) -(1/2)*f.(U[1:P,n-1]) ))
    end
    U[0,0:N] = U[1,0:N]

    return U
end

typealias RectangularDomain Tuple{Float64,Float64,Float64,Float64}
RectangularDomain(a,b,c,d) = RectangularDomain((a,b,c,d))

immutable Grid2D
    Ω :: RectangularDomain
    T :: Float64
    P :: Int64
    Q :: Int64
    N :: Int64
    Δx :: Float64
    Δy :: Float64
    Δt :: Float64
    #x :: LinSpace{Float64}
    x :: Vector{Float64}
    #y :: LinSpace{Float64}
    y :: Vector{Float64}
    t :: OVector{Float64}
end # Grid2D

function grid2D(Ω::RectangularDomain, T::Float64,
                    P::Int64, Q::Int64, N::Int64)
    a, b, c, d = Ω
    Δx = (b-a) / P
    Δy = (d-c) / Q
    Δt = T / N
    x = linspace(a+Δx/2, b-Δx/2, P)
    y = linspace(c+Δy/2, d-Δy/2, Q)
    t_ = linspace(0, T, N+1)
    t = OffsetArray(t_, 0:N)
    return Grid2D(Ω, T, P, Q, N, Δx, Δy, Δt, x, y, t)
end

function Crank_Nicolson_matrix(P::Int64, Q::Int64,
                               α::Float64, β::Float64)
    # Create the offdiagonals
    offdiagx = -1.0*ones(P-1)
    offdiagy = -1.0*ones(Q-1)

    # Create matrix Ax
    Ax = spdiagm((offdiagx,2*ones(P),offdiagx),(-1,0,1))
    Ax[1,1] = 1.0
    Ax[P,P] = 1.0
    Ax = α*Ax

    # Create matrix Ay
    Ay = spdiagm((offdiagy,2*ones(Q),offdiagy),(-1,0,1))
    Ay[1,1] = 1.0
    Ay[P,P] = 1.0
    Ay = β*Ay

    # Create the matrix A using kronecker products and sparse matricies
    Ix = speye(P)
    Iy = speye(Q)
    A = kron(Iy,Ax) + kron(Ay,Ix)

    return ( speye(P*Q) + A )
end

function semi_implicit(ϵ::Float64, gr::Grid2D, u0::Matrix{Float64})
    P, Q, N, Δx, Δy, Δt = gr.P, gr.Q, gr.N, gr.Δx, gr.Δy, gr.Δt
    @assert size(u0) == (P,Q)
    α = Δt/(2*(Δx^2))
    β = Δt/(2*(Δy^2))
    U_ = zeros(P*Q,N)
    U = OffsetArray(U_,1:P*Q,0:N-1)
    B = Crank_Nicolson_matrix(P,Q,α,β)
    #Binv = inv(full(B))
    C = Crank_Nicolson_matrix(P,Q,-α,-β)
    #BinvC = Binv*C
    Bfact = factorize(B)
    U[:,0] = vec(u0)
    final = zeros(P,Q,N)

    for n=0:N-2
      #U[1:P*Q,n+1] = B \ (C*U[1:P*Q,n]) # naive approach
      #U[1:P*Q,n+1] = BinvC*U[1:P*Q,n]
      U[1:P*Q,n+1] = Bfact \ (C*U[1:P*Q,n])
    end
    for n=1:N
      final[:,:,n] = reshape(U[1:P*Q,n-1],P,Q)
    end

    return final
end

end # Module
