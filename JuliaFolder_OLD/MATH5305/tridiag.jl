using ArgCheck

"""
    ℓ, d, e = trifact(β, α, ρ)
    L, U = trifact(A)

Computes the tridiagonal LU factorisation (without pivoting), if it exists.

    | α₁ ρ₁          |   | 1             | | d₁ e₁          |
    | β₁ α₂ ρ₂       |   | ℓ₁ 1          | |    d₂ e₂       |
    |    β₂ α₃ ρ₃    | = |    ℓ₂ 1       | |       d₃ e₃    |.
    |       β₃ α₄ ρ₄ |   |       ℓ₃ 1    | |          d₄ e₄ |
    |          β₄ α₅ |   |          ℓ₄ 1 | |             d₅ |

The low-level method passes and returns arrays.  For the high-level method,
`A` is of type `Tridiagonal`, whereas `L` and `U` are of type `Bidiagonal`.
"""
function trifact{T<:Number}(β::Vector{T}, α::Vector{T}, ρ::Vector{T})
    n = length(α)
    @argcheck length(β) == n-1
    @argcheck length(ρ) == n-1
    d = zeros(α)
    e = zeros(ρ)
    ℓ = zeros(β)
    d[1] = α[1]
    for j = 1:n-1
        ℓ[j] = β[j] / d[j]
        e[j] = ρ[j]
        d[j+1] = α[j+1] - ℓ[j] * e[j]
    end
    return ℓ, d, e
end

"""
    trifact!(β, α, ρ)
    trifact!(A)

In-place version of `trifact`.
"""
function trifact!{T<:Number}(β::Vector{T}, α::Vector{T}, ρ::Vector{T})
    n = length(α)
    @argcheck length(β) == n-1
    @argcheck length(ρ) == n-1
    for j = 1:n-1
        β[j] /= α[j]
        α[j+1] -= β[j] * ρ[j]
    end
    return nothing
end

function trifact{T<:Number}(A::Tridiagonal{T})
    ℓ, d, e, = trifact(A.dl, A.d, A.du)
    L = Bidiagonal(ones(d), ℓ, 'L')
    U = Bidiagonal(d, e, 'U')
    return L, U
end

function trifact!{T<:Number}(A::Tridiagonal{T})
    trifact!(A.dl, A.d, A.du)
    return nothing
end

"""
    x = trisolve(ℓ, d, e, b)
    x = trisolve(L, U, b)

Solve the tridiagonal linear system Ax=b given the vectors `ℓ`, `d`, `e`
from the LU factorisation of A, or the `Bidiagonal` factors `L` and `U` themselves.
"""
function trisolve{T<:Number}(ℓ::Vector{T}, d::Vector{T}, e::Vector{T}, 
                             b::Vector{T})
    n = length(d)
    @argcheck length(ℓ) == n-1
    @argcheck length(e) == n-1
    @argcheck length(b) == n
    y = zeros(b)
    y[1] = b[1]
    for j = 2:n
        y[j] = b[j] - ℓ[j-1] * y[j-1]
    end
    x = zeros(b)
    x[n] = y[n] / d[n]
    for j = n-1:-1:1
        x[j] = ( y[j] - e[j] * x[j+1] ) / d[j]
    end
    return x
end

function trisolve{T<:Number}(L::Bidiagonal{T}, U::Bidiagonal{T}, b::Vector{T})
    x = trisolve(L.ev, U.dv, U.ev, b)
    return x
end

"""
    trisolve!(ℓ, d, e, b)

In-place version of `trisolve`.
"""
function trisolve!{T<:Number}(ℓ::Vector{T}, d::Vector{T}, e::Vector{T}, 
                             b::Vector{T})
    n = length(d)
    @argcheck length(ℓ) == n-1
    @argcheck length(e) == n-1
    @argcheck length(b) == n
    for j = 2:n
        b[j] -= ℓ[j-1] * b[j-1]
    end
    b[n] /= d[n]
    for j = n-1:-1:1
        b[j] = ( b[j] - e[j] * b[j+1] ) / d[j]
    end
    return nothing
end

T = Float64
#T = Rational

β = T[1, 0, -2]
α = T[7, -9, 3, 10]
ρ = T[-2, 1, 1]
b = T[3, -2, 8, 5]

A = Tridiagonal(β, α, ρ)
L, U = trifact(A)
err = maxabs(A - L*U)
println("\nMax residual in LU factorisation = $err")

x = trisolve(L, U, b)
err = maxabs(b - A*x)
println("\nMax residual in solution x = $err")

println("\nA: ")
display(A)
println("\nL: ")
display(L)
println("\nU: ")
display(U)

keep_A = copy(A)
keep_b = copy(b)

trifact!(A)
println("\nIn-place factorisation:")
display(A)

trisolve!(A.dl, A.d, A.du, b)
err = maxabs(keep_b - keep_A * b)
println("\nMax residual for in-place solve = $err")
