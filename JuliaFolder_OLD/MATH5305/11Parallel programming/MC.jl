module MC

export Eκ, stdκ, diffusivity
export gridpoints, linear_system!, u_at_zero!

Eκ(x::Float64, L::Float64) = 1 + x/L
stdκ(s::Float64, C::Float64) = (C/2)*sqrt(zeta(2s)/3)

function κ(x::Float64, ω::Vector{Float64},
	   L::Float64, s::Float64, C::Float64)
    M = length(ω)
    val = 0.0
    for m = 1:M
        val += ω[m] * sin(m*π*x/L) / m^s
    end
    return Eκ(x, L) + C * val
end

function diffusivity(L::Float64, s::Float64, C::Float64)
    return (x,ω) -> κ(x, ω, L, s, C)
end

function κ(x::Vector{Float64}, ω::Vector{Float64},
	   L::Float64, s::Float64, C::Float64)
    P = length(x)
    v = Array(Float64, P)
    for p = 1:length(x)
        v[p] = κ(x[p], ω, L, s, C)
    end
    return v
end
    
function gridpoints(L::Float64, P::Int64)
    Δx = L / ( P + 0.5 )
    x = Float64[(p-1/2)*Δx for p = 1:P]
    return x, Δx
end

function linear_system!{T<:AbstractFloat}(
	                A::SymTridiagonal{T}, b::Vector{T}, 
			κ::Function, γ::T, Δx::T)
    P = size(b, 1)
    @assert size(A) == (P, P)
    for p = 1:P
	x = p * Δx
        A.dv[p] = κ(x)
    end
    b[1] = γ * Δx / κ(0.0)
    for p = 1:P-1
        A.ev[p] = -A.dv[p]
	b[p+1] = 0.0
    end
    for p = P:-1:2
        A.dv[p] += A.dv[p-1]
    end
    return nothing
end

function u_at_zero!{T<:AbstractFloat}(A::SymTridiagonal{T}, b::Vector{T}, 
	                      κ::Function, γ::T, Δx::Float64)
    linear_system!(A, b, κ, γ, Δx)
    F = ldltfact!(A)
    U = F \ b
    u0 = U[1] + γ * Δx / ( 2 * κ(0.0) )
    return u0
end

end # module
