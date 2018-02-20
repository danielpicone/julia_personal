module MaxEig

export maxeig

const MAXITS = 10000

function maxeig(N::Int64, α::Float64; output=true)

    A = Array(Float64, N, N)
    v = Array(Float64, N)
    w = Array(Float64, N)
    h = 1 / N
    c = h^α / gamma(α)
    if output
        println("Started: $N-by-$N matrix with α = $α")
    end
    start = time()
    for j=1:N, i=1:N
        @inbounds A[i,j] = c * ( abs(i-j) + 1 )^(α-1) 
    end
    λ, iter = power_method!(A, v, w)
    matvec!(w, A, v) # w =  Av
    w[:] .= w - λ*v
    resid = norm(w, 2)
    finish = time()
    if output
        println("Finished: $N-by-$N matrix with α = $α")
        @printf("\tMax eigenvalue = %0.14f\n", λ)
        @printf("\tUsed %0d iterations\n", iter)
        @printf("\t|| Av - λv ||₂ = %0.2e\n", resid)
        @printf("\tElapsed time = %0.4f seconds\n", finish - start)
    end 
end

function power_method!(A, v, w)
    rtol = 10 * eps(Float64)
    rand!(v)
    matvec!(w, A, v)
    λ = dot(w,v) / dot(v,v)
    v[:] .= w
    iter = 0
    prev_λ = 0.0
    for k = 1:MAXITS
        iter += 1
        matvec!(w, A, v)
        λ = dot(w,v) / dot(v,v)
        if abs(λ-prev_λ) < rtol * λ
            break
        end
        v[:] .= w / norm(w, 2)
        prev_λ = λ
    end
    return λ, iter
end

"""
  matvec!(y, A, x)

Matrix-vector product y = Ax
"""
function matvec!{T<:AbstractFloat}(y::Vector{T}, A::Matrix{T}, x::Vector{T})
    nrows = size(A, 1)
    ncols = size(A, 2)
    @assert size(x, 1) == ncols
    @assert size(y, 1) == nrows
    fill!(y, 0.0)
    for j = 1:ncols
        @simd for i = 1:nrows
            @inbounds y[i] += A[i,j] * x[j]
        end
    end
end

end # module
