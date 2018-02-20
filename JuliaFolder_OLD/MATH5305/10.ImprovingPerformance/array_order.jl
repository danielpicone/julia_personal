using ArgCheck

const M = 15_000
const N = 15_000

A = Float64[ 1.0 / ( 1.0 + abs(M-N) ) for m = 1:M, n = 1:N ]
x = ones(M)
y = rand(N)

iprod_builtin{T<:Real}(x::Vector{T}, A::Matrix{T}, y::Vector{T}) = dot(x, A*y)

function iprod_cols{T<:Real}(x::Vector{T}, A::Matrix{T}, y::Vector{T})
    M, N = size(A)
    @argcheck length(x) == M
    @argcheck length(y) == N
    double_sum = zero(T)
    for n = 1:N
        Σ = zero(T)
        for m = 1:M
	        Σ += x[m] * A[m,n]
	      end
	  double_sum += Σ * y[n]
    end
    return double_sum
end

function iprod_cols_opt{T<:Real}(x::Vector{T}, A::Matrix{T}, y::Vector{T})
    M, N = size(A)
    @argcheck length(x) == M
    @argcheck length(y) == N
    double_sum = zero(T)
    for n = 1:N
        Σ = zero(T)
        @simd for m = 1:M
	         @inbounds Σ += x[m] * A[m,n]
	         end
	  double_sum += Σ * y[n]
    end
    return double_sum
end

function iprod_rows{T<:Real}(x::Vector{T}, A::Matrix{T}, y::Vector{T})
    M, N = size(A)
    @argcheck length(x) == M
    @argcheck length(y) == N
    double_sum = zero(T)
    for m = 1:M
        Σ = zero(T)
        for n = 1:N
	    Σ += A[m,n] * y[n]
	end
	double_sum += x[m] * Σ
    end
    return double_sum
end

function iprod_rows_opt{T<:Real}(x::Vector{T}, A::Matrix{T}, y::Vector{T})
    M, N = size(A)
    @argcheck length(x) == M
    @argcheck length(y) == N
    double_sum = zero(T)
    for m = 1:M
        Σ = zero(T)
        @simd for n = 1:N
	    @inbounds Σ += A[m,n] * y[n]
	end
	double_sum += x[m] * Σ
    end
    return double_sum
end

A1 = rand(2,2)
x1 = rand(2)
y1 = rand(2)

for func in (iprod_rows, iprod_cols, iprod_rows_opt, iprod_cols_opt)
    func(x1, A1, y1) # compile
    secs = @elapsed for trial = 1:20
        func(x, A, y)
    end
    @printf("%15s  %8.4f\n", string(func), secs)
end
