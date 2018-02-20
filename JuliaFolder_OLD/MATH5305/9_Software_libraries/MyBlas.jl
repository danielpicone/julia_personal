module MyBlas

export daxpy!, fortran_daxpy!, c_daxpy!, impersonate_c_daxpy!

#const SYSTEM_BLAS = "/usr/lib64/libopenblas.so"
const SYSTEM_BLAS = "/usr/lib/libopenblas.so.0"

function daxpy!{T<:AbstractFloat}(n, a::T, x::Vector{T}, incx, 
                                           y::Vector{T}, incy)
    @assert length(x) <= 1 + (n-1)*incx
    @assert length(y) <= 1 + (n-1)*incy
    ix = 1
    iy = 1
    for i = 1:n
        y[ix] = a * x[ix] + y[iy]
        ix += incx
        iy += incy
    end
    return nothing
end

function fortran_daxpy!(n::Int32, a::Float64, 
                        x::Vector{Float64}, incx::Int32,
                        y::Vector{Float64}, incy::Int32)
    @assert length(x) <= 1 + (n-1)*incx
    @assert length(y) <= 1 + (n-1)*incy
    ccall((:__my_fblas_MOD_daxpy, "./libmyblas.so"), Void, 
          (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, 
           Ref{Cint}, Ref{Cdouble}, Ref{Cint}), 
           n, a, x, incx, y, incy)
    return nothing
end

function c_daxpy!(n::Int32, a::Float64, 
                  x::Vector{Float64}, incx::Int32,
                  y::Vector{Float64}, incy::Int32)
    @assert length(x) <= 1 + (n-1)*incx
    @assert length(y) <= 1 + (n-1)*incy
    ccall((:cblas_daxpy, SYSTEM_BLAS), Void, 
          (Cint, Cdouble, Ptr{Cdouble}, Cint, 
           Ptr{Cdouble}, Cint), 
           n, a, x, incx, y, incy)
    return nothing
end

function impersonate_c_daxpy!(n::Int32, a::Float64, 
                              x::Vector{Float64}, incx::Int32,
                              y::Vector{Float64}, incy::Int32)
    @assert length(x) <= 1 + (n-1)*incx
    @assert length(y) <= 1 + (n-1)*incy
    ccall((:impersonate_c_daxpy, "./libmyblas.so"), 
          Void, (Cint, Cdouble, Ptr{Cdouble}, 
                 Cint, Ptr{Cdouble}, Cint), 
           n, a, x, incx, y, incy)
    return nothing
end

end # module
