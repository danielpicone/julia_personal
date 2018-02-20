module MyStuff

export quadform

function quadform{T<:AbstractFloat}(A::Matrix{T}, v::Vector{T})
    n = size(A, 1)
    @assert size(A, 2) == n
    @assert size(v, 1) == n
    ans = 0.0
    for j = 1:n
        inner_sum = 0.0
        for i = 1:n
            inner_sum += v[i] * A[i,j] 
        end
        ans += inner_sum * v[j]
    end
    return ans
end # function quadform

end # module MyStuff
