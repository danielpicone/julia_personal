using ArgCheck

function first_expm(N::Int64, A::Matrix{Float64})
    @argcheck size(A,1) == size(A,2)
    C = one(A)
    S = one(A)
    for k = 1:N-1
        C *= (1/k) * A # = A^k / k!
        S += C
    end
    return S
end

function add!(S::Matrix{Float64}, Y::Matrix{Float64})
    
    m, n = size(S)
    @argcheck size(Y) == (m,n)
    for j = 1:n
        for i = 1:m
            @inbounds S[i,j] += Y[i,j]
        end
    end
    return nothing
end

function second_expm(N::Int64, A::Matrix{Float64})
    @argcheck size(A,1) == size(A,2)
    C = one(A)
    S = one(A)
    Y = zeros(A)
    for k = 1:N-1
        A_mul_B!(Y, C, A) # Y = C * A
        scale!(Y, 1/k)    # Y = C * A / k = A^k / k!
        add!(S, Y)        # S = S + Y
    end
    return S
end

function test_matrix(N)
    A = Float64[ 1 / ( 1 + abs(i-j) ) for i=1:N, j=1:N ]
    scale!(A, 1/(2*norm(A)))
    return A
end
