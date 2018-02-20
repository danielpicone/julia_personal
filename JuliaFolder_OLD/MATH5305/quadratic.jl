"""
    r1, r2 = quadroots(a, b, c)

Compute the roots of the quadratic equation ax^2+bx+c=0.
"""
function quadroots{T<:Real}(a::T, b::T, c::T)
    if a == zero(a) # not a quadratic equation
        throw(DomainError())
    end
    D = b^2 - 4*a*c
    if D < zero(D) # Complex conjugate roots
        s = sqrt(-D)
        r1 = complex(-b, +s) / (2a)
        r2 = complex(-b, -s) / (2a)
    else # Real roots
        s = sqrt(D)
        if b >= zero(D)
            r1 = (-b-D) / (2a)
            r2 = c / (a*r1)
        else
            r2 = (-b+D) / (2a)
            r1 = c / (a*r2)
        end
    end
    return r1, r2
end
