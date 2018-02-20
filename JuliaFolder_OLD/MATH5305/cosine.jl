# Use N terms of the Taylor expansion to approximate cos(x).
# The script assumes that you have already defined x and N in Julia's
# workspace.

a = x
s = a
xsqr = x^2
for k = 1:N-1
    a *= -xsqr / ((2k)*(2k-1)) # (k+1)th term
    s += a                     # sum of first k+1 terms
end

println("Taylor approx using $N terms:    $s")
println("Value given by library function: $(cos(x))")
