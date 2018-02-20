# Prints out the first N rows of Pascal's triangle.
# Assumes that N is already defined in Julia's workspace.

binom = zeros(Int64, N)
for n = 1:(N-1)
	for k = n:-1:1
		binom[1] = 1
    	binom[k+1] = binom[k] + binom[k+1]
    end
end
display(binom)
