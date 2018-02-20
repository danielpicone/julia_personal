# Prints out the first N rows of Pascal's triangle.
# Assumes that N is already defined in Julia's workspace.

open("Pascals_triangle.txt","w") do outfile
		binom = zeros(Int64, N, N)
		binom[1,1] = 1
		binom[N,1] = 1
		@printf(outfile,"%2d\n",binom[1,1])
		for n = 2:N
			binom[n,1] = 1
			@printf(outfile,"\n%2d ",binom[n,1])
    	for k = 2:n
    		binom[n,k] = binom[n-1,k-1] + binom[n-1,k]
				@printf(outfile,"%2d ",binom[n,k])
    	end
			@printf(outfile,"\n")
		end
		display(binom)
end
