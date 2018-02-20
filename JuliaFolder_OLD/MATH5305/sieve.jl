# This is the seive of Erathosthenes

p = 2*ones(Int32,N)
p[1] = 0
k = 1
i = 2
flag = 5
while (i <= N-1)&&(flag == 5)
	while (p[k] != 2)&&(k<=N)
		k += 1
	end
	if (k == N-1)
		flag = 0
	end
	@printf("%d\n",k)
	p[k] = 1
	for j=(i+2):(N)
		if j % k == 0
			p[j] = 0
		end
	end
	i += 1
end
