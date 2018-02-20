using MyStuff

const M = 1000

A = rand(M, M)
v = rand(M)

println("$M-by-$M quadratic form.")
ans1 = dot(v,A*v)
@printf("Julia instrinsics: ans = %20.14e\n", ans1)
ans2 = quadform(A, v)
@printf("quadform function: ans = %20.14e\n", ans2)
