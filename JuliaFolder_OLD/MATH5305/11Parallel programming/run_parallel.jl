using MaxEig
maxeig(5, 0.5; output=false) # trigger compilation
@time pmap(maxeig, [10000, 10000], [0.45, 0.65])
