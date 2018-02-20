# Jacobi Iterative solver for diagonally dominant linear systems




#A = [3 0 1;
#     1 8 2;
#     2 1 4]

numrow = 10000
numcol = 10000
A = rand(numrow,numcol)+50*eye(numrow)

#b = [2; 8; 3]
b = 10*rand(numrow,1)

x = zeros(numcol)

function itersum(A,x,i,numcol)
  sum = 0
  for j = 1:numcol
    if j==i
    elseif j!=i
      sum += A[i,j]*x[j]
    end
  end
  return sum
end


startjac = time()
iter = 100
for μ=1:iter
  for i=1:numcol
    x[i] = 1/A[i,i] * (b[i] - itersum(A,x,i,numcol))
  end
end
stopjac = time()
println("This took ",stopjac-startjac," seconds")

startsolv = time()
exactx = A\b
stopsolv = time()
println("Solving took ", stopsolv-startsolv," seconds")

startfact = time()
B = factorize(A)
exactxfact = B\b
stopfact = time()
println("Factorizing and solving took ", stopfact-startfact," seconds")


#println("The exact solution is ",exactx)
#println("The error is ",x-exactx)
# Did it converge?
sumerr = sum(abs(x-exactx))
if sumerr < 10
  println("Yes it converged!!")
else
  println("Nah that shit fucked up")
end
