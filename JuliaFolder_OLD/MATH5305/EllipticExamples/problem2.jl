const Lx = 1.0
const Ly = 1.0

f(x, y) = 32( x*(Lx-x) + y*(Ly-y) )
g(x, y) = 0.0
exact_u(x, y) = 16 * x * (Lx-x) * y * (Ly-y) + g(x, y)
