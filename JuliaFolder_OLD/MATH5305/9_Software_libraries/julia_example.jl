using MyBlas

x = [1.0; 0.0; -3.0; 2.0]
y = [-2.0; 1.0; 5.0; 0.0]
a = 1.5

function print_vector(v, label)
    print(label)
    for j = 1:length(v)
        @printf("%8.4f", v[j])
    end
    println()
end

println("a = $a")
print_vector(x, "x = ")
print_vector(y, "y = ")

println("Julia:")
y1 = copy(y)
daxpy!(length(x), a, x, 1, y1, 1)
print_vector(y1, "ax + y = ")

println("Fortran:")
y2 = copy(y)
n = Int32(length(x))
incx = Int32(1)
incy = Int32(1)
fortran_daxpy!(n, a, x, incx, y2, incy)
print_vector(y2, "ax + y = ")

println("C:")
y3 = copy(y)
c_daxpy!(n, a, x, incx, y3, incy)
print_vector(y3, "ax + y = ")

println("Fortran impersonating C:")
y4 = copy(y)
impersonate_c_daxpy!(n, a, x, incx, y4, incy)
print_vector(y4, "ax + y = ")
