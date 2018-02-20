function f(x,a)
  # This is a function to calculate the digits of the log function
  value = 0.0
  for i=0:4
    value += a[i+1]*x.^i
  end
  #return 1.0+6.0*x-7.0*x.^2+4.0*x.^3+3.0*x.^4
  return value
end
