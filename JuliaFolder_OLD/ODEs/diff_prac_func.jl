function g(x,a)
  # This is a function to calculate the digits of the log function
  value = 0.0
  for i=0:3
    value += i*a[i+1]*x^(i-1)
  end
  return value
end
