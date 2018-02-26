# This is a file to create csv files containing simple data

using PyPlot
using DataFrames

df = DataFrame()

xaxis = 2*2*π*(rand(1000)-0.5)
yaxis = sin.(xaxis).+0.3*randn(1000)

df[:x] = xaxis
df[:y] = yaxis

println(df)

plot(xaxis,yaxis,".")

using CSV
CSV.write("prac.csv",df)
