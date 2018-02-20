# This is the solution to Question 6
# We will use a string to determine which intial values to use


 #initial = "random"
 initial = "data"

include("AllenCahn.jl")
using AllenCahn
using PyPlot
using OffsetArrays

const Ω = (0.0,2.0)
const ϵ = 0.015
const c = 3 / (sqrt(2)*ϵ)
const T = 0.02
const P = 256
const N = 2500

gr = grid1D(Ω, T, P, N)

x, t = gr.x, gr.t

function h(x)
  if (0.2<x<0.4)
    h = 0.3
  elseif (0.8<x<1.1)
    h = -0.8
  elseif (1.4<x<1.5)
    h = 0.5
  else
    h = 0.0
  end
  return h
end


if initial=="data"
  u0 = h.(x)
elseif initial=="random"
  u0 = 2*(rand(P) - 0.5)*10.0^(-3)
else
  println("There is an error with which initial values to use.")
end

U = semi_implicit(ϵ, gr, u0)

# Find the indexes of t which correspond to the times given in the question
if initial=="data"
  t0(x) = x == 0.0
  t1(x) = x == 0.00012
  t2(x) = x == 0.00024
  t3(x) = x == 0.00036
  t4(x) = x == 0.00048

  indt0 = find(t0,t)
  indt1 = find(t1,t)
  indt2 = find(t2,t)
  indt3 = find(t3,t)
  indt4 = find(t4,t)
elseif initial=="random"
  t0(x) = x == 0.0
  t1(x) = x == 0.00024
  t2(x) = x == 0.00048
  t3(x) = x == 0.00072
  t4(x) = x == 0.00096

  indt0 = find(t0,t)
  indt1 = find(t1,t)
  indt2 = find(t2,t)
  indt3 = find(t3,t)
  indt4 = find(t4,t)
end


idx = [indt0; indt1; indt2; indt3; indt4]
#idx = 0:5:100
figure(1)
plot(x[1:P],U[1:P,idx],linewidth=1)
xlabel(L"$x$")
ylabel(L"$u$")
grid(false)
ax = gca()
ax[:set_xlim]([0,2])
if initial == "data"
  ax[:set_ylim]([-1,1])
  legend((L"$t=0.0$", L"$t=0.00012$", L"$t=0.00024$",
                  L"$t=0.00036$", L"$t=0.00048$"),loc = "lower right")
elseif initial == "random"
  ax[:set_ylim]([minimum(U[:,idx]),maximum(U[:,idx])])
  legend((L"$t=0.0$", L"$t=0.00024$", L"$t=0.00048$",
                  L"$t=0.00072$", L"$t=0.00096$"))
end
