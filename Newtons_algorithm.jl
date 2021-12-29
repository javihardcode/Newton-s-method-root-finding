##########################################
# Newtons' Algorithm 
##########################################

using ForwardDiff 
using Plots
D(f) = w -> ForwardDiff.derivative(f, w)

w = 0.5:0.01:2.5

# Parameters
lamd = 1/3
alph = 1/2
gamm = 3/4

# Draw Supply & Demand
l1(w) = (lamd/w).^(1/(1-lamd))
l2(w) = (alph./w).^(1/(1-alph))
l3(w) = (gamm./w).^(1/(1-gamm))
L(w)  = l1(w) + l2(w) + l3(w) 

plot(w, l1, label="l1")
plot!(w,l2, label="l2")
plot!(w,l3, label="l3",xlabel="wage",title="Supply & Demand")
plot!(w,L, label="Aggregate Demand")

s(w) = 3/2
plot!(w,s, label="Aggregate Supply")

# Now define the functions for the root finding problem
# f(w)-3/2 is the market clearing condition for w
f(w) = L(w)-3/2
L_prime = D(f)     # f´(x)

maxiter = 10000
iter = 1
normdiff = 1
w0 = 0.1 
w_old = w0 
tolerance = 1.0E-10

while normdiff > tolerance && iter <= maxiter
     w_new = w_old - f(w_old) / L_prime(w_old)
     normdiff = abs(w_new - w_old)
     w_old = w_new 
     iter = iter+1
end 

print(" |l(w)-x|= $normdiff in $iter iterations, solution $w_old" )

