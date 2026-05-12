using ForwardDiff, Optim

f(x) = sum(x.^2) + (x[1] - 1)*x[2]

x0 = [1.0, 2.0, 3.0]

# Value, gradient, Hessian at x0
f(x0)
ForwardDiff.gradient(f, x0)
ForwardDiff.hessian(f, x0)

# Optimize with BFGS (Optim.jl uses AD automatically if you use autodiff=:forward)
x0 = [1.0, 2.0, 0.0]
optRes = optimize(f, x0, BFGS(), autodiff = :forward)

optRes.minimizer
ForwardDiff.hessian(f, optRes.minimizer)