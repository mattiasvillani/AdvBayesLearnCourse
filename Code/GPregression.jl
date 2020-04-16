# Gaussian Process Regression in Julia using the GaussianProcesses.jl package

#Pkg.add("GaussianProcesses")
#Pkg.add("Random")
using GaussianProcesses, Random
using Plots

#%% Simulate some nonlinear regression data
Random.seed!(13579)
n = 10
x = 2π * rand(n)
y = sin.(x) + 0.05*randn(n)

#%% Set up the prior
mZero = MeanZero()
kern = SE(0.0,0.0)
logObsNoise = -1.0

#%% Fit GP
gp = GP(x, y, mZero, kern, logObsNoise) # Fits the GP
gp.kernel # Kernel with initial hyperparameters
gp.mll # log marginal likelihood at the initial hyperparameters
p1 = plot(gp)

#%% Prediction
x = 0:0.1:2π
μ, Σ = predict_y(gp,x, full_cov = false)

#%% Optimize hyperparameters by maximizing the log marginal likelihood
optimize!(gp)
gp.kernel # Kernel with optimized hyperparameters
gp.mll # log marginal likelihood at the optimal hyperparameters
p2 = plot(gp)
plot(p1,p2,fmt=:png)

savefig("PriorPosterior.png")
