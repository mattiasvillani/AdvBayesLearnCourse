## Bayesian analysis of the steady-state AR model with one lag

## Loading modules and settings
using Plots, Statistics, Distributions
using LaTeXStrings
using DelimitedFiles
using ColorSchemes
import ColorSchemes: Paired_12; colors = Paired_12
basicColors = Paired_12[[1,2,3,4,7,8]]
gr(legend = nothing)

## Setting up functions
"""
Simulate a AR(1) with steady-state
"""
function SimSteadyStateAR(nTime, μ, ϕ, σ, nLags = 1)
    if nLags != 1 error("Only NLags = 1 is implemented. For teaching") end
    y = zeros(nTime+1)
    for t in 1:nTime
        y[t+1] = μ + ϕ*(y[t]-μ) + σ*randn()
    end
    return(y[2:end])
end

ScaledInverseChiSq(ν,τ²) = InverseGamma(ν/2,ν*τ²/2)

"""
Gibbs sampling for the steady state AR(1) model

Samples `nIter` draws from the posterior of the model:
```math
yₜ = μ + ϕ(yₜ₋₁ - μ) + ε, where ε ∼ N(0,σ²)
```
with prior
```math
μ ∼ N(0,ψ)
ϕ | σ² ∼ N(0,σ²/κ₀)
σ² ∼ Inv-χ²(ν₀,σ₀^2)
```

"""
function GibbsSteadyStateAR(y, nIter, ψ, κ₀, ν₀, σ₀)

    # Set up storage and prelims
    n = length(y)
    postDraws = zeros(nIter,3)

    # Initial values
    x = [0;y[1:end-1]] # Set up the lag
    ϕ = x'x \ x'y
    μ = mean(y)
    σ² = var(y - ϕ*x)

    for i in 1:nIter

        # Tranforming the data a bit
        yTilde = y .- μ
        xTilde = [0;yTilde[1:end-1]] # Set up the lag

        # Draw ϕ
        ϕTilde = (xTilde'xTilde + κ₀) \ (x'yTilde)
        ϕ = rand(Normal(ϕTilde,√(σ²/(xTilde'xTilde + κ₀))))
        postDraws[i,1] = ϕ

        # Draw σ²
        νₙ = ν₀ + n
        τ² = ( (yTilde - ϕ*xTilde)'*(yTilde - ϕ*xTilde) + κ₀*ϕ^2 + ν₀*σ₀^2)/νₙ
        σ² = rand(ScaledInverseChiSq(νₙ,τ²))
        postDraws[i,2] = σ²

        # Draw μ
        yCheck = (y - ϕ*x)/(1 - ϕ)
        w = ((n/σ²)*(1-ϕ)^2)/((n/σ²)*(1-ϕ)^2 + 1/ψ^2)
        τₙ = √( 1/( (n/σ²)*(1-ϕ)^2 + 1/ψ^2 ) )
        μ = rand(Normal(w*mean(yCheck),τₙ))
        postDraws[i,3] = μ

    end

    return(postDraws)

end

## Simulating data and sampling from posterior

# Simulate data
μ = 1; ϕ = 0.5; σ = 1; n = 200;
y = SimSteadyStateAR(n, μ, ϕ, σ)
plot(y, lw = 1.5, xlab = "time", ylab = "y", color = basicColors[6]);
writedlm("timeseries.csv",  y, ' ')
z = readdlm("/home/mv/Dropbox/Teaching/AdvBayesLearnCourse/Labs/timeseries.csv")
# Set prior
ψ = 1; κ₀ = 8; ν₀ = 4; σ₀ = 1;

# Sample from posterior
nIter = 10000
postDraws = GibbsSteadyStateAR(y, nIter, ψ, κ₀, ν₀, σ₀)
mean(postDraws, dims = 1)
p1 = histogram(postDraws[:,1]);
p2 = histogram(postDraws[:,2]);
p3 = histogram(postDraws[:,3]);
plot(p1,p2,p3)
