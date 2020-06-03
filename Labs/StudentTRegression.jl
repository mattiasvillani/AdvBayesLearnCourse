# Gibbs sampling for the regression model with student-t errors
# Author: Mattias Villani, for Lab 4 of the course Advanced Bayesian Learning

using Plots
using Distributions
using LaTeXStrings
using DelimitedFiles
import ColorSchemes: Paired_12; colors = Paired_12
basicColors = Paired_12[[1,2,3,4,7,8]]

# Defining some distributions
ScaledInverseChiSq(ν,τ²) = InverseGamma(ν/2,ν*τ²/2)
TDistGeneral(μ,σ,ν) = LocationScale(μ,σ,TDist(ν))

"""
    Gibbs sampling for the regression model with student-t errors
        y = x'*β + ε, ε iid ~ t(0,σ,μ)
    with prior
        β | σ ~ N(0,τσ*I), log(σ) ∼ Uniform
"""
function GibbsTReg(y, X, ν, τ, nIter)
    n, q = size(X)
    β = X \ y
    ε = y - X*β
    σ = norm(ε)/√(n-q)
    postDraws = zeros(nIter,q+1)
    for i in 1:nIter

        # Update U
        U = [rand(ScaledInverseChiSq(ν+1,(ν*σ^2 + εi^2)/(ν+1))) for εi in ε]

        # Update β
        A = Diagonal(ones(n)./(.√U))
        yt = A*y
        Xt = A*X
        βhat = Xt \ yt # Least squares estimate
        Σpost = inv(Xt'*Xt + τ^(-2)*I(q))
        μpost = vec(Σpost*Xt'*Xt*βhat)
        β = rand(MvNormal(μpost,Σpost))
        postDraws[i,1:q] = β
        ε = y - X*β

        # Update σ
        σ = √rand(Gamma(n*ν/2, 1/((ν/2)*sum(ones(n)./U))))
        postDraws[i,q+1] = σ

    end # end Gibbs iterations
    return postDraws
end # end function GibbsStudentTReg


"""
    Predictive distribution for the regression model with student-t errors
        y = x'*β + ε, ε iid ~ t(0,σ,μ)
    with prior
        β | σ ~ N(0,τσ*I), log(σ) ∼ Uniform
"""
function PredTReg(X, postDraws, ν)
    nPreds,q = size(X)
    nSim = size(postDraws,1)
    yPreds = zeros(nSim,nPreds)
    for i = 1:nSim
        β = postDraws[i,1:q]
        σ = postDraws[i,end]
        yPreds[i,:] = [rand(TDistGeneral(μ, σ, ν)) for μ in X*β]
    end
    return yPreds
end

# Read data from file
data = readdlm("./Labs/regression.csv")
X = [ones(size(data,1),1) data[:,1]]
y = data[:,2]
scatter(X[:,2],y)

# Simulating from the posterior of β and σ
ν = 5
τ = 1
nIter = 10000
postDraws = GibbsTReg(y, X, ν, τ, nIter)
mean(postDraws, dims = 1)
p1 = histogram(postDraws[:,1], nbins = 100, title = L"\beta_0");
p2 = histogram(postDraws[:,2], nbins = 100, title = L"\beta_1");
p3 = histogram(postDraws[:,3].^2*(ν/(ν-2)), nbins = 100, title = L"\mathbb{V}(y)=\sigma^2\frac{\nu}{\nu-2}");
plot(p1,p2,p3)

# Simulating from the predictive distribution
xGrid = -1:0.1:1
XPred = [ones(length(xGrid),1) xGrid]
yPreds = PredTReg(XPred, postDraws, ν)
predBands = zeros(size(XPred,1),5)
for i in 1:size(yPreds,2)
    predBands[i,:] = quantile(yPreds[:,i], [0.005 0.025 0.5 0.975 0.995])
end
p = scatter(X[:,2],y, xlabel = L"x", ylabel = L"y");
plot!(xGrid,predBands[:,3], color = basicColors[6])
plot!(xGrid,predBands[:,1], color = basicColors[1])
plot!(xGrid,predBands[:,2], color = basicColors[2])
plot!(xGrid,predBands[:,4], color = basicColors[2])
plot!(xGrid,predBands[:,5], color = basicColors[1])
