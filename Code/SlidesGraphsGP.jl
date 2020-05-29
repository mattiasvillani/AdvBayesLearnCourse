using Plots, Statistics, GaussianProcesses
using LaTeXStrings
using ColorSchemes
import ColorSchemes: Paired_12; colors = Paired_12
basicColors = Paired_12[[1,2,3,4,7,8]]

default(legend = true)
figFolder = "/home/mv/Dropbox/Teaching/AdvBayesLearnCourse/Slides/Figures/"


## Sum of cosines
xGrid = -π:0.01:π
plot()
sumCos = zeros(size(xGrid))
sGrid = [0.10,0.25,0.5,1.0]
for (i,s) in enumerate(sGrid)
    a = ceil(0.5*rand(); digits = 1)
    global sumCos += a*cos.(2π*s*xGrid)
    plot!(xGrid, a*cos.(2π*s*xGrid); width = 1.5, color = basicColors[i],
        label = latexstring("$a \\pi\\cos($s x)"))
end
plot!(xGrid, sumCos; ylabel = latexstring("x(t)"), xlabel = latexstring("t"),
    width = 2, label = "sum of cosines", color = basicColors[6])
current()

savefig(figFolder*"SumCosines.svg")

## Simulate data from sum of cosines
plot(xGrid, sumCos + 0.2*randn(length(sumCos)); ylabel = latexstring("x(t)"), xlabel = latexstring("t"),
    width = 2, legend = nothing, color = basicColors[6])
current()

savefig(figFolder*"SumCosinesData.svg")

## Plotting ellipsoids
μ = [1,2,0.5]
Λ = Diagonal([0.01,0.5,3])
V = [-1 -1 2;1 -1 0;1 1 1]'
V = V*(V'*V)^-(1/2)
Σ = V*Λ*V'
MvNormal(μ, Σ)
X = rand(MvNormal(μ, Σ),1000)'
p1 = plot3d(label = "data", legend = nothing)
scatter3d!(X[:,1],X[:,2],X[:,3], xlab = latexstring("x_1"),
    ylab = latexstring("x_2"), zlab = latexstring("x_3"),
    markerstrokewidth = 0, color = basicColors[6],
    title = "One dominant direction", guidefont=font(8), tickfont=font(6))

Λ = Diagonal([0.01,0.5,3])
V = Diagonal(ones(3))
Σ = V*Λ*V'
MvNormal(μ, Σ)
X = rand(MvNormal(μ, Σ),1000)'
p2 = plot3d(label = "data", legend = nothing)
scatter3d!(X[:,1],X[:,2],X[:,3],
    xlab = latexstring("x_1"), ylab = latexstring("x_2"), zlab = latexstring("x_3"),
    markerstrokewidth = 0, color = basicColors[6],
    title = "No dominant direction", guidefont=font(8), tickfont=font(6))

l = @layout [a b]
plot(p2, p1, layout = l, size = (600, 250))

png(figFolder*"PCscatter")

Σ
A = ones(4,4)
eigen(A).values

μ = [1,2,0.5]
Λ = Diagonal([0.0,0.5,3])
V = [-1 -1 2;1 -1 0;1 1 1]'
V = V*(V'*V)^-(1/2)
Σ = V*Λ*V'
eigen(Σ)
cholesky(Σ)
MvNormal(μ, Σ)
X = rand(MvNormal(μ, Σ),10000)
cov(X')

eigen(zeros(3,3))

gr()
n = 100 # Number of training data points
x = rand(n)*3
trueMean = x -> sin.(2*x)
trueVar = x -> 0.1^2*exp.(0*x)
σ2 = 0.01
y = trueMean(x) + randn(n).*sqrt.(σ2);
xGrid = minimum(x):0.01:maximum(x);
p = plot(xlab = "x", ylims = [-2,1.2])
plot!(p, xGrid, trueMean(xGrid), lw = 2, color = basicColors[6],
    label = "True mean")
scatter!(p, x, y; ms = 4, xlab = "x", ylab = "y", markerstrokewidth = 0,
    label = "data", color = basicColors[2])

# Radial basis function
xGrid = 0:0.01:3
cGrid = 0:0.25:3
ℓ = 0.25
τ2 = 1/length(cGrid)
X = zeros(n,length(cGrid))
for (i,c) in enumerate(cGrid)
    X[:,i] = exp.(-(1/(2*ℓ^2))*(x.-c).^2)
end
w = inv(X'*X + (σ2/τ2)*Diagonal(ones(length(cGrid))))*X'*y
w = X \ y
Xgrid = zeros(length(xGrid),length(cGrid))
for (i,c) in enumerate(cGrid)
    Xgrid[:,i] = exp.(-(1/(2*ℓ^2))*(xGrid.-c).^2)
end
plot!(p, xGrid, Xgrid*w, lw = 2, label = "radial basis fit")


for (i,c) in enumerate(cGrid)
    plot!(p, xGrid,-2 .+ 1*exp.(-(1/(2*ℓ^2))*(xGrid.-c).^2), lw = 2, label = nothing)
end
p

png(figFolder*"SEInfiniteSplines")
